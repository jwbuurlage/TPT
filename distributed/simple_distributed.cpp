/*
 * FIXME:
 * - Can we group sends together? this should greatly reduce data size
 *   with a constant factor.
 */

#include <cassert>
#include <sstream>

#include "bulk/backends/mpi/mpi.hpp"
#include "bulk/bulk.hpp"

namespace bulk {
using namespace experimental;
}

#include "tomos/tomos.hpp"
#include "tomos/util/plotter.hpp"

#include "overlaps.hpp"

using namespace tomo;

using T = float;

// Assumes global volume has physical size [0, 1]^3
template <int G>
auto calculate_local_volume(bulk::rectangular_partitioning<3, G>& partitioning,
                            int s) {
    auto to_vec = [](auto array_like) -> math::vec3<T> {
        math::vec3<T> result;
        for (int d = 0; d < 3; ++d) {
            result[d] = array_like[d];
        }
        return result;
    };

    auto voxel_origin = to_vec(partitioning.origin(s));
    auto voxel_size = to_vec(partitioning.local_size(s));
    auto global_voxels = to_vec(partitioning.global_size());

    auto relative_origin = math::vec3<T>(voxel_origin / global_voxels);
    auto relative_size = math::vec3<T>(voxel_size / global_voxels);

    return volume<3_D, T>(voxel_size, relative_origin, relative_size);
}

struct shared_pixel {
    int target;
    int local_line;
    int remote_line;
};

std::pair<std::vector<shared_pixel>, std::vector<shared_pixel>>
compute_contributions(bulk::world& world,
                      const distributed::restricted_geometry<T>& geometry) {
    auto s = world.rank();
    auto p = world.active_processors();
    auto& gg = geometry.global_geometry();

    // block_size is number of local projections
    auto block_size = ((geometry.projection_count() - 1) / p) + 1;
    auto shadows = bulk::coarray<pod_shadow>(world, p * block_size);
    auto offsets = bulk::coarray<int>(world, p * block_size);

    // 1. send all `p` windows for all projections `i` to processor `i % p`
    for (int i = 0; i < geometry.projection_count(); ++i) {
        auto sh = geometry.local_shadow(i);
        shadows(i % p)[(i / p) * p + s] = to_pod(sh);
        offsets(i % p)[(i / p) * p + s] = geometry.offset(i);
    }
    world.sync();

    // 2. for each projection, for each pixel store the contributors
    auto global_shape = gg.projection_shape(0);
    auto pixel_count = math::product<2_D, int>(global_shape);
    auto pixels = std::vector<std::vector<int>>(pixel_count * block_size);
    for (int b = 0; b < block_size; ++b) {
        for (int t = 0; t < p; ++t) {
            auto w = shadows[b * p + t];
            for (int j = w.min_pt[1]; j <= w.max_pt[1]; ++j) {
                for (int i = w.min_pt[0]; i <= w.max_pt[0]; ++i) {
                    if (j * global_shape[0] + i >= pixel_count) {
                        std::cout << b << " [" << w.min_pt[0] << ", "
                                  << w.max_pt[0] << "] -> " << i << " " << j
                                  << " / " << pixel_count << " = "
                                  << global_shape[0] << " x " << global_shape[1]
                                  << "\n";
                        assert(false);
                    }
                    pixels[b * pixel_count + j * global_shape[0] + i].push_back(
                        t);
                }
            }
        }
    }

    // 3. decide who will be responsible for the pixel
    //    maybe (line_number % #contributors)
    auto local_idx = [](auto offset, auto lw, auto gw, auto idx) {
        auto x = idx % gw[0] - lw.min_pt[0];
        auto y = idx / gw[0] - lw.min_pt[1];
        auto idx_in_proj = y * (lw.max_pt[0] - lw.min_pt[0] + 1) + x;
        return offset + idx_in_proj;
    };

    auto my_contributions = bulk::queue<shared_pixel>(world);
    auto my_responsibilities = bulk::queue<shared_pixel>(world);

    for (int b = 0; b < block_size; ++b) {
        for (int i = 0; i < pixel_count; ++i) {
            auto pixel_idx = b * pixel_count + i;
            if (pixels[pixel_idx].size() == 0) {
                continue;
            }
            auto owner =
                pixels[pixel_idx][pixel_idx % pixels[pixel_idx].size()];
            // for contributors, we can compute the local index if we know the
            // offset.. lets just send them... aaaaand got 'em..

            auto owner_line =
                local_idx(offsets[b * p + owner], shadows[b * p + owner],
                          gg.projection_shape(0), i);

            for (auto j = 0u; j < pixels[pixel_idx].size(); ++j) {
                auto t = pixels[pixel_idx][j];
                if (t == owner) {
                    continue;
                }

                // 4. send to the owner the contributors, to the contributors
                // the owner
                auto contributor_line =
                    local_idx(offsets[b * p + t], shadows[b * p + t],
                              gg.projection_shape(0), i);
                my_contributions(t).send({owner, contributor_line, owner_line});
                my_responsibilities(owner).send(
                    {t, owner_line, contributor_line});
            }
        }
    }

    world.sync();

    // 5. store contributors and owners in a vector
    auto contributions = std::vector<shared_pixel>(my_contributions.size());
    auto responsibilities =
        std::vector<shared_pixel>(my_responsibilities.size());
    std::copy(my_contributions.begin(), my_contributions.end(),
              contributions.begin());
    std::copy(my_responsibilities.begin(), my_responsibilities.end(),
              responsibilities.begin());

    return {std::move(contributions), std::move(responsibilities)};
}

void communicate_contributions(bulk::world& world, projections<3_D, T>& projs,
                               const std::vector<shared_pixel>& contributions,
                               const std::vector<shared_pixel>& results) {
    auto q = bulk::queue<int, T>(world);
    auto share = [&](const auto& xs) {
        for (auto [target, local, remote] : xs) {
            q(target).send(remote, projs[local]);
        }
    };

    share(contributions);
    world.sync();
    for (auto [idx, value] : q) {
        projs[idx] += value;
    }

    share(results);
    world.sync();
    for (auto [idx, value] : q) {
        projs[idx] = value;
    }
}

void run(util::args options) {
    auto k = options.k;

    bulk::mpi::environment env;

    // this creates a volume [0, 1]^3, with k voxels in each axis
    auto global_volume = volume<3_D, T>(k);
    auto global_geometry =
        geometry::cone_beam<T>(global_volume, k, {2.0, 2.0}, {k, k}, 2.0, 2.0);

    auto processors = env.available_processors();
    auto partitioning = bulk::block_partitioning<3, 3>({k, k, k}, {2, 1, 4});

    // TODO: load partitioning result
    //    auto partitioning = bulk::block_partitioning<3, 3>(
    //        {k, k, k}, {2, 1, 4});

    env.spawn(processors, [&](auto& world) {
        auto s = world.rank();
        auto vs = calculate_local_volume(partitioning, world.rank());

        // Make scene with s = 0, connect to this scene from other procs
        util::ext_plotter<3_D, T> plotter("tcp://localhost:5555");
        auto scene_id = bulk::var<int>(world);
        if (s == 0) {
            plotter.make_scene("Distributed partial updates");
            plotter.send_geometry(global_geometry, global_volume);
            plotter.send_partition_information(
                partitioning, world.active_processors(), global_volume);
            scene_id.broadcast(plotter.scene_id());
        }
        world.sync();
        if (s != 0) {
            plotter.set_scene_id(scene_id);
        }
        plotter.group_request(world.active_processors());

        image<3_D, T> phantom(vs);
        fill_ellipsoids_(phantom, mshl_ellipsoids_<T>(), vs, global_volume);

        auto gs = distributed::restricted_geometry(global_geometry, vs);
        auto p = projections<3_D, T>(gs);

        using dimmer = dim::joseph<3_D, T>;
        auto fp = [](const auto& f, const auto& g, const auto& v, auto& q) {
            auto proj = dimmer(v);
            int line_number = 0;
            for (auto line : g) {
                for (auto elem : proj(line)) {
                    q[line_number] += f[elem.index] * elem.value;
                }
                ++line_number;
            }
        };
        fp(phantom, gs, vs, p);

        plotter.send_partial_projection_data(gs, p);

        auto bp = [](const auto& p_, const auto& g, const auto& v, auto& x_) {
            auto proj = dimmer(v);
            int line_number = 0;
            for (auto line : g) {
                for (auto elem : proj(line)) {
                    x_[elem.index] += p_[line_number] * elem.value;
                }
                ++line_number;
            }
        };

        auto [go_forth, and_back] = compute_contributions(world, gs);

        // communicate row and column sums
        auto invert = [](auto x) { return (T)1.0 / x; };
        auto invert_all = [&](auto& xs) {
            std::transform(xs.begin(), xs.end(), xs.begin(), invert);
        };
        auto r = projections<3_D, T>(gs);
        fp(image<3_D, T>(vs, (T)1.0), gs, vs, r);
        communicate_contributions(world, r, go_forth, and_back);
        invert_all(r.mutable_data());
        auto c = image<3_D, T>(vs);
        bp(projections<3_D, T>(gs, (T)1.0), gs, vs, c);
        invert_all(c.mutable_data());

        // buffer proj stack
        auto q = projections<3_D, T>(gs);
        auto z = image<3_D, T>(vs);
        auto x = image<3_D, T>(vs);
        for (int iter = 0; iter < 10; ++iter) {
            fp(x, gs, vs, q);
            for (int l = 0; l < gs.lines(); ++l) {
                q[l] = r[l] * (p[l] - q[l]);
            }
            communicate_contributions(world, q, go_forth, and_back);
            bp(q, gs, vs, z);
            for (int i = 0; i < vs.cells(); ++i) {
                x[i] += c[i] * z[i];
            }

            z.clear();
            q.clear();
        }
    });
}

int main(int argc, char* argv[]) {
    auto options = util::args(argc, argv);
    run(options);

    return 0;
}
