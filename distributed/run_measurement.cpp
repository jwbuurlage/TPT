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
#include "tomos/util/read_tiff.hpp"
#include "tomos/util/simple_args.hpp"
#include "tomos/util/trees.hpp"

#include "overlaps.hpp"

using namespace tomo;

using T = float;

// Assumes global volume has physical size [0, 1]^3
auto calculate_local_volume(bulk::rectangular_partitioning<3, 1>& partitioning,
                            int s) {

    auto to_vec = [](const auto& array_like) -> math::vec3<int> {
        math::vec3<int> result;
        for (int d = 0; d < 3; ++d) {
            result[d] = array_like[d];
        }
        return result;
    };

    auto voxel_origin = to_vec(partitioning.origin(s));
    auto voxel_size = to_vec(partitioning.local_size(s));
    auto global_voxels = to_vec(partitioning.global_size());

    auto relative_origin =
        math::vec3<T>(voxel_origin) / math::vec3<T>(global_voxels);
    auto relative_size =
        math::vec3<T>(voxel_size) / math::vec3<T>(global_voxels);

    return volume<3_D, T>(voxel_size, relative_origin, relative_size);
}

struct shared_pixel {
    int target;
    int local_line;
    int remote_line;
};

struct pixel_message {
    int local_line;
    int remote_line;
};

std::pair<std::vector<std::vector<pixel_message>>,
          std::vector<std::vector<pixel_message>>>
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
                        exit(-1);
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

    auto a_targets = std::vector<std::vector<int>>(p);
    auto a_locals = std::vector<std::vector<int>>(p);
    auto a_remotes = std::vector<std::vector<int>>(p);
    auto b_targets = std::vector<std::vector<int>>(p);
    auto b_locals = std::vector<std::vector<int>>(p);
    auto b_remotes = std::vector<std::vector<int>>(p);

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

                a_targets[t].push_back(owner);
                a_locals[t].push_back(contributor_line);
                a_remotes[t].push_back(owner_line);

                b_targets[owner].push_back(t);
                b_locals[owner].push_back(owner_line);
                b_remotes[owner].push_back(contributor_line);
            }
        }
    }

    pixels.clear();
    pixels.shrink_to_fit();

    auto my_contributions = bulk::queue<int[], int[], int[]>(world);
    auto my_responsibilities = bulk::queue<int[], int[], int[]>(world);
    for (int t = 0; t < p; ++t) {
        my_contributions(t).send(std::move(a_targets[t]), std::move(a_locals[t]), std::move(a_remotes[t]));
        my_responsibilities(t).send(std::move(b_targets[t]), std::move(b_locals[t]), std::move(b_remotes[t]));
    }
    world.sync();

    // 5. store contributors and owners in a vector
    auto contributions = std::vector<std::vector<pixel_message>>(p);
    auto responsibilities = std::vector<std::vector<pixel_message>>(p);

    for (auto& xs : my_contributions) {
        auto& owners = std::get<0>(xs);
        auto& locals = std::get<1>(xs);
        auto& remotes = std::get<2>(xs);
        for (auto i = 0u; i < owners.size(); ++i) {
            contributions[owners[i]].push_back({locals[i], remotes[i]});
        }
    }
    for (auto xs : my_responsibilities) {
        auto& owners = std::get<0>(xs);
        auto& locals = std::get<1>(xs);
        auto& remotes = std::get<2>(xs);
        for (auto i = 0u; i < owners.size(); ++i) {
            responsibilities[owners[i]].push_back({locals[i], remotes[i]});
        }
    }

    return {std::move(contributions), std::move(responsibilities)};
}

void communicate_contributions(
    bulk::world& world, projections<3_D, T>& projs,
    const std::vector<std::vector<pixel_message>>& contributions,
    const std::vector<std::vector<pixel_message>>& results) {
    std::vector<int> remotes_buf;
    std::vector<T> values_buf;
    int p = world.active_processors();
    auto q = bulk::queue<int[], T[]>(world);
    auto share = [&, p](const auto& xs) {
        for (int t = 0; t < p; ++t) {
            remotes_buf.clear();
            values_buf.clear();
            for (auto result : xs[t]) {
                auto local = result.local_line;
                auto remote = result.remote_line;
                remotes_buf.push_back(remote);
                values_buf.push_back(projs[local]);
            }
            q(t).send(remotes_buf, values_buf);
        }
        world.sync();
    };

    share(contributions);
    {
        auto& work = *q.begin();
        auto& remotes = std::get<0>(work);
        auto& values = std::get<1>(work);
        for (auto i = 0u; i < remotes.size(); ++i) {
            projs[remotes[i]] += values[i];
        }
    }
    share(results);
    {
        auto& work = *q.begin();
        auto& remotes = std::get<0>(work);
        auto& values = std::get<1>(work);
        for (auto i = 0u; i < remotes.size(); ++i) {
            projs[remotes[i]] = values[i];
        }
    }
}

void collect_orthos(bulk::world& world, tomo::image<3_D, T> x,
                    tomo::volume<3_D, T> global_volume,
                    bulk::rectangular_partitioning<3, 1>& part,
                    std::string name, std::string output_dir) {
    auto s = world.rank();

    auto voxels = global_volume.voxels();
    auto center = voxels / 2;
    auto xy = bulk::queue<int, int, T>(world);
    auto yz = bulk::queue<int, int, T>(world);
    auto xz = bulk::queue<int, int, T>(world);

    if (center[0] >= part.origin(s)[0] &&
        part.local_size(s)[0] > center[0] - part.origin(s)[0]) {
        for (int i = 0; i < part.local_size(s)[1]; ++i) {
            for (int j = 0; j < part.local_size(s)[2]; ++j) {
                yz(0).send(i + part.origin(s)[1], j + part.origin(s)[2],
                           x[x.index({center[0] - part.origin(s)[0], i, j})]);
            }
        }
    }

    if (center[1] >= part.origin(s)[1] &&
        part.local_size(s)[1] > center[1] - part.origin(s)[1]) {
        for (int i = 0; i < part.local_size(s)[0]; ++i) {
            for (int j = 0; j < part.local_size(s)[2]; ++j) {
                xz(0).send(i + part.origin(s)[0], j + part.origin(s)[2],
                           x[x.index({i, center[1] - part.origin(s)[1], j})]);
            }
        }
    }

    if (center[2] >= part.origin(s)[2] &&
        part.local_size(s)[2] > center[2] - part.origin(s)[2]) {
        for (int i = 0; i < part.local_size(s)[0]; ++i) {
            for (int j = 0; j < part.local_size(s)[1]; ++j) {
                xy(0).send(i + part.origin(s)[0], j + part.origin(s)[1],
                           x[x.index({i, j, center[2] - part.origin(s)[2]})]);
            }
        }
    }

    world.sync();

    if (s == 0) {
        auto slice_xy =
            tomo::image<2_D, T>(tomo::volume<2_D, T>({voxels[0], voxels[1]}));
        auto slice_xz =
            tomo::image<2_D, T>(tomo::volume<2_D, T>({voxels[0], voxels[2]}));
        auto slice_yz =
            tomo::image<2_D, T>(tomo::volume<2_D, T>({voxels[1], voxels[2]}));

        for (auto result : xy) {
            auto i = std::get<0>(result);
            auto j = std::get<1>(result);
            auto value = std::get<2>(result);
            slice_xy[slice_xy.index({i, j})] = value;
        }

        for (auto result : xz) {
            auto i = std::get<0>(result);
            auto j = std::get<1>(result);
            auto value = std::get<2>(result);
            slice_xz[slice_xz.index({i, j})] = value;
        }

        for (auto result : yz) {
            auto i = std::get<0>(result);
            auto j = std::get<1>(result);
            auto value = std::get<2>(result);
            slice_yz[slice_yz.index({i, j})] = value;
        }

        tomo::write_png(slice_xy, output_dir + name + "_slice_xy");
        tomo::write_png(slice_xz, output_dir + name + "_slice_xz");
        tomo::write_png(slice_yz, output_dir + name + "_slice_yz");
    }
}

void sirt(bulk::world& world,
          bulk::rectangular_partitioning<3, 1>& partitioning,
          tomo::volume<3_D, T> global_volume,
          geometry::trajectory<3_D, T>& global_geometry,
          tomo::util::report& table, std::string name, std::string column,
          std::string image_dir, int iters) {

    if (world.rank() == 0) {
        world.log("Running %s (%s)", name.c_str(), column.c_str());
    }

    auto vs = calculate_local_volume(partitioning, world.rank());
    image<3_D, T> phantom(vs);
    fill_ellipsoids_(phantom, mshl_ellipsoids_<T>(), vs, global_volume);

    auto gs = distributed::restricted_geometry<T>(global_geometry, vs);
    auto p = projections<3_D, T>(gs);

    using dimmer = dim::closest<3_D, T>;
    auto fp = [&](const auto& f, const auto& g, const auto& v, auto& q) {
        auto proj = dimmer(v);
        int line_number = 0;
        for (auto line : g) {
            for (auto elem : proj(line)) {
                q[line_number] += f[elem.index] * elem.value;
            }
            ++line_number;
        }
    };

    if (world.rank() == 0) {
        world.log("Running experiment %s (%s)", name.c_str(), column.c_str());
    }
    fp(phantom, gs, vs, p);

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
    if (world.rank() == 0) {
        world.log("Computing contributions %s (%s)", name.c_str(),
                  column.c_str());
    }
    auto result = compute_contributions(world, gs);
    auto& go_forth = result.first;
    auto& and_back = result.second;

    // communicate row and column sums
    if (world.rank() == 0) {
        world.log("Compute R and C %s (%s)", name.c_str(), column.c_str());
    }
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

    if (world.rank() == 0) {
        world.log("SIRT %s (%s)", name.c_str(), column.c_str());
    }
    auto stopwatch = bulk::util::timer();
    for (int iter = 0; iter < iters; ++iter) {
        fp(x, gs, vs, q);
        for (auto l = 0ull; l < gs.lines(); ++l) {
            q[l] = r[l] * (p[l] - q[l]);
        }
        communicate_contributions(world, q, go_forth, and_back);
        bp(q, gs, vs, z);
        for (auto i = 0ull; i < vs.cells(); ++i) {
            x[i] += c[i] * z[i];
        }

        z.clear();
        q.clear();
    }

    table.add_result(
        name, column,
        fmt::format("{:.2f}", stopwatch.get<std::milli>() / (1000.0 * iters)));

    collect_orthos(world, x, global_volume, partitioning, name + "_" + column,
                   image_dir);
}

void run(const std::vector<std::string>& geoms, std::string part_dir, int k,
         int iters, std::string outfile, std::string image_dir,
         tomo::util::report& table, bool trivial, bool bisected) {
    bulk::mpi::environment env;

    auto processors = env.available_processors();

    env.spawn(processors, [&](auto& world) {
        auto p = world.active_processors();

        for (auto geom_file : geoms) {
            auto name = fs::path(geom_file).stem().string();
            table.add_row(name);

            // this creates a volume [0, 1]^3, with k voxels in each axis
            auto problem = tomo::read_configuration<3_D, T>(geom_file, k);
            auto& global_geometry = *problem.acquisition_geometry;
            auto global_volume = problem.object_volume;

            auto tree_file = part_dir + "/" + std::to_string(processors) + "/" +
                             name + ".bsp";

            auto tree_partitioning =
                tomo::load_partitioning(tree_file, global_volume, log2(p));

            auto main_d = tree_partitioning->splits().root->value.d;

            // which dimension, the first split decides.....
            auto block_partitioning = bulk::block_partitioning<3, 1>(
                tomo::math::vec_to_array<3_D, int>(global_volume.voxels()), {p},
                {main_d});

            if (trivial) {
		    sirt(world, block_partitioning, global_volume,
			 (tomo::geometry::trajectory<3_D, T>&)global_geometry, table,
			 name, "trivial", image_dir, iters);
		}
if (bisected) {
            sirt(world, *tree_partitioning, global_volume,
                 (tomo::geometry::trajectory<3_D, T>&)global_geometry, table,
                 name, "bisected", image_dir, iters);
}
        }

        if (world.rank() == 0) {
            table.print();
            std::ofstream of(outfile, std::ios::out);
            table.print(of);
        }
    });
}

void usage(std::string program_name) {
    std::cout << "Usage: " << program_name
              << " --geom GEOMS --part PART_DIR "
                 "--out TABLE_FILE --images IMAGE_DIR [-k SIZE] [-i "
                 "ITERS]\n";
}

int main(int argc, char* argv[]) {
    auto opts = options{argc, argv};

    if (!opts.required_arguments({"--geom", "--part", "--out", "--images"})) {
        usage(argv[0]);
        return -1;
    }

    auto table = tomo::util::report("Runtimes", "geometry");
    table.add_column("trivial");
    table.add_column("bisected");

    run(opts.args("--geom"), opts.arg("--part"), opts.arg_as_or<int>("-k", -1),
        opts.arg_as_or<int>("-i", 1), opts.arg("--out"), opts.arg("--images"),
        table, opts.passed("--trivial"), opts.passed("--bisected"));

    return 0;
}
