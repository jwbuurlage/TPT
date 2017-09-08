/*
 * FIXME:
 * - Can we group sends together? this should greatly reduce data size
 *   with a constant factor.
 * ! we do it now because memory constraints
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

    if (s == 0) {
        world.log("coarrays");
    }
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
    if (s == 0) {
        world.log("pixels");
    }
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

    if (s == 0) {
        world.log("queues");
    }
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
    pixels.clear();

    if (s == 0) {
        world.log("clear pixels");
    }

    // 5. store contributors and owners in a vector
    if (s == 0) {
        world.log("store shared pixels, total pixel count %u + %u",
                  my_contributions.size(), my_responsibilities.size());
    }
    auto contributions = std::vector<shared_pixel>(my_contributions.size());
    auto responsibilities =
        std::vector<shared_pixel>(my_responsibilities.size());
    std::copy(my_contributions.begin(), my_contributions.end(),
              contributions.begin());
    std::copy(my_responsibilities.begin(), my_responsibilities.end(),
              responsibilities.begin());

    return {std::move(contributions), std::move(responsibilities)};
}

void count_messages(bulk::world& world,
                    bulk::rectangular_partitioning<3_D, 1>& partitioning,
                    geometry::trajectory<3_D, T>& global_geometry,
                    tomo::util::report& table, std::string name,
                    std::string column) {
    if (world.rank() == 0) {
        world.log("calculating for %s", name.c_str());
    }
    auto vs = calculate_local_volume(partitioning, world.rank());
    auto gs = distributed::restricted_geometry<T>(global_geometry, vs);

    auto result = compute_contributions(world, gs);
    auto& go_forth = result.first;
    auto& and_back = result.second;

    auto targets_to = std::set<int>();
    auto targets_fro = std::set<int>();
    for (auto x : go_forth) {
        targets_to.insert(x.target);
    }
    for (auto x : and_back) {
        targets_fro.insert(x.target);
    }

    auto x = bulk::var<int>(world, targets_to.size() + targets_fro.size());
    auto messages = bulk::foldl(x, [](int& lhs, int rhs) { lhs += rhs; });

    table.add_result(name, column, messages);
}

void run(const std::vector<std::string>& geoms, std::string part_dir,
         std::string outfile, tomo::util::report& table) {
    bulk::mpi::environment env;

    auto processors = env.available_processors();

    env.spawn(processors, [&](auto& world) {
        auto p = world.active_processors();
        if (world.rank() == 0) {
            world.log("Running message counting with %i processors", p);
        }

        for (auto geom_file : geoms) {
            auto name = fs::path(geom_file).stem().string();
            table.add_row(name);

            // this creates a volume [0, 1]^3, with k voxels in each axis
            auto problem = tomo::read_configuration<3_D, T>(geom_file);
            auto& global_geometry = *problem.acquisition_geometry;
            auto global_volume = problem.object_volume;

            auto tree_file = part_dir + "/" + std::to_string(processors) + "/" +
                             name + ".bsp";

            auto tree_partitioning =
                tomo::load_partitioning(tree_file, global_volume, log2(p));

            auto main_d = tree_partitioning->splits().root->value.d;

            // which dimension, the first split decides.....
            auto block_partitioning = bulk::block_partitioning<3_D, 1>(
                tomo::math::vec_to_array<3_D, int>(global_volume.voxels()), {p},
                {main_d});

            count_messages(world, block_partitioning,
                           (tomo::geometry::trajectory<3_D, T>&)global_geometry,
                           table, name, "trivial");
            count_messages(world, *tree_partitioning,
                           (tomo::geometry::trajectory<3_D, T>&)global_geometry,
                           table, name, "bisected");
        }

        if (world.rank() == 0) {
            table.print();
            std::ofstream of(outfile, std::ios::out);
            table.print(of);
        }
    });
}

void usage(std::string program_name) {
    std::cout << "Usage: " << program_name << " --geom GEOMS --part PART_DIR "
                                              "--out TABLE_FILE [-k SIZE] [-i "
                                              "ITERS]\n";
}

int main(int argc, char* argv[]) {
    auto opts = options{argc, argv};

    if (!opts.required_arguments({"--geom", "--part", "--out"})) {
        usage(argv[0]);
        return -1;
    }

    auto table = tomo::util::report("Messages", "geometry");
    table.add_column("trivial");
    table.add_column("bisected");

    run(opts.args("--geom"), opts.arg("--part"), opts.arg("--out"), table);

    return 0;
}
