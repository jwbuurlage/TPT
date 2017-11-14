/*
 * FIXME:
 * - Can we group sends together? this should greatly reduce data size
 *   with a constant factor.
 */

#include <cassert>
#include <sstream>

#include <bulk/backends/mpi/mpi.hpp>
#include <bulk/bulk.hpp>

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
        my_contributions(t).send(std::move(a_targets[t]),
                                 std::move(a_locals[t]),
                                 std::move(a_remotes[t]));
        my_responsibilities(t).send(std::move(b_targets[t]),
                                    std::move(b_locals[t]),
                                    std::move(b_remotes[t]));
    }
    world.sync();

    // 5. store contributors and owners in a vector
    auto contributions = std::vector<std::vector<pixel_message>>(p);
    auto responsibilities = std::vector<std::vector<pixel_message>>(p);

    for (auto& [owners, locals, remotes] : my_contributions) {
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

void count_messages(bulk::world& world,
          bulk::rectangular_partitioning<3, 1>& partitioning,
          tomo::volume<3_D, T> global_volume,
          geometry::trajectory<3_D, T>& global_geometry,
          tomo::util::report& table, std::string name, std::string column,
          std::string image_dir, int iters) {
(void)global_volume;
(void)image_dir;
(void)iters;
auto p = world.active_processors();

    if (world.rank() == 0) {
        world.log("Running %s (%s)", name.c_str(), column.c_str());
    }

    auto vs = calculate_local_volume(partitioning, world.rank());
    auto gs = distributed::restricted_geometry<T>(global_geometry, vs);

    auto result = compute_contributions(world, gs);
    auto& go_forth = result.first;
    auto& and_back = result.second;

    auto x = bulk::var<int>(world, 0);
for (int t = 0; t < p; ++t) {
	if (!go_forth[t].empty()) {
	x += 1;
	}
	if (!and_back[t].empty()) {
	x += 1;
        }
    }
    auto messages = bulk::foldl(x, [](int& lhs, int rhs) { lhs += rhs; });

    if (world.rank() == 0) {
        world.log("messages: %i", messages);
    }

    table.add_result(name, column, messages);
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
                count_messages(world, block_partitioning, global_volume,
                     (tomo::geometry::trajectory<3_D, T>&)global_geometry,
                     table, name, "trivial", image_dir, iters);
            }
            if (bisected) {
                count_messages(world, *tree_partitioning, global_volume,
                     (tomo::geometry::trajectory<3_D, T>&)global_geometry,
                     table, name, "bisected", image_dir, iters);
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
                 "ITERS] [--trivial] [--bisected]\n";
}

int main(int argc, char* argv[]) {
    auto opts = options{argc, argv};

    if (!opts.required_arguments({"--geom", "--part", "--out"})) {
        usage(argv[0]);
        return -1;
    }

    auto table = tomo::util::report("Runtimes", "geometry");
    table.add_column("trivial");
    table.add_column("trivial (com)");
    table.add_column("bisected");
    table.add_column("bisected (com)");

    run(opts.args("--geom"), opts.arg("--part"), opts.arg_as_or<int>("-k", -1),
        opts.arg_as_or<int>("-i", 1), opts.arg("--out"), "fake",
        table, opts.passed("--trivial"), opts.passed("--bisected"));

    return 0;
}
