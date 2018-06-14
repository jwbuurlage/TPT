#include "tomos/tomos.hpp"
#include "tomos/util/simple_args.hpp"
#include "tomos/util/trees.hpp"

#include <glm/gtx/string_cast.hpp>

using T = float;
constexpr tomo::dimension D = 3_D;

void run(std::string meta_file, std::string partition_file, int size, int p,
         tomo::util::report& table) {
    auto q = log2(p);

    auto problem = tomo::read_configuration<3_D, T>(meta_file, size);

    T scale = (T)size / problem.object_volume.voxels()[0];
    auto new_voxels = tomo::math::vec3<int>(
        scale * tomo::math::vec3<T>(problem.object_volume.voxels()));
    problem.object_volume.set_voxels(new_voxels);

    std::cout << "new_voxels = " << glm::to_string(new_voxels) << "\n";

    // read partitioning
    auto reloaded = tomo::deserialize_tree<T>(partition_file, q);
    auto tree = tomo::from_neutral_tree(*reloaded, problem.object_volume);

    auto processors = tree.size();

    auto part_bisected = bulk::tree_partitioning<D>(
        tomo::math::vec_to_array<D, int>(problem.object_volume.voxels()),
        processors, std::move(tree));

    // do trivial partitioning
    auto part_trivial = tomo::distributed::partition_trivial(
        *problem.acquisition_geometry, problem.object_volume, processors);

    // compare comvolume
    auto overlap_trivial = tomo::distributed::communication_volume<D, T>(
        *problem.acquisition_geometry, problem.object_volume, part_trivial);
    auto overlap_bisected = tomo::distributed::communication_volume<D, T>(
        *problem.acquisition_geometry, problem.object_volume, part_bisected);
    T imp = (T)0.0;
    if (overlap_trivial != 0)
        imp = (overlap_trivial - overlap_bisected) / (T)overlap_trivial;


    auto imba_trivial = tomo::distributed::load_imbalance<D, T>(
        problem.object_volume, part_trivial, *problem.acquisition_geometry);
    auto imba_bisected = tomo::distributed::load_imbalance<D, T>(
        problem.object_volume, part_bisected, *problem.acquisition_geometry);

    auto name = fs::path(meta_file).stem();
    table.add_row(name);
    table.add_result(name, "trivial", overlap_trivial);
    table.add_result(name, "binary", overlap_bisected);
    table.add_result(name, "trivial_e", imba_trivial);
    table.add_result(name, "binary_e", imba_bisected);
    table.add_result(name, "improvement", fmt::format("{:.1f}%", 100 * imp));
}

void usage(std::string program_name) {
    std::cout << "USAGE: " << program_name
              << " --in GEOMS --parts PART_DIR [-s SIZE]\n";
}

int main(int argc, char* argv[]) {
    auto opts = tomo::options{argc, argv};

    int p = 4;
    int leading_size = 256;
    if (opts.passed("-s")) {
        leading_size = opts.arg_as<int>("-s");
    }

    if (opts.passed("-p")) {
        p = opts.arg_as<int>("-p");
    }


    if (!opts.required_arguments({"--in", "--parts"})) {
        usage(argv[0]);
        return -1;
    }

    std::vector<std::pair<std::string, std::string>> ins_and_parts;
    if (opts.required_arguments({"--in", "--parts"})) {
        auto ins = opts.args("--in");
        auto part_dir = opts.arg("--parts");
        auto parts = ins;
        std::transform(parts.begin(), parts.end(), parts.begin(),
                       [&](std::string in) {
                           auto stem = fs::path(in).stem();
                           return part_dir + stem.string() + ".bsp";
                       });

        ins_and_parts = tomo::zip(ins, parts);
    }

    auto table = tomo::util::report(
        "Detector overlaps for geometric partitioning", "geometry");
    table.add_column("trivial");
    table.add_column("binary");
    table.add_column("trivial_e");
    table.add_column("binary_e");
    table.add_column("improvement");

    for (auto in_and_part : ins_and_parts) {
        run(in_and_part.first, in_and_part.second, leading_size, p, table);
    }

    table.print();

    return 0;
}
