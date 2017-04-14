#include "tomos/tomos.hpp"
#include "tomos/util/simple_args.hpp"
#include "tomos/util/trees.hpp"

#include <glm/gtx/string_cast.hpp>

using T = float;
constexpr tomo::dimension D = 3_D;

void run(std::string meta_file, std::string partition_file, int size,
         tomo::util::report& table) {
    auto problem = tomo::read_configuration<3_D, T>(meta_file);

    // rescale
    problem.acquisition_geometry->set_detector_shape({size, size});
    problem.acquisition_geometry->set_projections(size);

    T scale = (T)size / problem.object_volume.voxels()[0];
    auto new_voxels = tomo::math::vec3<int>(
        scale * tomo::math::vec3<T>(problem.object_volume.voxels()));
    problem.object_volume.set_voxels(new_voxels);

    std::cout << "new_voxels = " << glm::to_string(new_voxels) << "\n";

    // read partitioning
    auto reloaded = tomo::deserialize_tree<T>(partition_file);
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

    auto name = fs::path(meta_file).stem();
    table.add_row(name);
    table.add_result(name, "trivial", overlap_trivial);
    table.add_result(name, "binary", overlap_bisected);
    table.add_result(name, "improvement", fmt::format("{:.1f}%", 100 * imp));
}

int main(int argc, char* argv[]) {
    auto opts = tomo::options{argc, argv};

    int leading_size = 256;
    if (opts.passed("-s")) {
        leading_size = opts.arg_as<int>("-s");
    }

    if (!opts.required_arguments({"--in", "--parts"})) {
        std::cout << "Input and/or partitioning files not given.\n";
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
    table.add_column("improvement");

    for (auto in_and_part : ins_and_parts) {
        run(in_and_part.first, in_and_part.second, leading_size, table);
    }

    table.print();

    return 0;
}
