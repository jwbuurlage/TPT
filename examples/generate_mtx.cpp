#include <CLI/CLI.hpp>

#include "tpt/tpt.hpp"
using namespace tpt;

#include <filesystem>
#include <string>
namespace fst = std::filesystem;
using namespace std::string_literals;

int main(int argc, char** argv) {
    CLI::App app{"Generate a .mtx from a configuration"};

    using T = float;
    constexpr dimension D = 3_D;

    std::vector<std::string> geometry_files;

    app.add_option("--geometry", geometry_files,
                   "the input geometries (.toml files)")
        ->required();

    CLI11_PARSE(app, argc, argv);

    for (auto geometry_file : geometry_files) {
        auto name = fst::path(geometry_file).stem().string();
        auto problem = tpt::read_configuration<3_D, T>(geometry_file);

        auto v = problem.object_volume;
        auto& g = *problem.acquisition_geometry;
        auto kernel = dim::joseph<D, T>(v);

        {
            auto fout = std::ofstream(name + ".mtx");

            int nzs = 0;
            for (auto [idx, line] : g) {
                (void)idx;
                for (auto elem : kernel(line)) {
                    (void)elem;
                    nzs++;
                }
            }

            fout << "%%MatrixMarket matrix coordinate real general\n";
            fout << g.lines() << " " << v.cells() << " " << nzs << "\n";
            for (auto [idx, line] : g) {
                (void)idx;
                for (auto elem : kernel(line)) {
                    fout << idx + 1 << " " << elem.index + 1 << " "
                         << elem.value << "\n";
                }
            }
        }
    }
}
