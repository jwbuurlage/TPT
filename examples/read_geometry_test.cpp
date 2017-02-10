#include "tomo.hpp"
#include "util/plotter.hpp"

#include "fmt/format.h"

using T = float;
constexpr tomo::dimension D = 3_D;

int main(int argc, char* argv[]) {
    auto opt = tomo::util::args(argc, argv);

    try {
        auto plotter =
            tomo::ext_plotter<D, T>("tcp://localhost:5555", "Geometry spec");

        auto problem = tomo::read_configuration<D, T>(opt.data);

        auto proj = tomo::dim::joseph<D, T>(problem.object_volume);


        /*
         * We want to solve a reconstruction problem, the following things are
         * known from the file:
         * - geometry
         * - volume
         * - projections
         *
         * the other choices to make are:
         * - algorithm
         * - interpolation method
         * - result plotting/storing
         * - (for distributed) the partitioning to use
         *
         * construct a method that lets you choose these final 'policies'.
         */

        auto z = tomo::reconstruction::sirt(problem.object_volume,
                                            *problem.acquisition_geometry, proj,
                                            problem.projection_stack, opt.beta, opt.iterations);
        fmt::print("SIRT\n");
        tomo::ascii_plot(z);
        plotter.plot(z);
    } catch (const std::exception& e) {
        std::cout << "Reading configuration failed: " << e.what() << "\n";
    }

    return 0;
}
