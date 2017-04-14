#include "tomos/tomos.hpp"
#include "tomos/util/plotter.hpp"
#include "tomos/util/rescale.hpp"

#include "fmt/format.h"

using T = float;
constexpr tomo::dimension D = 3_D;

int main(int argc, char* argv[]) {
    auto opt = tomo::util::args(argc, argv);

    try {
        tomo::util::ext_plotter<D - 1, T> proj_plotter("tcp://localhost:5555", "Sinogram");
        tomo::util::ext_plotter<D, T> plotter ("tcp://localhost:5555", "Geometry spec");

        auto problem = tomo::read_configuration<D, T>(opt.data);

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

        // std::cout << "done loading, now should downscale\n";
        proj_plotter.plot(problem.projection_stack.get_projection(0));
        auto scaled_problem =
            tomo::rescale<D, T>(problem, {64, 64}, {64, 64, 64}, 64);

        auto proj = tomo::dim::joseph<D, T>(scaled_problem.object_volume);

        fmt::print("STARTING SIRT\n");
        auto z = tomo::reconstruction::sirt(
            scaled_problem.object_volume, *scaled_problem.acquisition_geometry,
            proj, scaled_problem.projection_stack, opt.beta, opt.iterations, {[&](tomo::image<D, T>& image) { plotter.plot(image); }});
        fmt::print("SIRT\n");
        tomo::ascii_plot(z);
        proj_plotter.plot(scaled_problem.projection_stack.get_projection(0));
        plotter.plot(z);
    } catch (const std::exception& e) {
        std::cout << "Reading configuration failed: " << e.what() << "\n";
    }

    return 0;
}
