#include "tomo.hpp"
#include "util/image_processing.hpp"
#include "util/plotter.hpp"
#include "util/rescale.hpp"
using tomo::img::operator-;
using tomo::img::operator+;
using tomo::img::operator/;
using tomo::img::operator*;

#include <chrono>
using namespace std::chrono_literals;
#include <limits>
#include <thread>

#include "fmt/format.h"

#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

using T = float;
constexpr tomo::dimension D = 3_D;

void show(std::string tag, const tomo::image<D - 1, T>& img,
          tomo::util::ext_plotter<D - 1, T>& plotter) {
    std::cout << tag << ":\n";
    std::this_thread::sleep_for(1s);
    plotter.plot(img);
}

int main(int argc, char* argv[]) {
    auto opt = tomo::util::args(argc, argv);

    try {
        tomo::util::ext_plotter<D, T> plotter("tcp://localhost:5555", "Geometry spec");
        tomo::util::ext_plotter<D - 1, T> proj_plotter("tcp://localhost:5555", "Sinogram");

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
        auto scaled_problem =
            tomo::rescale<D, T>(problem, {64, 64}, {64, 64, 64}, 64);
        auto& projs = scaled_problem.projection_stack;
        show("proj", projs.get_projection(0), proj_plotter);

        auto proj = tomo::dim::joseph<D, T>(scaled_problem.object_volume);

        auto root_dir = fs::path(opt.data).remove_filename().string();
        auto& g = scaled_problem.acquisition_geometry;

        // preprocessing
        auto ref = tomo::tiff_to_image<2_D, T>(root_dir + "0um-100s.tif");
        ref = tomo::img::clamp(ref, (T)10, std::numeric_limits<T>::max());

        auto ref1 = tomo::tiff_to_image<2_D, T>(root_dir + "0um-100s-2.tif");
        ref1 = tomo::img::clamp(ref, (T)10, std::numeric_limits<T>::max());

        auto ref2 = tomo::tiff_to_image<2_D, T>(root_dir + "0um-100s-3.tif");
        ref2 = tomo::img::clamp(ref2, (T)10, std::numeric_limits<T>::max());

        auto flat_field = ref + ref1 + ref2;
        flat_field = tomo::downscale_(flat_field, g->detector_shape());
        show("flatfield", flat_field, proj_plotter);

        for (int i = 0; i < g->projection_count(); ++i) {
            auto flat_cor = (T)300 * (projs.get_projection(i) / flat_field);
            flat_cor = -tomo::img::log(flat_cor);
            projs.set_projection(i, flat_cor);
        }

        fmt::print("STARTING SIRT\n");
        auto z = tomo::reconstruction::sirt(
            scaled_problem.object_volume, *scaled_problem.acquisition_geometry,
            proj, scaled_problem.projection_stack, opt.beta, opt.iterations,
            {[&](tomo::image<D, T>& image) { plotter.plot(image); }});
        fmt::print("SIRT\n");
        tomo::ascii_plot(z);
        proj_plotter.plot(scaled_problem.projection_stack.get_projection(0));
        plotter.plot(z);
    } catch (const std::exception& e) {
        std::cout << "Reading configuration failed: " << e.what() << "\n";
    }

    return 0;
}
