#include "tomo.hpp"

using T = float;
constexpr tomo::dimension D = 3_D;

int main(int argc, char* argv[]) {
    auto opt = tomo::util::args(argc, argv);

    try {
        auto geom_and_volume =
            tomo::read_configuration<D, T>("data/geometries/parallel.toml");

        auto vol = geom_and_volume.second;
        auto& geom = *geom_and_volume.first;

        auto phantom = tomo::modified_shepp_logan_phantom<T>(vol);
        auto proj = tomo::dim::joseph<D, T>(vol);
        auto sino = tomo::forward_projection<D, T>(phantom, geom, proj);

        tomo::ascii_plot(phantom, 3, 1);
    } catch (const std::exception& e) {
        std::cout << "Reading configuration failed: " << e.what() << "\n";
    }

    return 0;
}
