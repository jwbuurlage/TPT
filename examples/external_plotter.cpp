#include <cmath>
#include <iostream>
#include <memory>

#include <glm/gtx/string_cast.hpp>

#include "tomo.hpp"
#include "util/plotter.hpp"

using T = double;
constexpr tomo::dimension D = 3_D;

int main(int argc, char* argv[]) {
    auto opt = tomo::util::args(argc, argv);

    auto v = tomo::volume<D, T>(opt.k);
    auto recon = std::make_shared<tomo::util::dummy_reconstructor<T>>(v);

    // request a plotter scene
    tomo::util::ext_plotter<D, T> plotter("tcp://localhost:5555",
                                          "Sequential serve-and-plot");

	plotter.subscribe("tcp://localhost:5556");
	plotter.serve();

    plotter.set_reconstructor(recon);
    recon->reconstruct(opt.k);

    return 0;
}
