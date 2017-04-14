#include <cmath>
#include <iostream>
#include <memory>

#include <glm/gtx/string_cast.hpp>

#include "tomos/tomos.hpp"
#include "tomos/util/plotter.hpp"

using T = float;
constexpr tomo::dimension D = 3_D;

int main(int argc, char* argv[]) {
    auto opt = tomo::util::args(argc, argv);

    int size = opt.k;
    auto v = tomo::volume<D, T>(opt.k);
    auto g = tomo::geometry::cone_beam<T>(v, size, {(T)2.0, (T)2.0},
                                          {size, size}, (T)4.0, (T)4.0);
    auto f = tomo::modified_shepp_logan_phantom<T>(v);
    auto kernel = tomo::dim::joseph<3_D, T>(v);
    auto ps = tomo::forward_projection<3_D, T>(f, g, kernel);

    auto recon = std::make_shared<tomo::util::dummy_reconstructor<T>>(v, kernel,
                                                                      g, f, ps);

    // request a plotter scene
    tomo::util::ext_plotter<D, T> plotter("tcp://localhost:5555",
                                          "Sequential serve-and-plot");

    plotter.send_projection_data(g, ps, v);

    plotter.subscribe("tcp://localhost:5556");
    plotter.serve();

    plotter.set_reconstructor(recon);
    recon->reconstruct();

    return 0;
}
