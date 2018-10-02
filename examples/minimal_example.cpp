#include "tomos/tomos.hpp"

using namespace tomo;

int main() {
    using T = float;
    constexpr dimension D = 2_D;

    int size = 128;
    auto v = tomo::volume<D, T>(size);
    auto g = tomo::geometry::fan_beam<T>(v, size, math::vec<D, T>(2.5),
                                         math::vec<D, int>(size), 4.0, 4.0);
    auto k = tomo::dim::joseph<D, T>(v);
    auto f = tomo::modified_shepp_logan_phantom<T>(v);
    auto p = tomo::forward_projection<D, T>(f, g, k);

    auto x = tomo::reconstruction::sirt(v, g, k, p);
    tomo::ascii_plot(x);
}
