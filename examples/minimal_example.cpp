#include "tpt/tpt.hpp"

using namespace tpt;

int main() {
    using T = float;
    constexpr dimension D = 2_D;

    int size = 128;
    auto v = volume<D, T>(size);
    auto g = geometry::fan_beam<T>(v, size, math::vec<D, T>(2.5),
                                         math::vec<D, int>(size), 4.0, 4.0);
    auto k = dim::joseph<D, T>(v);
    auto f = modified_shepp_logan_phantom<T>(v);
    auto p = forward_projection<D, T>(f, g, k);

    auto x = reconstruction::sirt(v, g, k, p);
    ascii_plot(x);
}
