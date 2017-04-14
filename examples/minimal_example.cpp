#include "tomos/tomos.hpp"

int main() {
    using T = float;
    constexpr tomo::dimension D = 2_D;

    int size = 128;
    auto v = tomo::volume<D, T>(size);
    auto g = tomo::geometry::parallel<D, T>(v, size, size);
    auto f = tomo::modified_shepp_logan_phantom<T>(v);
    auto k = tomo::dim::joseph<D, T>(v);
    auto p = tomo::forward_projection<D, T>(f, g, k);

    auto x = tomo::reconstruction::sirt(v, g, k, p);
    tomo::ascii_plot(x);
}
