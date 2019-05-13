#include "tpt/tpt.hpp"
using namespace tpt;

int main() {
    using T = float;
    constexpr dimension D = 2_D;

    int size = 16;
    auto v = volume<D, T>(size);
    auto g = geometry::parallel<D, T>(v, size / 2);
    auto f = modified_shepp_logan_phantom<T>(v);
    auto k = dim::joseph<D, T>(v);
    auto p = forward_projection<D, T>(f, g, k);

    int nzs = 0;
    for (auto [idx, line] : g) {
        (void)idx;
        for (auto elem : k(line)) {
            (void)elem;
            nzs++;
        }
    }

    std::cout << "%%MatrixMarket matrix coordinate real general\n";
    std::cout << g.lines() << " " << v.cells() << " " << nzs << "\n";
    for (auto [idx, line] : g) {
        (void)idx;
        for (auto elem : k(line)) {
            std::cout << idx + 1 << " " << elem.index + 1 << " 1.0\n";
        }
    }

    //    auto x = reconstruction::sirt(v, g, k, p);
    //    ascii_plot(x);
}
