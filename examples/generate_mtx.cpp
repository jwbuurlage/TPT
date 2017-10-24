#include "tomos/tomos.hpp"
using namespace tomo;

int main() {
    using T = float;
    constexpr dimension D = 3_D;

    int size = 5;
    auto v = tomo::volume<D, T>(size);
    auto g = tomo::geometry::parallel<D, T>(v, 4);
    auto f = tomo::modified_shepp_logan_phantom<T>(v);
    auto k = tomo::dim::joseph<D, T>(v);
    auto p = tomo::forward_projection<D, T>(f, g, k);

    int nzs = 0;
    for (auto [idx, line] : g) {
        (void)idx;
        for (auto elem : k(line)) {
            (void)elem;
            nzs++;
        }
    }

    std::cout << "%%MatrixMarket matrix coordinate integer general\n";
    std::cout << g.lines() << " " << v.cells() << " " << nzs << "\n";
    for (auto [idx, line] : g) {
        (void)idx;
        for (auto elem : k(line)) {
            std::cout << idx << " " << elem.index << " 0\n";
        }
    }

    //    auto x = tomo::reconstruction::sirt(v, g, k, p);
    //    tomo::ascii_plot(x);
}
