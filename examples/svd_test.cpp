#include "tomos/tomos.hpp"
#include "tomos/util/simple_args.hpp"
#include "tomos/util/svd.hpp"
#include "tomos/util/image_processing.hpp"

#include <string>

using namespace std::string_literals;
using namespace tomo;
// enable image operations
using namespace tomo::img;

int main(int argc, char* argv[]) {
    using T = float;
    constexpr dimension D = 2_D;

    // CLI options
    auto opts = tomo::options{argc, argv};
    auto size = opts.arg_as_or<int>("-s", 32);
    auto rank = opts.arg_as_or<int>("-r", 10);

    auto v = tomo::volume<D, T>(size);
    auto g = tomo::geometry::fan_beam<T>(v, size, math::vec<D, T>(2.5),
                                         math::vec<D, int>(size), 4.0, 4.0);
    auto k = tomo::dim::joseph<D, T>(v);

    auto[U, S, V] = tomo::util::truncated_svd(g, k, v, rank);

    (void)V;
    for (int i = 0; i < rank; ++i) {
        std::cout << "sigma_i: " << S[i] << "\n";
        tomo::write_png(V[i], "svd_"s + std::to_string(i));
    }

    for (auto& s : S) {
        s = (T)1 / s;
    }

    // next, we try to do a reconstruction in this basis
    auto f = tomo::modified_shepp_logan_phantom<T>(v);
    auto b = tomo::forward_projection<D, T>(f, g, k);

    tomo::image<D, T> result(v);
    for (int i = 0; i < rank; ++i) {
        result = result + ((S[i] * math::dot(b, U[i])) * V[i]);
    }

    tomo::write_png(result, "svd_recon"s);
    tomo::ascii_plot(result);

    auto err = math::norm(f);
    auto derr = math::norm(f - result);
    std::cout << "Error: " << derr / err << "\n";
}
