#include "tomos/tomos.hpp"
#include "tomos/util/column_iterator.hpp"
#include "tomos/util/simple_args.hpp"
using namespace tomo;

#include <random>

using namespace std::string_literals;

int main(int argc, char* argv[]) {
    using T = float;

    auto opts = tomo::options{argc, argv};

    constexpr dimension D = 2_D;

    auto size = opts.arg_as_or<int>("-s", 32);
    auto sweeps = opts.arg_as_or<int>("-i", 10);
    auto beta = opts.arg_as_or<T>("--beta", (T)1);

    auto volume = tomo::volume<D, T>(size);
    auto phantom = tomo::modified_shepp_logan_phantom<T>(volume);
    tomo::ascii_plot(phantom);

    auto geometry = tomo::geometry::fan_beam<T>(
        volume, size, math::vec<D, T>(2.5), math::vec<D, int>(size), 5.0, 5.0);
    auto kernel = tomo::dim::joseph<D, T>(volume);
    auto b = tomo::forward_projection(phantom, geometry, kernel);

    //    auto p = tomo::forward_projection<D, T>(phantom, geometry, kernel);
    //    auto x = tomo::reconstruction::sirt(volume, geometry, kernel, p);
    //
    //    tomo::ascii_plot(x);
    //    tomo::write_png(x, "fan_beam");

    auto column = tomo::column<D, T>(geometry, kernel);

    auto cs = tomo::image<D, T>(volume);
    for (auto j = 0u; j < volume.cells(); ++j) {
        auto voxel_index = volume.unroll(j);
        for (auto[line_idx, value] : column(voxel_index)) {
            (void)line_idx;
            cs[j] += value * value;
        }
    }

    auto x = tomo::image<D, T>(volume);
    auto r = b;

    auto js = std::vector<uint64_t>(volume.cells());
    std::iota(js.begin(), js.end(), 0);

//    std::random_device rd;
//    std::mt19937 g(rd());
//    std::shuffle(js.begin(), js.end(), g);

    // sequential cyclic
    for (auto k = 0; k < sweeps; ++k) {
        for (auto q = 0u; q < volume.cells(); ++q) {
            auto j = js[q];
            auto voxel_index = volume.unroll(j);
            if (cs[j] < math::epsilon<T>) {
                continue;
            }
            auto delta = (T)0;
            for (auto[line_idx, value] : column(voxel_index)) {
                delta += value * r[line_idx];
            }
            delta /= cs[j];
            delta *= beta;
            for (auto[line_idx, value] : column) {
                r[line_idx] -= value * delta;
            }
            x[j] += delta;
        }
        tomo::write_png(x, "fan_beam_"s + std::to_string(k));
        std::cout << "residual norm (" << k << "): " << math::norm(r) << "\n";
    }

    // print result
    tomo::ascii_plot(x);
}
