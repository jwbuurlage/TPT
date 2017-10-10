#include "tomos/algorithms/column_action.hpp"
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
    auto count = opts.arg_as_or<uint64_t>("-c", size / 4);

    auto volume = tomo::volume<D, T>(size);
    auto phantom = tomo::modified_shepp_logan_phantom<T>(volume);
    tomo::ascii_plot(phantom);

    auto geometry = tomo::geometry::fan_beam<T>(
        volume, size, math::vec<D, T>(2.5), math::vec<D, int>(size), 4.0, 4.0);
    auto kernel = tomo::dim::joseph<D, T>(volume);
    auto b = tomo::forward_projection(phantom, geometry, kernel);

    if (opts.passed("--cyclic")) {
        // sequential cyclic
        auto x = tomo::reconstruction::column_action_cyclic(
            volume, geometry, kernel, b, beta, sweeps);
        tomo::ascii_plot(x);
    }

    if (opts.passed("--block")) {
        // std::vector<uint64_t> indices(volume.cells());
        // auto block_size = volume.cells() / count;
        // std::iota(indices.begin(), indices.end(), 0);

        // std::random_device rd;
        // std::mt19937 g(rd());

        // for (auto i = 0u; i < block_size; ++i) {
        //     std::shuffle(indices.begin() + i * count,
        //                  indices.begin() + (i + 1) * count, g);
        // }

        // auto block = [&](uint64_t index) -> std::vector<uint64_t> {
        //     auto result = std::vector<uint64_t>(block_size);
        //     for (auto i = 0u; i < block_size; ++i) {
        //         result[i] = indices[index + i * count];
        //     }
        //     return result;
        // };

        auto block_size = volume.cells() / count;

        auto block = [&](uint64_t index) -> std::vector<uint64_t> {
            auto result = std::vector<uint64_t>(block_size);
            for (auto i = 0u; i < block_size; ++i) {
                result[i] = i + index * block_size;
            }
            return result;
        };

        // block
        auto y = tomo::reconstruction::column_action_block(
            volume, geometry, kernel, b, block, count, beta, sweeps);
        tomo::ascii_plot(y);
    }
}
