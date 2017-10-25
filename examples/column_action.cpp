#include "tomos/algorithms/column_action.hpp"
#include "tomos/tomos.hpp"
#include "tomos/util/column_iterator.hpp"
#include "tomos/util/matrix_sums.hpp"
#include "tomos/util/simple_args.hpp"

using namespace tomo;
using namespace tomo::img;

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
    auto output_base = opts.arg_or("--output", "fan_beam_");

    auto volume = tomo::volume<D, T>(size);
    auto phantom = tomo::modified_shepp_logan_phantom<T>(volume);
    tomo::ascii_plot(phantom);

    auto geometry = tomo::geometry::fan_beam<T>(
        volume, size, math::vec<D, T>(2.5),
        math::vec<D, int>((int)(1.5 * size)), 4.0, 4.0);
    auto kernel = tomo::dim::joseph<D, T>(volume);
    auto b = tomo::forward_projection(phantom, geometry, kernel);

    // sequential cyclic
    auto rs = tomo::row_sums(geometry, kernel);
    auto cs = tomo::column_sums(geometry, kernel);

    index_space* idxs = nullptr;
    if (opts.passed("--idxs-reverse")) {
        idxs = new reverse_index_space{(int)volume.cells()};
    } else if (opts.passed("--idxs-back-forth")) {
        idxs = new back_forth_index_space{(int)volume.cells()};
    } else if (opts.passed("--idxs-random")) {
        idxs = new random_index_space{(int)volume.cells()};
    } else if (opts.passed("--idxs-hilbert")) {
        idxs = new hilbert_index_space(size);

        std::vector<int> results;
        for (auto i = 0u; i < volume.cells(); ++i) {
            results.push_back((*idxs)(0, i));
        }
        std::sort(results.begin(), results.end());
        std::vector<int> results2(volume.cells());
        std::iota(results2.begin(), results2.end(), 0);
        assert(results == results2);
    } else {
        idxs = new index_space();
    }

    auto t = tomo::math::dot(rs, b) / tomo::math::dot(rs, rs);
    auto x0 =
        tomo::image<D, T>(volume, opts.passed("--init-constant") ? t : (T)0);
    if (opts.passed("--init-sirt")) {
        for (auto& r : rs) {
            r = (math::abs(r) > math::epsilon<T>) ? ((T)1.0 / r) : (T)0.0;
        }
        for (auto& bc : cs) {
            bc = (math::abs(bc) > math::epsilon<T>) ? ((T)1.0 / bc) : (T)0.0;
        }
        tomo::reconstruction::iterate::sirt(x0, geometry, kernel, b, rs, cs);
    }

    if (opts.passed("--cyclic")) {
        auto err = math::norm(phantom);
        auto rfout = std::ofstream(output_base + "residuals.md");
        auto efout = std::ofstream(output_base + "rel_errors.md");
        auto x = tomo::reconstruction::column_action_cyclic(
            volume, geometry, kernel, b, beta, sweeps, {x0}, {idxs},
            {[&](const image<D, T>& xk, int k, const projections<D, T>& rk) {
                tomo::write_png(xk, output_base + std::to_string(k));
                rfout << k << " " << math::norm(rk) << "\n";
                auto derr = math::norm(phantom - xk);
                efout << k << " " << derr / err << "\n";
            }});
        tomo::ascii_plot(x);
    }

    if (opts.passed("--block")) {
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

    if (opts.passed("--generate-ata")) {
        int n = (int)volume.cells();
        auto ata = tomo::image<2_D, T>({n});

        // do something
        for (auto[idx, line] : geometry) {
            (void)idx;
            for (auto i : kernel(line)) {
                for (auto j : kernel) {
                    ata({(*idxs)(0, i.index), (*idxs)(0, j.index)}) +=
                        i.value * j.value;
                }
            }
        }

        tomo::ascii_plot(ata);
        tomo::write_png(ata, output_base + "ata");
    }

    delete idxs;
}
