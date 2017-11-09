#include "tomos/tomos.hpp"
#include "tomos/util/image_processing.hpp"
#include "tomos/util/matrix_sums.hpp"
#include "tomos/util/simple_args.hpp"
#include "tomos/util/svd.hpp"

#include <string>

using namespace std::string_literals;
using namespace tomo;
// enable image operations
using namespace tomo::img;

int main(int argc, char* argv[]) {
    using T = double;
    constexpr dimension D = 2_D;
    using DIM = tomo::dim::joseph<D, T>;

    // ------------------------------------------------------------------------
    // SETUP THE PROBLEM
    // CLI options
    auto opts = tomo::options{argc, argv};
    auto size = opts.arg_as_or<int>("-s", 32);
    auto iter = opts.arg_as_or<int>("-i", 10);
    auto output_base = opts.arg_or("--output", "roi_");

    auto v = tomo::volume<D, T>(size);
    //    auto g = tomo::geometry::fan_beam<T>(
    //        v, size, math::vec<D, T>(2.5), math::vec<D, int>((int)((T)1.5 *
    //        size)),
    //        4.0, 4.0);
    auto g = tomo::geometry::parallel<D, T>(v, size);

    auto k = DIM(v);
    // next, we try to do a reconstruction in this basis
    auto f = tomo::modified_shepp_logan_phantom<T>(v);
    auto b = tomo::forward_projection<D, T>(f, g, k);

    // LOCAL VOLUME
    // 1. MAKE LOCAL VOLUME
    auto lv = tomo::volume<D, T>(math::vec<D, int>(size / 2, size / 4),
                                 math::vec<D, T>((T)0.25, (T)0.1),
                                 math::vec<D, T>((T)0.5, (T)0.25));
    (void)lv;

    // 2. MASK THE RECONSTRUCTION ON THE LOCAL VOLUME
    auto mask = [&](tomo::image<D, T> img) -> tomo::image<D, T> {
        auto imgsize = img.get_volume().voxels()[0];
        for (int i = lv.origin()[0] * imgsize;
             i < (lv.origin()[0] + lv.physical_lengths()[0]) * imgsize; ++i) {
            for (int j = lv.origin()[1] * imgsize;
                 j < (lv.origin()[1] + lv.physical_lengths()[1]) * imgsize;
                 ++j) {
                img({i, j}) = (T)0;
            }
        }
        return img;
    };

    auto circle = [&](tomo::image<D, T> img) -> tomo::image<D, T> {
        auto center = math::vec<D, T>(lv.voxels() / 2);
        for (int i = 0; i < lv.voxels()[0]; ++i) {
            for (int j = 0; j < lv.voxels()[1]; ++j) {
                auto pos = math::vec<D, T>(i, j);
                if (math::norm<D, T>(pos - center) >
                    math::max(lv.voxels()[0] / 2, lv.voxels()[1] / 2)) {
                    img({i, j}) = (T)0;
                }
            }
        }
        return img;
    };

    auto lk = DIM(lv);

    // ------------------------------------------------------------------------
    // BASELINE: IGNORING OUTSIDE INFLUENCE

    if (opts.passed("--baseline")) {
        auto lx = tomo::reconstruction::sirt(lv, g, lk, b, 1.0, iter, {}, true,
                                             (T)-1, (T)1);
        tomo::ascii_plot(circle(lx));
        tomo::write_png(circle(lx), output_base + "baseline"s);
    }

    // -----------------------------------------------------------------------
    // CONSTANT APPROXIMATION TO OUTSIDE

    if (opts.passed("--constant")) {
        auto rs = tomo::row_sums(g, k);
        auto t = tomo::math::dot(rs, b) / tomo::math::dot(rs, rs);
        auto xout = mask(tomo::image<D, T>(v, t));
        auto lb = b - tomo::forward_projection<D, T>(xout, g, k);

        auto lx = tomo::reconstruction::sirt(lv, g, lk, lb, 1.0, iter, {}, true,
                                             (T)-1, (T)1);
        tomo::ascii_plot(circle(lx));
        tomo::write_png(circle(lx), output_base + "constant"s);
    }

    // ------------------------------------------------------------------------
    // SVD APPROACH
    if (opts.passed("--svd")) {
        auto rank = opts.arg_as_or<int>("-r", 10);

        auto[U, S, V] = tomo::util::truncated_svd(g, k, v, rank);

        for (auto& s : S) {
            s = (T)1 / s;
        }

        tomo::image<D, T> result(v);
        for (int i = 0; i < rank; ++i) {
            result = result + ((S[i] * math::dot(b, U[i])) * V[i]);
        }

        if (opts.passed("--save-svd")) {
            for (int i = 0; i < rank; ++i) {
                tomo::write_png(V[i],
                                output_base + "svd_"s + std::to_string(i));
            }
            tomo::ascii_plot(result);
            tomo::write_png(result, output_base + "svd_recon"s);
        }

        auto err = math::norm(f);
        auto derr = math::norm(f - result);
        std::cout << "Rel Error: " << derr / err << "\n";
        auto residual =
            math::norm(b - tomo::forward_projection<D, T>(result, g, k));
        std::cout << "Rel Residual: " << residual / math::norm(b) << "\n";

        auto lb = b - tomo::forward_projection<D, T>(mask(result), g, k);
        auto lx = tomo::reconstruction::sirt(lv, g, lk, lb, 1.0, iter, {}, true,
                                             (T)-1, (T)1);
        tomo::ascii_plot(circle(lx));
        tomo::write_png(circle(lx), output_base + "with_svd"s);
    }
    // -------------------------------------------------------------------------
    // MULTIGRID APPROACH

    if (opts.passed("--multi-grid")) {
        auto coarse_size = opts.arg_as_or<int>("-cs", size / 8);
        auto cv = tomo::volume<D, T>(coarse_size);
        auto ck = DIM(cv);
        auto cx = tomo::reconstruction::sirt(cv, g, ck, b, 1.0, iter);
        tomo::write_png(cx, output_base + "coarse"s);

        auto lb = b - tomo::forward_projection<D, T>(mask(cx), g, ck);
        auto lx = tomo::reconstruction::sirt(lv, g, lk, lb, 1.0, iter, {}, true,
                                             (T)-1, (T)1);
        tomo::ascii_plot(circle(lx));
        tomo::write_png(circle(lx), output_base + "with_multigrid"s);
    }

    // -------------------------------------------------------------------------
    // FBP APPROCH

    if (opts.passed("--fbp")) {
        // TODO
    }
}
