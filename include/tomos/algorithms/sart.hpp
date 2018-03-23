#pragma once

#include <iostream>
#include <vector>

#include "../logging.hpp"

namespace tomo {
namespace reconstruction {

/**
 * The Simultaneous Algebraic Reconstruction Technique (SART), is a tomographic
 * reconstruction method based on the classic Kaczmarz method for inverse
 * problems, but performed in blocks.
 *
 * \tparam D the dimension of the problem
 * \tparam T the scalar type in use
 * \tparam Geometry the type of geometry of the problem
 * \tparam Projector the discrete integration method to use
 *
 * \param v the volume of the imaged object
 * \param g the geometry of the problem
 * \param p the measurements (projections)
 * \param beta (optional) a relaxation parameter
 * \param iterations (optional) the number of iterations to perform
 *
 * \returns An image object representing the reconstructed object.
 */
template <dimension D, typename T>
image<D, T> sart(const volume<D, T>& v, const tomo::geometry::base<D, T>& g,
                 tomo::dim::base<D, T>& kernel, const projections<D, T>& p,
                 double beta = 0.5, int iterations = 10,
                 std::function<void(image<D, T>&, int)> callback = {}) {
    image<D, T> f(v);

    // the size of a single block
    int k = math::reduce<D - 1>(g.projection_shape(0));

    // compute $w_i \cdot w_i$
    std::vector<T> w_norms(g.lines());
    for (auto[line_number, line] : g) {
        for (auto elem : kernel(line)) {
            w_norms[line_number] += elem.value * elem.value;
        }
    }

    auto f_next = f;

    for (int iter = 0; iter < iterations; ++iter) {
        int s = k;
        int t = 0;

        for (auto[row, line] : g) {
            if (s == k) {
                // we now update the image
                if (t > 0) {
                    f = f_next;
                }
                ++t;
                s = 0;
            }

            if (w_norms[row] > math::epsilon<T>) {
                T alpha = 0.0;
                for (auto elem : kernel(line)) {
                    alpha += f[elem.index] * elem.value;
                }

                auto factor = beta * ((p[row] - alpha) / w_norms[row]);
                for (auto elem : kernel(line))
                    f_next[elem.index] += factor * elem.value;
            }

            ++s;
        }

        if (callback) {
            callback(f_next, iter);
        }
    }

    return f;
}

// SART with CG as precon
template <dimension D, typename T>
image<D, T> sart_cg(const volume<D, T>& v, const tomo::geometry::base<D, T>& g,
                    tomo::dim::base<D, T>& kernel, const projections<D, T>& b,
                    int iterations = 10, T stepsize = 1.0,
                    std::function<void(image<D, T>&, int)> callback = {}) {
    std::vector<T> w_norms(g.lines());
    for (auto[line_number, line] : g) {
        for (auto elem : kernel(line)) {
            w_norms[line_number] += elem.value * elem.value;
        }
    }

    auto x = image<D, T>(v);
    int block_size = math::reduce<D - 1>(g.projection_shape(0));

    auto xb = image<D, T>(v);
    auto yb = std::vector<T>(block_size);
    auto zb = std::vector<T>(block_size);

    auto fp_block = [&](int block, image<D, T>& z) -> std::vector<T>& {
        std::fill(yb.begin(), yb.end(), 0);

        for (auto it = g.iter_proj(block); it != g.iter_proj(block + 1); ++it) {
            auto[row, line] = *it;
            for (auto elem : kernel(line)) {
                yb[row] += z[elem.index] * elem.value;
            }
        }

        return yb;
    };
    auto bp_block = [&](int block, std::vector<T>& y) -> image<D, T>& {
        std::fill(xb.begin(), xb.end(), 0);

        for (auto it = g.iter_proj(block); it != g.iter_proj(block + 1); ++it) {
            auto[row, line] = *it;
            for (auto elem : kernel(line)) {
                xb[elem.index] += y[row] * elem.value;
            }
        }

        return xb;
    };

    auto cg = [&](int block, std::vector<T>& bb, int initer) {
        auto y = std::vector<T>(block_size);
        auto r = bb;
        auto p = r;

        for (int k = 0; k < initer; ++k) {
            auto gamma = math::dot(r, r);
            auto bpp = bp_block(block, p);
            auto alpha = gamma / math::dot(bpp, bpp);
            for (int i = 0; i < block_size; ++i) {
                y[i] += alpha * p[i];
            }
            for (int i = 0; i < block_size; ++i) {
                auto t = fp_block(block, bpp);
                r[i] -= alpha * t[i];
            }
            auto beta = math::dot(r, r) / gamma;
            for (int i = 0; i < block_size; ++i) {
                p[i] = r[i] + beta * p[i];
            }
        }

        return y;
    };

    for (int k = 0; k < iterations; ++k) {
        for (int block = 0; block < g.projection_count(); ++block) {
            // compute residual for this block
            auto& alphas = fp_block(block, x);
            for (int i = 0; i < block_size; ++i) {
                alphas[i] = b[block * block_size + i] - alphas[i];
            }
            // perform some iterations of CG (i.e. finding how to adjust x to
            // reduce the residual for this block)
            auto yl = cg(block, alphas, 1);

            // Update the image
            x += stepsize * bp_block(block, yl);
        }

        if (callback) {
            callback(x, k);
        }
    }

    return x;
}

} // namespace reconstruction
} // namespace tomo
