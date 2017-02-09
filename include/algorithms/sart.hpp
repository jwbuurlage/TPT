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
                 double beta = 0.5, int iterations = 10) {
    image<D, T> f(v);

    // the size of a single block
    int k = g.groups()[0];

    // compute $w_i \cdot w_i$
    std::vector<T> w_norms(g.lines());
    int line_number = 0;
    for (auto line : g) {
        for (auto elem : kernel(line)) {
            w_norms[line_number] += elem.value * elem.value;
        }
        ++line_number;
    }

    auto f_next = f;

    for (int iter = 0; iter < iterations; ++iter) {
        int s = k;
        int t = 0;

        int row = 0;
        for (auto line : g) {
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

            ++row;
            ++s;
        }
    }

    return f;
}

} // namespace reconstruction
} // namespace tomo
