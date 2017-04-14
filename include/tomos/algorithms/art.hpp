#pragma once

#include "../common.hpp"
#include "../geometry.hpp"
#include "../image.hpp"
#include "../math.hpp"
#include "../projections.hpp"
#include "../volume.hpp"

namespace tomo {
namespace reconstruction {

/**
 * The Algebraic Reconstruction Technique (ART), is a tomographic reconstruction
 * method based on the classic Kaczmarz method for inverse problems.
 *
 * \tparam D the dimension of the problem
 * \tparam T the scalar type in use
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
image<D, T> art(const volume<D, T>& v, const tomo::geometry::base<D, T>& g,
                tomo::dim::base<D, T>& kernel, const projections<D, T>& p,
                double beta = 0.5, int iterations = 10) {
    image<D, T> f(v);

    // compute $w_i \cdot w_i$
    std::vector<T> w_norms(g.lines());
    int line_number = 0;
    for (auto line : g) {
        for (auto elem : kernel(line)) {
            w_norms[line_number] += elem.value * elem.value;
        }
        ++line_number;
    }

    for (int k = 0; k < iterations; ++k) {
        int row = 0;
        for (auto line : g) {
            T alpha = 0.0;
            for (auto elem : kernel(line))
                alpha += f[elem.index] * elem.value;

            auto factor = beta * ((p[row] - alpha) / w_norms[row]);
            for (auto elem : kernel)
                f[elem.index] += factor * elem.value;

            ++row;
        }
    }

    return f;
}

} // namespace reconstruction
} // namespace tomo
