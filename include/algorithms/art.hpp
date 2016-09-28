#pragma once

#include "../common.hpp"
#include "../geometry.hpp"
#include "../image.hpp"
#include "../math.hpp"
#include "../sinogram.hpp"
#include "../volume.hpp"

namespace tomo {
namespace reconstruction {

/**
 * The Algebraic Reconstruction Technique (ART), is a tomographic reconstruction
 * method based on the classic Kaczmarz method for inverse problems.
 *
 * \tparam D the dimension of the problem
 * \tparam T the scalar type in use
 * \tparam Geometry the type of geometry of the problem
 * \tparam Projector the discrete integration method to use
 *
 * \param v the volume of the imaged object
 * \param g the geometry of the problem
 * \param p the measurements (sinogram)
 * \param beta (optional) a relaxation parameter
 * \param iterations (optional) the number of iterations to perform
 *
 * \returns An image object representing the reconstructed object.
 */
template <dimension D, typename T, class Geometry, class Projector>
image<D, T> art(const volume<D>& v, const Geometry& g,
                const sinogram<D, T, Geometry, Projector>& p, double beta = 0.5,
                int iterations = 10) {
    image<D, T> f(v);
    Projector proj(v);

    // compute $w_i \cdot w_i$
    std::vector<T> w_norms(g.lines());
    int line_number = 0;
    for (auto line : g) {
        for (auto elem : proj(line)) {
            w_norms[line_number] += elem.value * elem.value;
        }
        ++line_number;
    }

    for (int k = 0; k < iterations; ++k) {
        int row = 0;
        for (auto line : g) {
            T alpha = 0.0;
            for (auto elem : proj(line))
                alpha += f[elem.index] * elem.value;

            auto factor = beta * ((p[row] - alpha) / w_norms[row]);
            for (auto elem : proj)
                f[elem.index] += factor * elem.value;

            ++row;
        }
    }

    return f;
}

} // namespace reconstruction
} // namespace tomo
