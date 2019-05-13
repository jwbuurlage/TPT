#pragma once

#include "../common.hpp"
#include "../geometry.hpp"
#include "../image.hpp"
#include "../math.hpp"
#include "../projections.hpp"
#include "../volume.hpp"

namespace tpt {
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
image<D, T> art(const volume<D, T>& v, const geometry::base<D, T>& g,
                dim::base<D, T>& kernel, const projections<D, T>& p,
                double beta = 0.5, int iterations = 10,
                std::function<void(image<D, T>&, int)> callback = {}) {
    image<D, T> f(v);

    // compute $w_i \cdot w_i$
    std::vector<T> w_norms(g.lines());
    for (auto[line_number, line] : g) {
        for (auto elem : kernel(line)) {
            w_norms[line_number] += elem.value * elem.value;
        }
    }

    for (int k = 0; k < iterations; ++k) {
        for (auto[row, line] : g) {
            T alpha = 0.0;
            for (auto elem : kernel(line))
                alpha += f[elem.index] * elem.value;

            auto factor = beta * ((p[row] - alpha) / w_norms[row]);
            for (auto elem : kernel)
                f[elem.index] += factor * elem.value;
        }

        if (callback) {
            callback(f, k); 
        }
    }

    return f;
}

} // namespace reconstruction
} // namespace tpt
