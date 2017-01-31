#pragma once

#include <functional>
#include <vector>

namespace tomo {
namespace reconstruction {

/**
 * The Simultaneous Iterative Reconstruction Technique (SIRT), is a tomographic
 * reconstruction method based on the iterative scheme:
 *
 * \f[ \vec{x}_{i + 1} = \vec{x}_i + C W^T R(\vec{p} - W \vec{x}_i). \f]
 *
 * Here:
 * - \f$\vec{x}_i\f$ represents the reconstructed image in the i-th iteration.
 * - \f$W\f$ is the projection matrix defined by the geometry.
 * - \f$C, R\f$ are the inverse column and row sums respectively.
 * - \f$\vec{p}\f$ is the sinogram.
 *
 * \tparam D the dimension of the problem
 * \tparam T the scalar type in use
 *
 * \param v the volume of the imaged object
 * \param g the geometry of the problem
 * \param p the measurements (sinogram)
 * \param beta (optional) a relaxation parameter
 * \param iterations (optional) the number of iterations to perform
 *
 * \returns An image object representing the reconstructed object.
 */
template <dimension D, typename T>
image<D, T> sirt(const volume<D>& v, const tomo::geometry::base<D, T>& g,
                 tomo::dim::base<D, T>& kernel, const sinogram<D, T>& p,
                 double beta = 0.5, int iterations = 10,
                 std::function<void(image<D, T>&)> callback = {}) {
    image<D, T> f(v);

    // first we compute R and C
    std::vector<T> R(g.lines());
    std::vector<T> bC(v.cells());

    int line_idx = 0;
    for (auto line : g) {
        for (auto elem : kernel(line)) {
            R[line_idx] += elem.value;
            bC[elem.index] += elem.value;
        }
        ++line_idx;
    }

    for (auto& r : R)
        r = (math::abs(r) > math::epsilon<T>) ? ((T)1.0 / r) : (T)0.0;
    for (auto& bc : bC)
        bc = (math::abs(bc) > math::epsilon<T>) ? ((T)beta / bc) : (T)0.0;

    sinogram<D, T> s1(g);
    image<D, T> s2(v);
    for (int k = 0; k < iterations; ++k) {
        s1.clear();

        // compute Wx
        int line_number = 0;
        for (auto line : g) {
            for (auto elem : kernel(line)) {
                s1[line_number] += f[elem.index] * elem.value;
            }
            ++line_number;
        }

        // compute R(p - Wx)
        for (int j = 0; j < g.lines(); ++j) {
            s1[j] = (p[j] - s1[j]) * R[j];
        }

        // zero the image
        s2.clear();

        // multiply with W^T
        int row = 0;
        for (auto line : g) {
            for (auto elem : kernel(line)) {
                s2[elem.index] += elem.value * s1[row];
            }
            ++row;
        }

        // update image while scaling with beta * C
        for (int j = 0; j < v.cells(); ++j) {
            f[j] += bC[j] * s2[j];
        }

        if (callback) {
            callback(f);
        }
    }

    return f;
}

} // namespace reconstruction
} // namespace tomo
