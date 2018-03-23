#pragma once

#include <functional>
#include <optional>
#include <vector>

#include "../geometry.hpp"
#include "../projector.hpp"
#include "../util/matrix_sums.hpp"

namespace tomo {
namespace reconstruction {

/**
 * The Simultaneous Iterative Reconstruction Technique (SIRT), is a
 * tomographic
 * reconstruction method based on the iterative scheme:
 *
 * \f[ \vec{x}_{i + 1} = \vec{x}_i + C W^T R(\vec{p} - W \vec{x}_i). \f]
 *
 * Here:
 * - \f$\vec{x}_i\f$ represents the reconstructed image in the i-th
 * iteration.
 * - \f$W\f$ is the projection matrix defined by the geometry.
 * - \f$C, R\f$ are the inverse column and row sums respectively.
 * - \f$\vec{p}\f$ is the projections.
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
image<D, T> sirt(const volume<D, T>& v, const tomo::geometry::base<D, T>& g,
                 tomo::dim::base<D, T>& kernel, const projections<D, T>& p,
                 double beta = 1.0, int iterations = 10,
                 std::function<void(image<D, T>&, int)> callback = {},
                 bool box_constraint = false, T box_min = -1, T box_max = 1) {
    image<D, T> f(v);

    // first we compute R and C
    auto rs = tomo::row_sums<D, T>(g, kernel);
    auto bcs = tomo::column_sums<D, T>(g, kernel);

    for (auto& r : rs) {
        r = (math::abs(r) > math::epsilon<T>) ? ((T)1.0 / r) : (T)0.0;
    }

    for (auto& bc : bcs) {
        bc = (math::abs(bc) > math::epsilon<T>) ? ((T)beta / bc) : (T)0.0;
    }

    projections<D, T> s1(g);
    image<D, T> s2(v);
    for (int k = 0; k < iterations; ++k) {
        // compute Wx
        for (auto[idx, line] : g) {
            for (auto elem : kernel(line)) {
                s1[idx] += f[elem.index] * elem.value;
            }
        }

        // compute R(p - Wx)
        for (auto j = 0u; j < g.lines(); ++j) {
            s1[j] = (p[j] - s1[j]) * rs[j];
        }

        // multiply with W^T
        for (auto[idx, line] : g) {
            for (auto elem : kernel(line)) {
                s2[elem.index] += elem.value * s1[idx];
            }
        }

        // update image while scaling with beta * C
        for (auto j = 0u; j < v.cells(); ++j) {
            f[j] += bcs[j] * s2[j];
        }

        s1.clear();
        s2.clear();

        if (box_constraint) {
            math::box(f, box_min, box_max);
        }

        if (callback) {
            callback(f, k);
        }
    }

    return f;
}

template <dimension D, typename T>
image<D, T> landweber(const volume<D, T>& v, const tomo::geometry::base<D, T>& g,
                 tomo::dim::base<D, T>& kernel, const projections<D, T>& p,
                 double beta = 1.0, int iterations = 10,
                 std::function<void(image<D, T>&, int)> callback = {},
                 bool box_constraint = false, T box_min = -1, T box_max = 1) {
    image<D, T> f(v);

    projections<D, T> s1(g);
    image<D, T> s2(v);
    for (int k = 0; k < iterations; ++k) {
        // compute Wx
        for (auto[idx, line] : g) {
            for (auto elem : kernel(line)) {
                s1[idx] += f[elem.index] * elem.value;
            }
        }

        // compute p - Wx
        for (auto j = 0u; j < g.lines(); ++j) {
            s1[j] = p[j] - s1[j];
        }

        // multiply with W^T
        for (auto[idx, line] : g) {
            for (auto elem : kernel(line)) {
                s2[elem.index] += elem.value * s1[idx];
            }
        }

        // update image
        for (auto j = 0u; j < v.cells(); ++j) {
            f[j] += s2[j] * beta;
        }

        s1.clear();
        s2.clear();

        if (box_constraint) {
            math::box(f, box_min, box_max);
        }

        if (callback) {
            callback(f, k);
        }
    }

    return f;
}



} // namespace reconstruction
} // namespace tomo
