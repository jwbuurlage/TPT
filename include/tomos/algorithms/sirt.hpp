#pragma once

#include <functional>
#include <optional>
#include <vector>

#include "../geometry.hpp"
#include "../projector.hpp"
#include "../util/matrix_sums.hpp"

namespace tomo {
namespace reconstruction {

namespace iterate {

template <dimension D, typename T>
tomo::image<D, T>&
sirt(tomo::image<D, T>& x, const tomo::geometry::base<D, T>& geometry,
     tomo::dim::base<D, T>& kernel, const tomo::projections<D, T>& b,
     tomo::projections<D, T>& rs, tomo::image<D, T>& bcs,
     tomo::projections<D, T>* b1 = nullptr, tomo::image<D, T>* b2 = nullptr) {
    auto v = kernel.get_volume();

    bool clean_up_1 = false;
    bool clean_up_2 = false;
    if (!b1) {
        b1 = new tomo::projections<D, T>(geometry, (T)0);
        clean_up_1 = true;
    } else {
        b1->clear();
    }
    if (!b2) {
        b2 = new tomo::image<D, T>(v, (T)0);
        clean_up_2 = true;
    } else {
        b2->clear();
    }

    // compute Wx
    for (auto[proj, idx, line] : geometry) {
        (void)proj;
        for (auto elem : kernel(line)) {
            (*b1)[idx] += x[elem.index] * elem.value;
        }
    }

    // compute R(p - Wx)
    for (auto j = 0u; j < geometry.lines(); ++j) {
        (*b1)[j] = (b[j] - (*b1)[j]) * rs[j];
    }

    // multiply with W^T
    for (auto[proj, idx, line] : geometry) {
        (void)proj;
        for (auto elem : kernel(line)) {
            (*b2)[elem.index] += elem.value * (*b1)[idx];
        }
    }

    // update image while scaling with beta * C
    for (auto j = 0u; j < v.cells(); ++j) {
        x[j] += bcs[j] * (*b2)[j];
    }

    if (clean_up_1) {
        delete b1;
    }
    if (clean_up_2) {
        delete b2;
    }

    return x;
}

} // namespace iterate

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
                 double beta = 0.5, int iterations = 10,
                 std::function<void(image<D, T>&)> callback = {}) {
    image<D, T> f(v);

    // first we compute R and C
    auto R = tomo::row_sums<D, T>(g, kernel);
    auto bC = tomo::column_sums<D, T>(g, kernel);

    for (auto& r : R) {
        r = (math::abs(r) > math::epsilon<T>) ? ((T)1.0 / r) : (T)0.0;
    }

    for (auto& bc : bC) {
        bc = (math::abs(bc) > math::epsilon<T>) ? ((T)beta / bc) : (T)0.0;
    }

    projections<D, T> s1(g);
    image<D, T> s2(v);
    for (int k = 0; k < iterations; ++k) {
        iterate::sirt(f, g, kernel, p, R, bC, &s1, &s2);

        if (callback) {
            callback(f);
        }
    }

    return f;
}

} // namespace reconstruction
} // namespace tomo
