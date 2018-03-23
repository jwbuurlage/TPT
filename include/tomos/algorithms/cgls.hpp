#pragma once

#include <functional>
#include <optional>
#include <vector>

#include "../geometry.hpp"
#include "../projector.hpp"

namespace tomo {
namespace reconstruction {

// CGLS
// ----
//
// x0 = 0
// d0 = b
// r0 = A^T b
// p0 = r0
// t0 = A p0
//
// for k : 1..
// (1) \alpha_k = ||r_{k - 1} ||^2 / ||t_{k - 1}||^2
// (2) x_k = x_{k - 1} + \alpha_k p_{k - 1}
// (3) d_k = d_{k - 1} - \alpha_k t_{k - 1}
// (4) r_k = A^T d_k
// (5) \beta_k = ||r_k||^2 / ||r_{k - 1}||^2
// (6) p_k = r_k + \beta_k * p_{k - 1}
// (7) t_k = A p_k
//
// -- done --

/**
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
image<D, T> cgls(const volume<D, T>& v, const tomo::geometry::base<D, T>& g,
                 tomo::dim::base<D, T>& kernel, const projections<D, T>& b,
                 int iterations = 10,
                 std::function<void(image<D, T>&, int)> callback = {}) {
    using namespace tomo::img;

    // x0 = 0
    image<D, T> x(v);

    // d0 = b
    auto d = b;

    // r0 = A^T b
    auto r = tomo::back_projection(b, g, kernel, v);

    // p0 = r0
    auto p = r;

    // t0 = A p0

    // for k : 1..
    for (int k = 0; k < iterations; ++k) {
        // (..7) t_k = A p_k
        auto t = tomo::forward_projection(p, g, kernel);

        // (1) \alpha_k = ||r_{k - 1} ||^2 / ||t_{k - 1}||^2
        auto r_norm = math::norm(r);
        auto t_norm = math::norm(t);
        auto alpha = (r_norm * r_norm) / (t_norm * t_norm);

        // (2) x_k = x_{k - 1} + \alpha_k p_{k - 1}
        x = x + alpha * p;

        // (3) d_k = d_{k - 1} - \alpha_k t_{k - 1}
        d -= alpha * t;

        // (4) r_k = A^T d_k
        r = tomo::back_projection(d, g, kernel, v);

        // (5) \beta_k = ||r_k||^2 / ||r_{k - 1}||^2
        auto rk_norm = math::norm(r);
        auto beta = (rk_norm * rk_norm) / (r_norm * r_norm);

        // (6) p_k = r_k + \beta_k * p_{k - 1}
        p = r + beta * p;

        if (callback) {
            callback(x, k);
        }
    }

    return x;
}

} // namespace reconstruction
} // namespace tomo
