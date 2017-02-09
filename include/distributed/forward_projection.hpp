#pragma once

#include <string>
#include <vector>

#include "common.hpp"
#include "image.hpp"
#include "partitioned_projection_stack.hpp"
#include "partitioned_image.hpp"
#include "projector.hpp"
#include "projectors/linear.hpp"
#include "projections.hpp"
#include "volume.hpp"

#include "bulk/bulk.hpp"

namespace tomo {
namespace distributed {

/**
 * Perform a forward-projection of a given image.
 *
 * TODO: alternatively, we can separate the 'geometry communication' info, and
 * the partitioned projection stack itself, so that we can still return it from this
 * function.
 * */
template <dimension D, dimension G, typename T>
void forward_project(tomo::distributed::partitioned_image<D, G, T>& f,
                     const tomo::geometry::base<D, T>& g, tomo::dim::base<D, T>& proj,
                     partitioned_projection_stack<D, T>& ps) {
    int line_number = 0;
    for (auto line : g) {
        for (auto elem : proj(line)) {
            ps[line_number] += f[elem.index] * elem.value;
        }
        ++line_number;
    }

    ps.harmonize();
}

/**
 * Perform a back-projection of a given image.
 *
 * TODO: alternatively, we can separate the 'geometry communication' info, and
 * the partitioned projection_stack itself, so that we can still return it from this
 * function.
 * */
template <dimension D, dimension G, typename T>
void back_project(tomo::distributed::partitioned_image<D, G, T>& f,
                  const tomo::geometry::base<D, T>& g, tomo::dim::base<D, T>& proj,
                  partitioned_projection_stack<D, T>& ps) {
    int line_number = 0;
    for (auto line : g) {
        for (auto elem : proj(line)) {
            f[elem.index] += ps[line_number] * elem.value;
        }
        ++line_number;
    }
}

} // namespace distributed
} // namespace tomo
