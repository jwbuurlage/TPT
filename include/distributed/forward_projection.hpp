#pragma once

#include <string>
#include <vector>

#include "common.hpp"
#include "image.hpp"
#include "partitioned_sinogram.hpp"
#include "partitioned_image.hpp"
#include "projector.hpp"
#include "projectors/linear.hpp"
#include "sinogram.hpp"
#include "volume.hpp"

#include "bulk/bulk.hpp"

namespace tomo {
namespace distributed {

/**
 * Perform a forward-projection of a given image.
 *
 * TODO: alternatively, we can separate the 'geometry communication' info, and
 * the partitioned sinogram itself, so that we can still return it from this
 * function.
 * */
template <dimension D, typename T, class Geometry, class Projector>
void forward_project(tomo::distributed::partitioned_image<D, T>& f,
                     const Geometry& g, Projector& proj,
                     partitioned_sinogram<D, T, Geometry>& sino) {
    int line_number = 0;
    for (auto line : g) {
        for (auto elem : proj(line)) {
            sino[line_number] += f[elem.index] * elem.value;
        }
        ++line_number;
    }

    sino.harmonize();
}

/**
 * Perform a back-projection of a given image.
 *
 * TODO: alternatively, we can separate the 'geometry communication' info, and
 * the partitioned sinogram itself, so that we can still return it from this
 * function.
 * */
template <dimension D, typename T, class Geometry, class Projector>
void back_project(tomo::distributed::partitioned_image<D, T>& f,
                  const Geometry& g, Projector& proj,
                  partitioned_sinogram<D, T, Geometry>& sino) {
    int line_number = 0;
    for (auto line : g) {
        for (auto elem : proj(line)) {
            f[elem.index] += sino[line_number] * elem.value;
        }
        ++line_number;
    }
}

} // namespace distributed
} // namespace tomo
