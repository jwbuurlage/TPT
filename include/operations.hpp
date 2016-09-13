#pragma once

#include <string>

#include "common.hpp"
#include "image.hpp"
#include "projector.hpp"
#include "projectors/linear.hpp"
#include "sinogram.hpp"
#include "volume.hpp"

namespace tomo {

/** Perform a forward-projection of a given image. */
template <dimension D, typename T, class Geometry,
          class Projector = linear_projector<D, T>>
sinogram<D, T, Geometry, Projector>
forward_projection(const image<D, T>& f, const Geometry& g, Projector& proj) {
    auto sino = sinogram<D, T, Geometry, Projector>(g);

    int line_number = 0;
    for (auto line : g) {
        proj.reset(line);
        for (auto elem : proj) {
            sino[line_number] += f[elem.index] * elem.value;
        }
        ++line_number;
    }

    return sino;
}

/** Perform a back-projection of a given sinogram. */
template <dimension D, typename T, typename Geometry, class Projector>
image<D, T> backward_projection(const sinogram<D, T, Geometry, Projector>& p,
                                const Geometry& g);

} // namespace tomo
