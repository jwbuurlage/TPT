#pragma once

#include <string>

#include "common.hpp"
#include "image.hpp"
#include "projector.hpp"
#include "projectors/linear.hpp"
#include "sinogram.hpp"
#include "volume.hpp"

namespace tomo {

/**
 * Perform a forward-projection of a given image.
 * FIXME image should be const ref.
 * */
template <dimension D, typename T, class Geometry, class Image, class Projector>
sinogram<D, T, Geometry, Projector>
forward_projection(Image& f, const Geometry& g, Projector& proj) {
    auto sino = sinogram<D, T, Geometry, Projector>(g);

    int line_number = 0;
    for (auto line : g) {
        for (auto elem : proj(line)) {
            sino[line_number] += f[elem.index] * elem.value;
        }
        ++line_number;
    }

    return sino;
}

/** Perform a back-projection of a given sinogram. */
template <dimension D, typename T, typename Geometry, class Projector>
image<D, T> back_projection(const sinogram<D, T, Geometry, Projector>& sino,
                                const Geometry& g, Projector& proj, volume<D> v) {
    auto f = image<D, T>(v);

    int line_number = 0;
    for (auto line : g) {
        for (auto elem : proj(line)) {
            f[elem.index] += sino[line_number] * elem.value;
        }
        ++line_number;
    }

    return f;
}

} // namespace tomo
