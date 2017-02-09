#pragma once

#include <string>

#include "common.hpp"
#include "image.hpp"
#include "projector.hpp"
#include "projectors/linear.hpp"
#include "projections.hpp"
#include "volume.hpp"

namespace tomo {

/**
 * Perform a forward-projection of a given image.
 * TODO: image should be const ref.
 * */
template <dimension D, typename T, class Image>
projections<D, T> forward_projection(Image& f, const geometry::base<D, T>& g,
                                  dim::base<D, T>& proj) {
    auto sino = projections<D, T>(g);

    int line_number = 0;
    for (auto line : g) {
        for (auto elem : proj(line)) {
            sino[line_number] += f[elem.index] * elem.value;
        }
        ++line_number;
    }

    return sino;
}

/** Perform a back-projection of the given projections. */
template <dimension D, typename T>
image<D, T> back_projection(const projections<D, T>& sino, const geometry::base<D, T>& g,
                            dim::base<D, T>& proj, volume<D, T> v) {
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
