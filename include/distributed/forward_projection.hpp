#pragma once

#include <string>

#include "common.hpp"
#include "image.hpp"
#include "projector.hpp"
#include "projectors/linear.hpp"
#include "sinogram.hpp"
#include "volume.hpp"

#include "bulk/bulk.hpp"

namespace tomo {
namespace distributed {

template <tomo::dimension D>
class sinogram_exchanges {
  public:
    struct exchange {
        int global_line_number;
        int owner;
        int remote_index;
    };

    sinogram_exchanges(bulk::partitioning<D, 1>& part, int s) {
        // perform a forward projection and mark all the non-local partitionings
    }

  private:
};

template <dimension D, typename T, class Geometry, class Image, class Projector,
          class World>
class partitioned_sinogram {
  public:
    partitioned_sinogram(Geometry& geometry, sinogram_exchanges<D>& exchanges)
        : exchanges_(exchanges) {}

    void harmonize() {
        for (auto exchange : exchanges_) {
            //..
        }
    }

  private:
    sinogram_exchanges<D> exchanges_;
    bulk::coarray<T, World> data_;
};

/**
 * Perform a forward-projection of a given image.
 *
 * TODO: alternatively, we can separate the 'geometry communication' info, and
 * the partitioned sinogram itself, so that we can still return it from this
 * function.
 * */
template <dimension D, typename T, class Geometry, class Image, class World,
          class Projector = dim::linear<D, T>>
auto forward_project(
    Image& f, const Geometry& g, Projector& proj,
    partitioned_sinogram<D, T, Geometry, Image, Projector, World> sino) {
    int line_number = 0;
    for (auto line : g) {
        for (auto elem : proj(line)) {
            sino[line_number] += f[elem.index] * elem.value;
        }
        ++line_number;
    }

    return sino;
}

} // namespace distributed
} // namespace tomo
