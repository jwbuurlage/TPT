#pragma once

#include <memory>

#include "../projections.hpp"
#include "../geometry.hpp"
#include "../volume.hpp"
#include "../common.hpp"

namespace tomo {

template <tomo::dimension D, typename T>
struct reconstruction_problem {
    std::unique_ptr<tomo::geometry::base<D, T>> acquisition_geometry;
    tomo::volume<D, T> object_volume;
    tomo::projections<D, T> projection_stack;
};

} // namespace tomo
