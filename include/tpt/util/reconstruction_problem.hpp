#pragma once

#include <memory>

#include "../projections.hpp"
#include "../geometry.hpp"
#include "../volume.hpp"
#include "../common.hpp"

namespace tpt {

template <dimension D, typename T>
struct reconstruction_problem {
    std::unique_ptr<geometry::base<D, T>> acquisition_geometry;
    volume<D, T> object_volume;
    std::unique_ptr<projections<D, T>> projection_stack;
};

} // namespace tpt
