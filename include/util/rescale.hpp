#pragma once

#include "reconstruction_problem.hpp"

namespace tomo {

template <tomo::dimension D, typename T>
reconstruction_problem<D, T>
rescale(reconstruction_problem<D, T>& problem,
        math::vec<D - 1, int> scaled_detector_shape,
        math::vec<D, int> scaled_voxels, int scaled_projection_count) {
    // volume
    auto vol = problem.object_volume;
    vol.set_voxels(scaled_voxels);

    auto g = std::move(problem.acquisition_geometry);
    auto old_projection_count = g->projection_count();
    g->set_detector_shape(scaled_detector_shape);
    g->set_projections(scaled_projection_count);

    auto projs = tomo::projections<D, T>(*g);
    auto projection_step = (T)old_projection_count / scaled_projection_count;
    for (int i = 0; i < g->projection_count(); ++i) {
        auto current_projection = (int)(i * projection_step);
        auto proj = problem.projection_stack.get_projection(current_projection);
        auto f = tomo::downscale_(proj, scaled_detector_shape);
        projs.set_projection(i, f);
    }

    // geometry
    // easiest way is 'set detector', take every 'x-th' projection only

    // the easiest way is 'get_projection' -> 'downscale_image' ->
    // 'set_projection'

    return reconstruction_problem<D, T>{std::move(g), vol, std::move(projs)};
}

} // namespace tomo
