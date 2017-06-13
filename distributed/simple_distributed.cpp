#include "bulk/backends/mpi/mpi.hpp"
#include "bulk/bulk.hpp"

namespace bulk {
using namespace experimental;
}

#include "tomos/tomos.hpp"
#include "tomos/util/plotter.hpp"

using namespace tomo;

using T = float;

struct shadow {
    math::vec2<int> min_pt;
    math::vec2<int> max_pt;
};

auto pixels(shadow s) {
    return (s.max_pt.x - s.min_pt.x) * (s.max_pt.y - s.min_pt.y);
}

math::vec<2_D, int> compute_pixel_intersection(math::vec3<T> point,
                                               geometry::projection<3_D, T> p) {
    auto geometric_point = math::ray_plane_intersection<T>(
        p.source_location, point - p.source_location, p.detector_location,
        math::cross<T>(p.detector_tilt[0], p.detector_tilt[1]));
    assert(geometric_point);
    auto pt = geometric_point.value();
    // detector location is the *center* of the detector, so we adjust
    auto offset = pt - (p.detector_location -
                        (T)0.5 * (p.detector_size[0] * p.detector_tilt[0] +
                                  p.detector_size[1] * p.detector_tilt[1]));
    auto dx = math::dot<3_D, T>(offset, p.detector_tilt[0]);
    auto dy = math::dot<3_D, T>(offset, p.detector_tilt[1]);
    return {(int)(dx / p.detector_size[0] * p.detector_shape[0]),
            (int)(dy / p.detector_size[1] * p.detector_shape[1])};
}

shadow compute_shadow(std::vector<math::vec3<T>> points,
                      geometry::projection<3_D, T> p, math::vec2<int> shape) {
    shadow s = {
        {std::numeric_limits<int>::max(), std::numeric_limits<int>::max()},
        {-1, -1}};
    for (auto pt : points) {
        auto intersection = compute_pixel_intersection(pt, p);
        intersection = math::min(intersection, shape - math::vec2<int>{1, 1});
        intersection = math::max(intersection, {0, 0});
        s.min_pt = math::min(s.min_pt, intersection);
        s.max_pt = math::max(s.max_pt, intersection);
    }
    return s;
}

class restricted_geometry : public geometry::base<3, T> {
  public:
    restricted_geometry(geometry::trajectory<3_D, T>& geometry,
                        volume<3_D, T> local_volume)
        : geometry::base<3, T>(geometry.projection_count(), {-1, -1}),
          geometry_(geometry), local_volume_(local_volume),
          pixels_(geometry.projection_count() + 1) {
        project();
        std::transform(shadows_.begin(), shadows_.end(), pixels_.begin() + 1,
                       [&](auto shadow) { return pixels(shadow); });
        std::partial_sum(pixels_.begin(), pixels_.end(), pixels_.begin());
        this->set_line_count(pixels_[pixels_.size() - 1]);
    }

    void project() {
        for (int i = 0; i < geometry_.projection_count(); ++i) {
            shadows_.push_back(compute_shadow(local_volume_.corners(),
                                              geometry_.get_projection(i),
                                              geometry_.detector_shape()));
        }
    }

    virtual math::ray<3, T> get_line(int i) const override {
        // 1 compute projection (FIXME how to identify in O(1)?, need to rethink
        // the whole 'get_line' system, making the iterators locally can boost
        // performance greatly)
        auto p = (std::find_if(pixels_.begin(), pixels_.end(),
                               [=](auto x) { return x > i; }) -
                  pixels_.begin()) -
                 1;
        assert(p >= 0);
        assert(p < geometry_.projection_count());
        // 2 compute global index in projection
        auto j = i - pixels_[p];
        auto x = j % (shadows_[p].max_pt.x - shadows_[p].min_pt.x);
        auto y = j / (shadows_[p].max_pt.x - shadows_[p].min_pt.x);
        x += shadows_[p].min_pt.x;
        y += shadows_[p].min_pt.y;
        int global_index = p * geometry_.detector_pixel_count() +
                           y * geometry_.detector_shape()[0] + x;
        assert(global_index >= 0 && global_index < geometry_.lines());
        // 3 return get_line of global index
        return geometry_.get_line(global_index);
    }

  private:
    geometry::trajectory<3_D, T>& geometry_;
    volume<3_D, T> local_volume_;
    std::vector<shadow> shadows_;
    std::vector<int> pixels_;
};

// Assumes global volume has physical size [0, 1]^3
auto calculate_local_volume(bulk::rectangular_partitioning<3, 1>& partitioning,
                            int s) {
    auto to_vec = [](auto array_like) -> math::vec3<T> {
        math::vec3<T> result;
        for (int d = 0; d < 3; ++d) {
            result[d] = array_like[d];
        }
        return result;
    };

    auto voxel_origin = to_vec(partitioning.origin(s));
    auto voxel_size = to_vec(partitioning.local_size(s));
    auto global_voxels = to_vec(partitioning.global_size());

    auto relative_origin = math::vec3<T>(voxel_origin / global_voxels);
    auto relative_size = math::vec3<T>(voxel_size / global_voxels);

    return volume<3_D, T>(voxel_size, relative_origin, relative_size);
}

void run(util::args options) {
    bulk::mpi::environment env;

    // this creates a volume [0, 1]^3, with k voxels in each axis
    auto global_volume = volume<3_D, T>(options.k);
    auto global_geometry =
        geometry::cone_beam<T>(global_volume, options.k, {1.5, 1.5},
                               {options.k, options.k}, 10.0, 2.0);

    auto processors = env.available_processors();
    auto partitioning = bulk::block_partitioning<3, 1>(
        {options.k, options.k, options.k}, {processors}, {1});

    env.spawn(processors, [&](auto& world) {
        auto local_volume =
            calculate_local_volume(partitioning, world.processor_id());
        (void)world;

        util::ext_plotter<3_D, T> plotter(
            "tcp://localhost:5555", "Distributed restricted geometry test");

        image<3_D, T> phantom(local_volume);
        fill_ellipsoids_(phantom, mshl_ellipsoids_<T>(), local_volume,
                         global_volume);

        auto local_geometry =
            restricted_geometry(global_geometry, local_volume);
        auto local_proj_stack = projections<3_D, T>(local_geometry);

        // This is now a normal FP
        auto kernel = dim::joseph<3_D, T>(local_volume);
        int line_number = 0;
        for (auto line : local_geometry) {
            for (auto elem : kernel(line)) {
                local_proj_stack[line_number] +=
                    phantom[elem.index] * elem.value;
            }
            ++line_number;
        }

        /* // we need to compute overlaps and so on.. TODO
        auto sinos = bulk::gather_all(rsino.data()); */

        plotter.plot(phantom);
    });
}

int main(int argc, char* argv[]) {
    auto options = util::args(argc, argv);
    run(options);

    return 0;
}
