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

math::vec<2_D, int> compute_pixel_intersection(math::vec3<T> x,
                                               geometry::projection<3_D, T> p) {
    auto geometric_point = math::ray_plane_intersection<T>(
        p.source_location, x - p.source_location, p.detector_location,
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

shadow compute_shadow(std::vector<math::vec3<T>> xs,
                      geometry::projection<3_D, T> p) {
    shadow s = {
        {std::numeric_limits<int>::max(), std::numeric_limits<int>::max()},
        {-1, -1}};
    for (auto x : xs) {
        auto intersection = compute_pixel_intersection(x, p);
        intersection =
            math::min(intersection, p.detector_shape - math::vec2<int>{1, 1});
        intersection = math::max(intersection, {0, 0});
        s.min_pt = math::min(s.min_pt, intersection);
        s.max_pt = math::max(s.max_pt, intersection);
    }
    // FIXME enlarge by one
    return s;
}

class restricted_geometry : public geometry::base<3, T> {
  public:
    restricted_geometry(geometry::trajectory<3_D, T>& geometry,
                        volume<3_D, T> local_volume)
        : geometry::base<3, T>(geometry.projection_count(), false),
          geometry_(geometry), local_volume_(local_volume),
          pixels_(geometry.projection_count() + 1) {
        project();
        std::transform(shadows_.begin(), shadows_.end(), pixels_.begin() + 1,
                       [&](auto shadow) { return pixels(shadow); });
        std::partial_sum(pixels_.begin(), pixels_.end(), pixels_.begin());
        this->compute_lines_();
    }

    void project() {
        for (int i = 0; i < geometry_.projection_count(); ++i) {
            shadows_.push_back(compute_shadow(local_volume_.corners(),
                                              geometry_.get_projection(i)));
        }
    }

    math::vec<2_D, int> projection_shape(int i) const override {
        return {shadows_[i].max_pt.x - shadows_[i].min_pt.x,
                shadows_[i].max_pt.y - shadows_[i].min_pt.y};
    }

    math::vec<3_D, T> detector_corner(int i) const override {
        return geometry_.detector_corner(i) +
               (T)shadows_[i].min_pt[0] * geometry_.projection_delta(i)[0] +
               (T)shadows_[i].min_pt[1] * geometry_.projection_delta(i)[1];
    }

    math::vec<3_D, T> source_location(int i) const override {
        return geometry_.source_location(i);
    }

    std::array<math::vec<3_D, T>, 2> projection_delta(int i) const override {
        return geometry_.projection_delta(i);
    }

  private:
    geometry::trajectory<3_D, T>& geometry_;
    volume<3_D, T> local_volume_;
    std::vector<shadow> shadows_;
    std::vector<int> pixels_;
};

// Assumes global volume has physical size [0, 1]^3
template <int G>
auto calculate_local_volume(bulk::rectangular_partitioning<3, G>& partitioning,
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
    auto partitioning = bulk::block_partitioning<3, 3>(
        {options.k, options.k, options.k}, {2, 2, 2});

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
        plotter.send_projection_data(local_geometry, local_proj_stack,
                                     local_volume);
    });
}

int main(int argc, char* argv[]) {
    auto options = util::args(argc, argv);
    run(options);

    return 0;
}
