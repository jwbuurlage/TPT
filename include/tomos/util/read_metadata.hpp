#pragma once

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>

#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

#include <cpptoml.h>

#include "../common.hpp"
#include "../geometries/cone.hpp"
#include "../geometries/dual_axis_parallel.hpp"
#include "../geometries/helical_cone_beam.hpp"
#include "../geometries/laminography.hpp"
#include "../geometries/parallel.hpp"
#include "../geometries/tomosynthesis.hpp"
#include "../geometry.hpp"
#include "../projections.hpp"
#include "../volume.hpp"
#include "read_tiff.hpp"
#include "reconstruction_problem.hpp"

namespace tomo {

class invalid_geometry_config_error : public std::runtime_error {
    using runtime_error::runtime_error;
};

// ---------------------------------
// HELPER FUNCTIONS FOR READING VECS

template <tomo::dimension D, typename T, typename S>
tomo::math::vec<D, T> stdvec_to_tomovec(std::vector<S> in) {
    assert(in.size() == D);

    auto out = tomo::math::vec<D, T>();
    for (int d = 0; d < D; ++d) {
        out[d] = in[d];
    }

    return out;
}

template <tomo::dimension D, typename T, typename R>
tomo::math::vec<D, T> read_vec_as_(std::shared_ptr<cpptoml::table> parameters,
                                   std::string name) {
    return stdvec_to_tomovec<D, T>(*parameters->get_array_of<R>(name));
}

template <tomo::dimension D, typename T>
tomo::math::vec<D, T> read_vec_(std::shared_ptr<cpptoml::table> parameters,
                                std::string name, std::true_type) {
    return read_vec_as_<D, T, int64_t>(parameters, name);
}

template <tomo::dimension D, typename T>
tomo::math::vec<D, T> read_vec_(std::shared_ptr<cpptoml::table> parameters,
                                std::string name, std::false_type) {
    return read_vec_as_<D, T, double>(parameters, name);
}

template <tomo::dimension D, typename T>
tomo::math::vec<D, T> read_vec(std::shared_ptr<cpptoml::table> parameters,
                               std::string name) {
    return read_vec_<D, T>(parameters, name, std::is_integral<T>());
}

template <tomo::dimension D, typename T, typename R>
std::vector<tomo::math::vec<D, T>>
read_vecs_as_(std::shared_ptr<cpptoml::table> parameters, std::string name) {
    std::vector<tomo::math::vec<D, T>> result;
    auto entries = parameters->get_array_of<cpptoml::array>(name);
    for (auto entry : *entries) {
        result.push_back(stdvec_to_tomovec<D, T>(*entry->get_array_of<R>()));
    }
    return result;
}

template <tomo::dimension D, typename T>
std::vector<tomo::math::vec<D, T>>
read_vecs_(std::shared_ptr<cpptoml::table> parameters, std::string name,
           std::true_type) {
    return read_vecs_as_<D, T, int64_t>(parameters, name);
}

template <tomo::dimension D, typename T>
std::vector<tomo::math::vec<D, T>>
read_vecs_(std::shared_ptr<cpptoml::table> parameters, std::string name,
           std::false_type) {
    return read_vecs_as_<D, T, double>(parameters, name);
}

template <tomo::dimension D, typename T>
std::vector<tomo::math::vec<D, T>>
read_vecs(std::shared_ptr<cpptoml::table> parameters, std::string name) {
    return read_vecs_<D, T>(parameters, name, std::is_integral<T>());
}

// ---------------------------------

template <typename T>
std::unique_ptr<tomo::geometry::cone_beam<T>>
read_circular_cone_beam_geometry(std::shared_ptr<cpptoml::table> parameters,
                                 tomo::volume<3_D, T> v, int k = -1) {
    auto projection_count =
        (int)(*parameters->get_as<int64_t>("projection-count"));

    auto detector_size = read_vec<2_D, T>(parameters, "detector-size");
    auto detector_shape = read_vec<2_D, int>(parameters, "detector-shape");
    auto source_position = read_vec<3_D, T>(parameters, "source-position");
    auto detector_position = read_vec<3_D, T>(parameters, "detector-position");

    auto detector_tilt_vec = read_vecs<3_D, T>(parameters, "detector-tilt");
    if (detector_tilt_vec.size() != 2) {
        throw invalid_geometry_config_error(
            "detector tilt does not consist of two vectors");
    }

    if (k > 0) {
        detector_shape = {k, k};
        projection_count = k;
    }

    return std::make_unique<tomo::geometry::cone_beam<T>>(
        v, projection_count, detector_size, detector_shape, source_position,
        detector_position, std::array<tomo::math::vec<3_D, T>, 2>{
                               detector_tilt_vec[0], detector_tilt_vec[1]});
}

template <typename T>
std::unique_ptr<tomo::geometry::helical_cone_beam<T>>
read_helical_cone_beam_geometry(std::shared_ptr<cpptoml::table> parameters,
                                tomo::volume<3_D, T> v, int k = -1) {
    auto projection_count =
        (int)(*parameters->get_as<int64_t>("projection-count"));

    auto detector_size = read_vec<2_D, T>(parameters, "detector-size");
    auto detector_shape = read_vec<2_D, int>(parameters, "detector-shape");
    auto source_position = read_vec<3_D, T>(parameters, "source-position");
    auto detector_position = read_vec<3_D, T>(parameters, "detector-position");
    auto rotations = (T)(*parameters->get_as<double>("rotations"));

    auto detector_tilt_vec = read_vecs<3_D, T>(parameters, "detector-tilt");
    if (detector_tilt_vec.size() != 2) {
        throw invalid_geometry_config_error(
            "detector tilt does not consist of two vectors");
    }

    if (k > 0) {
        detector_shape = {k, k};
        projection_count = k;
    }

    return std::make_unique<tomo::geometry::helical_cone_beam<T>>(
        v, projection_count, detector_size, detector_shape, source_position,
        detector_position, rotations,
        std::array<tomo::math::vec<3_D, T>, 2>{detector_tilt_vec[0],
                                               detector_tilt_vec[1]});
}

template <typename T>
std::unique_ptr<tomo::geometry::cone_beam<T>>
read_dynamic_cone_beam_geometry(std::shared_ptr<cpptoml::table> parameters,
                                tomo::volume<3_D, T> v, int k = -1) {
    return nullptr;
}

template <typename T>
std::unique_ptr<tomo::geometry::laminography<T>>
read_laminography_geometry(std::shared_ptr<cpptoml::table> parameters,
                           tomo::volume<3_D, T> v, int k = -1) {
    auto projection_count =
        (int)(*parameters->get_as<int64_t>("projection-count"));

    auto detector_size = read_vec<2_D, T>(parameters, "detector-size");
    auto detector_shape = read_vec<2_D, int>(parameters, "detector-shape");
    auto source_position = read_vec<3_D, T>(parameters, "source-position");
    auto detector_position = read_vec<3_D, T>(parameters, "detector-position");
    auto source_radius = (T)(*parameters->get_as<double>("source-radius"));
    auto detector_radius = (T)(*parameters->get_as<double>("detector-radius"));

    if (k > 0) {
        detector_shape = {k, k};
        projection_count = k;
    }

    return std::make_unique<tomo::geometry::laminography<T>>(
        v, projection_count, detector_size, detector_shape, source_position,
        detector_position, source_radius, detector_radius);
}

template <typename T>
std::unique_ptr<tomo::geometry::tomosynthesis<T>>
read_tomosynthesis_geometry(std::shared_ptr<cpptoml::table> parameters,
                            tomo::volume<3_D, T> v, int k = -1) {

    auto projection_count =
        (int)(*parameters->get_as<int64_t>("projection-count"));

    auto detector_size = read_vec<2_D, T>(parameters, "detector-size");
    auto detector_shape = read_vec<2_D, int>(parameters, "detector-shape");
    auto source_position = read_vec<3_D, T>(parameters, "source-position");
    auto detector_position = read_vec<3_D, T>(parameters, "detector-position");
    auto source_arc = (T)(*parameters->get_as<double>("source-arc"));

    if (k > 0) {
        detector_shape = {k, k};
        projection_count = k;
    }

    return std::make_unique<tomo::geometry::tomosynthesis<T>>(
        v, projection_count, detector_size, detector_shape, source_position,
        detector_position, source_arc);
    return nullptr;
}

template <typename T>
std::unique_ptr<tomo::geometry::cone_beam<T>>
read_double_helix_geometry(std::shared_ptr<cpptoml::table> parameters,
                           tomo::volume<3_D, T> v, int k = -1) {
    return nullptr;
}

template <tomo::dimension D, typename T>
std::unique_ptr<tomo::geometry::parallel<D, T>>
read_parallel_geometry(std::shared_ptr<cpptoml::table> parameters,
                       tomo::volume<D, T> v, int k = -1) {
    /* TODO:
     * - Extend parallel geometry support
     * - Detector shape is square
     * - Many fields ignored:
     *     tilt, size, position, source position, detector shape, and so on
     */
    auto angle_count = (int)(*parameters->get_as<int64_t>("projection-count"));

    if (k > 0) {
        angle_count = k;
        v.set_voxels({k, k, k});
    }
    return std::make_unique<tomo::geometry::parallel<D, T>>(v, angle_count);
}

template <tomo::dimension D, typename T>
std::unique_ptr<tomo::geometry::dual_axis_parallel<D, T>>
read_dual_parallel_geometry(std::shared_ptr<cpptoml::table> parameters,
                            tomo::volume<D, T> v, int k = -1) {
    /* TODO:
     * - Extend parallel geometry support
     * - Detector shape is square
     * - Many fields ignored:
     *     tilt, size, position, source position, detector shape, and so on
     */
    auto angle_count = (int)(*parameters->get_as<int64_t>("projection-count"));
    if (k > 0) {
        angle_count = k / 2;
        v.set_voxels({k, k, k});
    }

    return std::make_unique<tomo::geometry::dual_axis_parallel<D, T>>(
        v, angle_count);
}

template <tomo::dimension D, typename T>
std::unique_ptr<tomo::geometry::base<D, T>>
read_geometry(std::string kind, std::shared_ptr<cpptoml::table> parameters,
              tomo::volume<D, T> v, int k = -1) {
    if (kind == "parallel") {
        return read_parallel_geometry<D, T>(parameters, v, k);
    } else if (kind == "dual-parallel") {
        return read_dual_parallel_geometry<D, T>(parameters, v, k);
    } else if (D != 3_D) {
        throw invalid_geometry_config_error("Only parallel available in 2D");
    } else if (kind == "circular-cone-beam") {
        return read_circular_cone_beam_geometry<T>(parameters, v, k);
    } else if (kind == "helical-cone-beam") {
        return read_helical_cone_beam_geometry<T>(parameters, v, k);
    } else if (kind == "laminography") {
        return read_laminography_geometry<T>(parameters, v, k);
    } else if (kind == "tomosynthesis") {
        return read_tomosynthesis_geometry<T>(parameters, v, k);
    } else {
        throw invalid_geometry_config_error(
            "Invalid or unsupported 'type' supplied for geometry");
    }
    return nullptr;
}

template <tomo::dimension D, typename T>
tomo::volume<D, T> read_volume(std::shared_ptr<cpptoml::table> parameters, int k = -1) {
    auto voxels = read_vec<3_D, int>(parameters, "voxels");

    // cpptoml only supports double precision and 64 bit integer parsing,
    // we cast down later
    auto min_array = parameters->get_array_of<double>("min-point");
    auto max_array = parameters->get_array_of<double>("max-point");
    if ((*min_array).size() != D || (*max_array).size() != D) {
        throw invalid_geometry_config_error(
            "required field missing or invalid for volume, please supply "
            "appropriately sized arrays for 'min-point', 'max-point'");
    }

    auto min_point = stdvec_to_tomovec<D, T>(*min_array);
    auto max_point = stdvec_to_tomovec<D, T>(*max_array);

    if (k > 0) {
        voxels = {k, k, k};
    }

    return tomo::volume<D, T>(voxels, min_point, max_point - min_point);
}

template <tomo::dimension D, typename T>
std::unique_ptr<tomo::projections<D, T>>
read_projection_stack(std::shared_ptr<cpptoml::table> parameters,
                      tomo::geometry::base<D, T>& g, fs::path root_directory) {
	return nullptr;
    try {
        auto proj_count = *parameters->get_qualified_as<int64_t>(
            "parameters.projection-count");
        auto filename_pattern =
            *parameters->get_as<std::string>("projection-filename-pattern");
        return std::make_unique<tomo::projections<D, T>>(
            std::move(tomo::tiff_stack_to_projections<D, T>(
                g, filename_pattern, proj_count, root_directory)));
    } catch (const std::exception& e) {
        return nullptr;
    }
}

template <tomo::dimension D, typename T>
reconstruction_problem<D, T> read_configuration(std::string file, int k = -1) {
    using namespace std::string_literals;

    auto config = cpptoml::parse_file(file);

    auto root_directory = fs::path(file).remove_filename();

    auto specifies = config->get_as<std::string>("specifies");
    if (*specifies != "geometry"s) {
        throw invalid_geometry_config_error(
            "'specifies' field is not equal to 'geometry'");
    }

    auto dim = config->get_as<int64_t>("dimension");
    if (*dim != D) {
        throw invalid_geometry_config_error("dimension of the specification "
                                            "does not match the required "
                                            "dimension");
    }

    auto v = read_volume<D, T>(config->get_table("volume"), k);

    auto kind = config->get_as<std::string>("type");
    auto g = read_geometry<D, T>(*kind, config->get_table("parameters"), v, k);

    assert(g->lines() > 0);

    auto projs = read_projection_stack<D, T>(config, *g, root_directory);

    return reconstruction_problem<D, T>{std::move(g), v, std::move(projs)};
}

} // namespace tomo
