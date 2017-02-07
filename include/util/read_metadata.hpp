#pragma once

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>

#include <cpptoml.h>

#include "geometry.hpp"
#include "volume.hpp"
#include "geometries/parallel.hpp"

namespace tomo {

class invalid_geometry_config_error : public std::runtime_error {
    using runtime_error::runtime_error;
};

template <tomo::dimension D, typename T, typename S>
tomo::math::vec<D, T> stdvec_to_tomovec(std::vector<S> in) {
    assert(in.size() == D);

    auto out = tomo::math::vec<D, T>();
    for (int d = 0; d < D; ++d) {
        out[d] = in[d];
    }

    return out;
}

template <tomo::dimension D, typename T>
std::unique_ptr<tomo::geometry::parallel<D, T>>
read_parallel_geometry(std::shared_ptr<cpptoml::table> parameters,
                       tomo::volume<D, T> v) {
    /* TODO:
     * - Extend parallel geometry support
     * - Detector shape is square
     * - Many fields ignored:
     *     tilt, size, position, source position, and so on
     */
    auto angle_count = (int)(*parameters->get_as<int64_t>("projection-count"));
    auto detector_shape = stdvec_to_tomovec<D - 1, int>(
        *parameters->get_array_of<int64_t>("detector-shape"));
    return std::make_unique<tomo::geometry::parallel<D, T>>(v, angle_count,
                                                            detector_shape[0]);
}

template <tomo::dimension D, typename T>
std::unique_ptr<tomo::geometry::base<D, T>>
read_geometry(std::string kind, std::shared_ptr<cpptoml::table> parameters,
              tomo::volume<D, T> v) {
    if (kind == "parallel") {
        return read_parallel_geometry<D, T>(parameters, v);
    } else {
        throw invalid_geometry_config_error(
            "invalid or unsupported 'type' supplied for geometry");
    }
    return nullptr;
}

template <tomo::dimension D, typename T>
tomo::volume<D, T> read_volume(std::shared_ptr<cpptoml::table> parameters,
                               tomo::math::vec<D - 1, int> detector_size) {

    auto voxels =
        tomo::math::extend<D, int>(detector_size, D - 1, detector_size[0]);

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

    return tomo::volume<D, T>(voxels, min_point, max_point - min_point);
}

template <tomo::dimension D, typename T>
std::pair<std::unique_ptr<tomo::geometry::base<D, T>>, tomo::volume<D, T>>
read_configuration(std::string file) {
    using namespace std::string_literals;

    auto config = cpptoml::parse_file(file);

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

    // First read the detector count from the parameter list
    // just read it separately
    auto detector_size_array =
        config->get_qualified_array_of<int64_t>("parameters.detector-shape");
    auto detector_size = stdvec_to_tomovec<D - 1, int>(*detector_size_array);

    auto v = read_volume<D, T>(config->get_table("volume"), detector_size);

    auto kind = config->get_as<std::string>("type");
    auto g = read_geometry<D, T>(*kind, config->get_table("parameters"), v);

    return std::make_pair(std::move(g), v);
}

} // namespace tomo
