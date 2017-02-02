#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>

#include <cpptoml.h>

#include "tomo.hpp"

class invalid_geometry_config_error : public std::runtime_error {
    using runtime_error::runtime_error;
};

template <tomo::dimension D, typename T>
std::unique_ptr<tomo::geometry::parallel<D, T>>
read_parallel_geometry(std::shared_ptr<cpptoml::table> parameters,
                       tomo::volume<D, T> v) {
    (void)parameters;
    return std::make_unique<tomo::geometry::parallel<D, T>>(v, 10, 10);
}

template <tomo::dimension D, typename T>
std::unique_ptr<tomo::geometry::base<D, T>>
read_geometry(std::string kind, std::shared_ptr<cpptoml::table> parameters,
              tomo::volume<D, T> v) {
    if (kind == "parallel") {
        read_parallel_geometry<D, T>(parameters, v);
    } else {
        throw invalid_geometry_config_error(
            "invalid or unsupported 'kind' supplied for geometry");
    }
    return nullptr;
}

template <tomo::dimension D, typename T>
tomo::volume<D, T> read_volume(std::shared_ptr<cpptoml::table> parameters) {
    // cpptoml only supports double precision and 64 bit integer parsing,
    // we cast down later
    auto origin_array = parameters->get_array_of<double>("origin");
    auto lengths_array = parameters->get_array_of<double>("lengths");
    auto voxels_array = parameters->get_array_of<int64_t>("voxels");
    if ((*origin_array).size() != D || (*lengths_array).size() != D ||
        (*voxels_array).size() != D) {
        throw invalid_geometry_config_error(
            "required field missing or invalid for volume, please supply "
            "appropriately sized arrays for 'origin', 'lengths' and 'voxels'");
    }
    return tomo::volume<D, T>(10);
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

    auto dimension = config->get_as<int64_t>("dimension");
    if (*dimension != D) {
        throw invalid_geometry_config_error("dimension of the specification "
                                            "does not match the required "
                                            "dimension");
    }

    auto v = read_volume<D, T>(config->get_table("volume"));

    auto kind = config->get_as<std::string>("kind");
    auto g = read_geometry<D, T>(*kind, config->get_table("parameters"), v);

    return std::make_pair(std::move(g), v);
}

int main() {
    try {
        auto geom_and_volume =
            read_configuration<3_D, float>("data/geometries/parallel.toml");
    } catch (const std::exception& e) {
        std::cout << "Reading configuration failed: " << e.what() << "\n";
    }

    return 0;
}
