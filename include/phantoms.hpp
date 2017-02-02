#pragma once

#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "common.hpp"
#include "image.hpp"
#include "math.hpp"
#include "volume.hpp"

namespace tomo {

using namespace tomo::literals;

/**
 * An ellipse is defined completely by 6 reals, these are:
 * [I, a, b, x0, y0, phi].
 */
template <typename T>
using ellipse = std::array<T, 6>;

/**
 * An (3D) ellipsoid is defined completely by 10 reals, these are:
 * [I, a, b, c, x0, y0, z0, phi_1, phi_2, phi_3].
 */
template <typename T>
using ellipsoid = std::array<T, 10>;

/** Load a phantom from file. */
// template <dimension D, typename T>
// image<D, T> phantom(std::string filename) {}

template <typename ImageLike, typename T>
void fill_ellipses_(ImageLike& image, const std::vector<ellipse<T>>& ellipses,
                    volume<2_D, T> image_volume,
                    volume<2_D, T> phantom_volume) {
    for (auto e : ellipses) {
        auto& intensity = e[0];
        auto aa = e[1] * e[1];
        auto bb = e[2] * e[2];
        auto& x0 = e[3];
        auto& y0 = e[4];
        auto phi = e[5] * M_PI / (T)180.0;

        auto cos_p = cos(phi);
        auto sin_p = sin(phi);

        auto origin = image_volume.origin();
        for (int i = 0; i < image_volume.voxels()[0]; ++i) {
            for (int j = 0; j < image_volume.voxels()[1]; ++j) {
                // we normalize the vector so that it lies in [0,1]
                math::vec2<T> idx{i, j};
                // component-wise division
                idx /= image_volume.voxels();

                // we see where it lies with respect to the phantom, which
                // depends on:
                // image size, phantom size, image origin, phantom origin
                math::vec2<T> v_x = (origin + idx * image_volume.physical_lengths());
                v_x -= phantom_volume.origin();
                v_x /= phantom_volume.physical_lengths();

                // v_x is now relative in the subimage of the volume

                // compute the intensity increase
                auto x = ((T)2.0 * v_x[0] - (T)1.0) - x0;
                auto y = ((T)2.0 * v_x[1] - (T)1.0) - y0;

                auto t1 = x * cos_p + y * sin_p;
                auto t2 = y * cos_p - x * sin_p;
                if ((t1 * t1 / aa + t2 * t2 / bb) <= (T)1.0)
                    image({i, j}) += intensity;
            }
        }
    }
}

template <typename ImageLike, typename T>
void fill_ellipsoids_(ImageLike& image,
                      const std::vector<ellipsoid<T>>& ellipsoids,
                      volume<3_D, T> image_volume,
                      volume<3_D, T> phantom_volume) {
    // FIXME phi2 and phi3 are ignored
    for (auto e : ellipsoids) {
        auto& intensity = e[0];
        math::vec<3_D, T> a = {e[1], e[2], e[3]};
        math::vec<3_D, T> x_0 = {e[4], e[5], e[6]};
        auto phi = e[7] * M_PI / (T)180.0;
        auto e_1 = math::standard_basis<3_D, T>(1);

        math::vec3<T> origin = image_volume.origin();
        for (int i = 0; i < image_volume.voxels()[0]; ++i) {
            for (int j = 0; j < image_volume.voxels()[1]; ++j) {
                for (int k = 0; k < image_volume.voxels()[2]; ++k) {
                    // we normalize the vector so that it lies in [0,1]
                    math::vec3<T> idx{i, j, k};
                    // component-wise division
                    idx /= image_volume.voxels();

                    // we see where it lies with respect to the phantom, which
                    // depends on:
                    // image size, phantom size, image origin, phantom origin
                    math::vec3<T> v_x = (origin + idx * image_volume.physical_lengths());
                    v_x -= phantom_volume.origin();
                    v_x /= phantom_volume.physical_lengths();

                    auto x = ((T)2.0 * v_x - (T)1.0) - x_0;
                    auto t = math::rotate<T>(x, e_1, phi);
                    if (math::norm<3_D, T>(t / a) <= (T)1.0) {
                        image({i, j, k}) += intensity;
                    }
                }
            }
        }
    }
}

/** A helper function that creates an image from a collection of ellipses. */
template <typename T>
image<2_D, T>
construct_phantom_from_ellipses_(const std::vector<ellipse<T>>& ellipses,
                                 const volume<2_D, T>& v) {
    image<2_D, T> f(v);
    fill_ellipses_(f, ellipses, v, v);
    return f;
}

/** A helper function that creates an image from a collection of ellipses. */
template <typename T>
image<3_D, T>
construct_phantom_from_ellipsoids_(const std::vector<ellipsoid<T>>& ellipsoids,
                                   const volume<3_D, T>& v) {
    image<3_D, T> f(v);
    fill_ellipsoids_(f, ellipsoids, v, v);
    return f;
}

/** Obtain the Shepp-Logan phantom. */
// shepp logan phantoms taken from:
// https://github.com/ismrmrd/ismrmrd-python-tools/blob/master/
// > smrmrdtools/simulation.py
template <typename T>
auto shl_ellipses_() {
    return std::vector<ellipse<T>>{{2, .69, .92, 0, 0, 0},
                                   {-.98, .6624, .8740, 0, -.0184, 0},
                                   {-.02, .1100, .3100, .22, 0, -18},
                                   {-.02, .1600, .4100, -.22, 0, 18},
                                   {.01, .2100, .2500, 0, .35, 0},
                                   {.01, .0460, .0460, 0, .1, 0},
                                   {.02, .0460, .0460, 0, -.1, 0},
                                   {.01, .0460, .0230, -.08, -.605, 0},
                                   {.01, .0230, .0230, 0, -.606, 0},
                                   {.01, .0230, .0460, .06, -.605, 0}};
}

template <typename T>
auto mshl_ellipses_() {
    return std::vector<ellipse<T>>{{1, .69, .92, 0, 0, 0},
                                   {-.80, .6624, .8740, 0, -.0184, 0},
                                   {-.20, .1100, .3100, .22, 0, -18},
                                   {-.20, .1600, .4100, -.22, 0, 18},
                                   {.10, .2100, .2500, 0, .35, 0},
                                   {.10, .0460, .0460, 0, .1, 0},
                                   {.10, .0460, .0460, 0, -.1, 0},
                                   {.10, .0460, .0230, -.08, -.605, 0},
                                   {.10, .0230, .0230, 0, -.606, 0},
                                   {.10, .0230, .0460, .06, -.605, 0}};
}

// these parameters are taken from the ODL library:
// >
// https://github.com/odlgroup/odl/blob/master/odl/phantom/transmission.py
template <typename T>
auto mshl_ellipsoids_() {
    return std::vector<ellipsoid<T>>{
        {1.00, 0.6900, 0.9200, 0.810, 0.0000, 0.0000, 0.00, 0.0, 0.0, 0.0},
        {-0.8, 0.6624, 0.8740, 0.780, 0.0000, -0.0184, 0.00, 0.0, 0.0, 0.0},
        {-0.2, 0.1100, 0.3100, 0.220, 0.2200, 0.0000, 0.00, -18.0, 0.0, 0.0},
        {-0.2, 0.1600, 0.4100, 0.280, -0.2200, 0.0000, 0.00, 18.0, 0.0, 0.0},
        {0.1, 0.2100, 0.2500, 0.410, 0.0000, 0.3500, 0.00, 0.0, 0.0, 0.0},
        {0.1, 0.0460, 0.0460, 0.050, 0.0000, 0.1000, 0.00, 0.0, 0.0, 0.0},
        {0.1, 0.0460, 0.0460, 0.050, 0.0000, -0.1000, 0.00, 0.0, 0.0, 0.0},
        {0.1, 0.0460, 0.0230, 0.050, -0.0800, -0.6050, 0.00, 0.0, 0.0, 0.0},
        {0.1, 0.0230, 0.0230, 0.020, 0.0000, -0.6060, 0.00, 0.0, 0.0, 0.0},
        {0.1, 0.0230, 0.0460, 0.020, 0.0600, -0.6050, 0.00, 0.0, 0.0, 0.0}};
}

template <typename T>
image<2_D, T> shepp_logan_phantom(const volume<2_D, T>& v) {
    return construct_phantom_from_ellipses_<T>(shl_ellipses_<T>(), v);
}

/** Obtain a modified Shepp-Logan phantom with higher contrast. */
template <typename T>
image<2_D, T> modified_shepp_logan_phantom(const volume<2_D, T>& v) {
    return construct_phantom_from_ellipses_<T>(mshl_ellipses_<T>(), v);
}

/** Obtain a modified Shepp-Logan phantom with higher contrast. */
template <typename T>
image<3_D, T> modified_shepp_logan_phantom(const volume<3_D, T>& v) {

    return construct_phantom_from_ellipsoids_<T>(mshl_ellipsoids_<T>(), v);
}

} // namespace tomo
