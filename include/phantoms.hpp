#pragma once

#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "common.hpp"
#include "image.hpp"
#include "volume.hpp"

namespace tomo {

using namespace tomo::literals;

/**
 * An ellipse is defined completely by 6 reals, these are:
 * [I, a, b, x0, y0, phi].
 */
template <typename T>
using ellipse = std::array<T, 6>;

/** Load a phantom from file. */
// template <dimension D, typename T>
// image<D, T> phantom(std::string filename) {}

/** A helper function that creates an image from a collection of ellipses. */
template <typename T>
image<2_D, T>
construct_phantom_from_ellipses_(const std::vector<ellipse<T>>& ellipses,
                                 const volume<2_D>& v) {
    image<2_D, T> f(v);
    for (auto e : ellipses) {
        f[f.index({{1, 2}})] = 1;
        f({{1, 2}}) = 1;

        auto& intensity = e[0];
        auto aa = e[1] * e[1];
        auto bb = e[2] * e[2];
        auto& x0 = e[3];
        auto& y0 = e[4];
        auto phi = e[5] * M_PI / (T)180.0;

        auto cos_p = cos(phi);
        auto sin_p = sin(phi);

        for (int i = 0; i < v.x(); ++i) {
            for (int j = 0; j < v.y(); ++j) {
                auto x = (2.0 * ((double)i / v.x()) - 1) - x0;
                auto y = (2.0 * ((double)j / v.y()) - 1) - y0;
                auto t1 = x * cos_p + y * sin_p;
                auto t2 = y * cos_p - x * sin_p;
                if ((t1 * t1 / aa + t2 * t2 / bb) <= (T)1.0)
                    f({{i, j}}) += intensity;
            }
        }
    }

    return f;
}

/** Obtain the Shepp-Logan phantom. */
// shepp logan phantoms taken from:
// https://github.com/ismrmrd/ismrmrd-python-tools/blob/master/
// > smrmrdtools/simulation.py
template <typename T>
image<2_D, T> shepp_logan_phantom(const volume<2_D>& v) {
    std::vector<ellipse<T>> ellipses = {{{2, .69, .92, 0, 0, 0}},
                                        {{-.98, .6624, .8740, 0, -.0184, 0}},
                                        {{-.02, .1100, .3100, .22, 0, -18}},
                                        {{-.02, .1600, .4100, -.22, 0, 18}},
                                        {{.01, .2100, .2500, 0, .35, 0}},
                                        {{.01, .0460, .0460, 0, .1, 0}},
                                        {{.02, .0460, .0460, 0, -.1, 0}},
                                        {{.01, .0460, .0230, -.08, -.605, 0}},
                                        {{.01, .0230, .0230, 0, -.606, 0}},
                                        {{.01, .0230, .0460, .06, -.605, 0}}};

    return construct_phantom_from_ellipses_<T>(ellipses, v);
}

/** Obtain a modified Shepp-Logan phantom with higher contrast. */
template <typename T>
image<2_D, T> modified_shepp_logan_phantom(const volume<2_D>& v) {
    std::vector<ellipse<T>> ellipses = {{{1, .69, .92, 0, 0, 0}},
                                        {{-.80, .6624, .8740, 0, -.0184, 0}},
                                        {{-.20, .1100, .3100, .22, 0, -18}},
                                        {{-.20, .1600, .4100, -.22, 0, 18}},
                                        {{.10, .2100, .2500, 0, .35, 0}},
                                        {{.10, .0460, .0460, 0, .1, 0}},
                                        {{.10, .0460, .0460, 0, -.1, 0}},
                                        {{.10, .0460, .0230, -.08, -.605, 0}},
                                        {{.10, .0230, .0230, 0, -.606, 0}},
                                        {{.10, .0230, .0460, .06, -.605, 0}}};

    return construct_phantom_from_ellipses_<T>(ellipses, v);
}

} // namespace tomo
