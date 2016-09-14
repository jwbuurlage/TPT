#pragma once

#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "common.hpp"
#include "geometries/list.hpp"
#include "image.hpp"
#include "sinogram.hpp"

namespace tomo {

/** A helper function that downscales the image for plotting in the terminal. */
template <typename T>
image<2_D, T> downscale_(const image<2_D, T>& f, volume<2_D> new_volume) {
    assert(f.size(0) == f.size(1));
    image<2_D, T> g(new_volume);

    double scale = ((double)new_volume[0] * new_volume[0]) /
                   ((double)f.size(0) * f.size(0));

    for (int i = 0; i < f.size(0); ++i) {
        for (int j = 0; j < f.size(1); ++j) {
            // contributes to cell:
            auto cell_x = (int)(((double)i / f.size(0)) * new_volume[0]);
            auto cell_y = (int)(((double)j / f.size(1)) * new_volume[1]);
            g({cell_x, cell_y}) += f({i, j}) * scale;
        }
    }

    return g;
}

/** Print the image as ascii art to stdout. */
template <typename T>
void ascii_plot(const image<2_D, T>& f) {
    int limit = 40;
    if (f.size(0) > limit || f.size(1) > limit) {
        ascii_plot_output(downscale_(f, volume<2_D>(limit, limit)),
                          {limit, limit});
    } else {
        ascii_plot_output(f, {f.size(0), f.size(1)});
    }
}

/** Print the sinogram as ascii art to stdout. */
template <dimension D, typename T, class G, class P>
void ascii_plot(const sinogram<D, T, G, P>& sino) {
    ascii_plot_output(sino, sino.geometry().groups());
}

/** Output an image-like object to the standard output. */
template <class ImageLike>
void ascii_plot_output(const ImageLike& image, math::vec2<int> dimensions) {
    using namespace std::string_literals;
    auto chars = "$@B%8&WM#*oahkbdpqwmZO0QLCJUYXzcvunxrjft/"
                 "\\|()1{}[]?-_+~<>i!lI;:,\"^`'. "s;
    std::reverse(chars.begin(), chars.end());

    double max = 0;
    for (int k = 0; k < dimensions[0] * dimensions[1]; ++k)
        if (image[k] > max)
            max = image[k];

    int cur = 0;
    for (int i = 0; i < dimensions[0]; ++i) {
        for (int j = 0; j < dimensions[1]; ++j) {
            auto idx =
                (unsigned int)(math::sqrt(image[cur++] / max) * chars.size());
            if (idx >= chars.size())
                idx = chars.size() - 1;
            std::cout << chars[idx] << " ";
        }
        std::cout << "\n";
    }
}

/**
 * Create a random list geometry in D dimensions.
 *
 * \tparam D the dimension of the scanning volume
 * \tparam T the scalar type to use
 *
 * \param k the number of lines to generate
 * \param v the volume to be scanned
 */
template <dimension D, typename T>
auto random_list_geometry(int k, volume<D> v) {
    // Seed with a real random value, if available
    std::random_device r;
    std::default_random_engine e(r());
    T precision = (T)100.0;

    std::uniform_int_distribution<int> rand_int(0, 1234567);
    std::normal_distribution<T> gaussian(0, precision);
    std::uniform_real_distribution<T> rand_pos(0.5, v[0] - 0.5);

    std::vector<math::vec<D, T>> random_origin(k);
    std::vector<int> chosen_ds;
    std::vector<int> chosen_planes;
    for (auto& origin : random_origin) {
        // The origin should be on the boundary of the volume
        // Note that a D-dimensional volume has 2 * D (hyper)sides. We choose
        // such a side
        // for our origin to lie on.

        // First we draw a random dimension
        int chosen_d = rand_int(e) % D;
        chosen_ds.push_back(chosen_d);

        // Next we choose one of the two planes
        int plane = rand_int(e) % 2;
        chosen_planes.push_back(plane);

        // Now we generate a random origin, here we make sure it is not 'too
        // much' to the sides
        for (int d = 0; d < D; ++d) {
            if (d == chosen_d) {
                origin[d] = (plane == 0) ? (T)0.0 : (T)v[0];
            } else {
                origin[d] = rand_pos(e);
            }
        }
    }

    // Next we generate random directions
    // Here we use that the Guassian distribution is spherically symmetric, so
    // that we can generate a point on the sphere using D random gaussian
    // variables.
    std::vector<math::vec<D, T>> random_direction(k);
    int j = 0;
    for (auto& direction : random_direction) {
        do {
            for (int i = 0; i < D; ++i) {
                direction[i] = gaussian(e);
            }
        } while (direction == math::vec<D, T>((T)0.0));
        direction = math::normalize(direction);
        // now we fix the direction such that it always points inward
        if ((chosen_planes[j] == 0 && (direction[chosen_ds[j]] < 0)) ||
            (chosen_planes[j] == 1 && (direction[chosen_ds[j]] > 0))) {
            direction[chosen_ds[j]] = -direction[chosen_ds[j]];
        }
        ++j;
    }

    // initialize random set of lines, this implies \phi for now
    std::vector<line<D, T>> line_list(k);
    for (int i = 0; i < k; ++i) {
        line_list[i].origin = random_origin[i] + (T)0.5 * random_direction[i];
        line_list[i].delta = random_direction[i];
    }

    auto g = geometry::list<D, T>(std::move(line_list));
    g.set_dimensions({k, 1});
    return g;
}

} // namespace tomo
