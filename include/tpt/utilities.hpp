#pragma once

#include <algorithm>
#include <cassert>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "common.hpp"
#include "geometries/list.hpp"
#include "image.hpp"
#include "projections.hpp"

namespace tpt {

template <typename T>
image<2_D, T> slice(const image<3_D, T>& f, int slice, int along_axis = 0) {
    std::array<int, 2> axes;
    int idx = 0;
    for (int d = 0; d < 3; ++d) {
        if (d == along_axis)
            continue;
        axes[idx++] = d;
    }

    auto large_volume = f.get_volume();

    volume<2_D, T> v(math::restrict<3_D, T>(large_volume.voxels(), along_axis));
    image<2_D, T> slice_image(v);

    math::vec<3_D, int> image_idx;
    image_idx[along_axis] = slice;
    for (int i = 0; i < large_volume.voxels()[axes[0]]; ++i) {
        image_idx[axes[0]] = i;
        for (int j = 0; j < large_volume.voxels()[axes[1]]; ++j) {
            image_idx[axes[1]] = j;
            slice_image({i, j}) = f(image_idx);
        }
    }

    return slice_image;
}

/** A helper function that downscales the image for plotting in the
 * terminal. */
template <dimension D, typename T>
image<D, T> downscale(const image<D, T>& f, math::vec<D, int> new_size) {
    volume<D, T> new_volume(new_size);
    image<D, T> g(new_volume);

    auto voxels = f.get_volume().voxels();
    auto stride = math::vec3<T>(voxels) / math::vec3<T>(new_size);

    for (auto i = 0u; i < new_volume.cells(); ++i) {
        g[i] = f(math::vec3<T>(new_volume.unroll(i)) * stride);
    }

    return g;
}

/** A helper function that gets the data for a given 2D slice of a 3D volume */
template <typename T>
image<2_D, T> slice_of_image(const image<3_D, T>& f, math::slice<T> s) {
    volume<2_D, T> new_volume(f.get_volume().voxels()[0]);
    image<2_D, T> g(new_volume);

    // s.base_point;
    for (int i = 0; i < new_volume.voxels()[0]; ++i) {
        for (int j = 0; j < new_volume.voxels()[1]; ++j) {
            // get voxel index of the following vector (in [-1, 1]^3 space)
            auto point = s.base + ((T)i / new_volume.voxels()[0]) * s.b_x +
                         ((T)j / new_volume.voxels()[1]) * s.b_y;
            auto normalized_point = (T)0.5 * (point + math::vec3<T>(1));
            auto voxel_coordinate = math::vec3<int>(
                normalized_point * math::vec3<T>(f.get_volume().voxels()));
            if (math::valid_index(voxel_coordinate, f.get_volume())) {
                g({i, j}) = f(voxel_coordinate);
            }
        }
    }

    return g;
}

/** A helper function that downscales the image for plotting in the
 * terminal. */
template <typename T>
image<2_D, T> downscale_(const image<2_D, T>& f, math::vec<2_D, int> new_size) {
    volume<2_D, T> new_volume(new_size);
    image<2_D, T> g(new_volume);
    // TODO: This should be the area of the one in the other..
    auto scale = ((T)new_size[0]) / f.get_volume().voxels()[0];

    for (int i = 0; i < f.size(0); ++i) {
        for (int j = 0; j < f.size(1); ++j) {
            // contributes to cell:
            auto cell_x = (int)(((T)i / f.size(0)) * new_volume.voxels()[0]);
            auto cell_y = (int)(((T)j / f.size(1)) * new_volume.voxels()[1]);
            g({cell_x, cell_y}) += f({i, j}) * scale;
        }
    }

    return g;
}

/** Print the image as ascii art to stdout. */
template <typename T>
void ascii_plot(const image<2_D, T>& f, T max = -1) {
    const int limit = 40;
    int max_size = f.size(0) > f.size(1) ? f.size(0) : f.size(1);

    if (max_size > limit) {
        T scale = (T)limit / (T)max_size;
        auto downscaled_img =
            downscale_(f, {f.size(0) * scale, f.size(1) * scale});
        ascii_plot_output(downscaled_img,
                          {f.size(0) * scale, f.size(1) * scale}, max);
    } else {
        ascii_plot_output(f, {f.size(0), f.size(1)}, max);
    }
}

/** Print the 3d image in slices as ascii art to stdout. */
template <typename T>
void ascii_plot(const image<3_D, T>& f, int slices = 4, int axis = 0) {
    T max = (T)0;
    for (auto k = 0u; k < f.get_volume().cells(); ++k)
        if (f[k] > max)
            max = f[k];

    for (int i = 0; i < slices; ++i) {
        auto delta = (T)f.get_volume().voxels()[axis] / (slices + 1);
        auto slice_idx = (int)(delta * (i + 1));
        auto g = slice(f, slice_idx, axis);
        std::cout << "Slice @ x_" << axis << " = " << slice_idx << "\n";
        ascii_plot(g, max);
    }
}

/** Print the projections as ascii art to stdout. */
template <dimension D, typename T>
void ascii_plot(const projections<D, T>& sino) {
    ascii_plot_output(sino, sino.geometry().groups());
}

/** Output an image-like object to the standard output. */
template <class ImageLike>
void ascii_plot_output(ImageLike& image, math::vec2<int> dimensions,
                       typename ImageLike::value_type max = -1) {
    using T = typename ImageLike::value_type;
    using namespace std::string_literals;

    auto chars = "$@B%8&WM#*oahkbdpqwmZO0QLCJUYXzcvunxrjft/"
                 "\\|()1{}[]?-_+~<>i!lI;:,\"^`'. "s;
    std::reverse(chars.begin(), chars.end());

    T min = std::numeric_limits<T>::max();
    for (auto k = 0u; k < image.get_volume().cells(); ++k) {
        if (image[k] < min) {
            min = image[k];
        }
    }

    if (max < 0) {
        for (auto k = 0u; k < image.get_volume().cells(); ++k) {
            if (image[k] > max) {
                max = image[k];
            }
        }
    }

    std::cout << min << " (min) \n";
    std::cout << max << " (max) \n";

    max -= min;

    int cur = 0;
    for (int j = 0; j < dimensions[1]; ++j) {
        for (int i = 0; i < dimensions[0]; ++i) {
            auto idx = (std::size_t)(
                math::sqrt(std::max(image[cur++] - min, (T)0.0) / max) *
                chars.size());
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
auto random_list_geometry(int k, volume<D, T> v) {
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
        // Note that a D-dimensional volume has 2 * D (hyper)sides. We
        // choose
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
    // Here we use that the Guassian distribution is spherically symmetric,
    // so
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
    std::vector<math::ray<D, T>> line_list(k);
    for (int i = 0; i < k; ++i) {
        line_list[i].source = random_origin[i] - random_direction[i];
        line_list[i].detector = line_list[i].source +
                                random_direction[i] * math::sqrt<T>(D) *
                                    (T)2.0 * random_direction[i];
    }

    auto g = geometry::list<D, T>(std::move(line_list));
    g.set_dimensions({k, 1});
    return g;
}

/** Compute the 'zip' of two vectors. */
template <typename T, typename U>
std::vector<std::pair<T, U>> zip(std::vector<T> xs, std::vector<U> ys) {
    std::vector<std::pair<T, U>> result;
    for (size_t i = 0; i < (xs.size() < ys.size() ? xs.size() : ys.size());
         ++i) {
        result.push_back(make_pair(xs[i], ys[i]));
    }
    return result;
}

} // namespace tpt
