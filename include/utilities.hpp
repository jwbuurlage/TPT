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

template <dimension D, typename T>
void plot(const image<D, T>& f, std::string filename);

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

template <dimension D, typename T, class G, class P>
void ascii_plot(const sinogram<D, T, G, P>& sino) {
    ascii_plot_output(sino, sino.geometry().groups());
}

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

template <dimension D, typename T>
auto random_list_geometry(tomo::volume<D> v, int k) {
    static_assert(D == 2_D, "only the generation of 2D lines is supported");

    // initialize random set of lines, this implies \phi for now
    std::vector<tomo::line<D, T>> line_list;

    // Seed with a real random value, if available
    std::random_device r;

    // Choose a random mean between 1 and 6
    std::default_random_engine e(r());
    int precision = 123456;
    std::uniform_int_distribution<int> rand(0, precision);

    for (int i = 0; i < k; ++i) {
        if (rand(e) % 2 == 0) {
            // place left
            double offset = v.y() * ((double)rand(e) / precision);
            double angle = M_PI * ((double)rand(e) / precision) - M_PI / 2;
            line_list.push_back(
                {{0, offset},
                 {tomo::math::cos(angle), tomo::math::sin(angle)}});
        } else {
            // place bottom
            double offset = v.x() * ((double)rand(e) / precision);
            double angle = M_PI * ((double)rand(e) / precision);
            line_list.push_back(
                {{offset, 0},
                 {tomo::math::cos(angle), tomo::math::sin(angle)}});
        }
    }

    auto g = tomo::list_geometry<D>(std::move(line_list));
    g.set_dimensions({k, 1});
    return g;
}

} // namespace tomo
