#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "common.hpp"
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
            auto idx = (unsigned int)(math::sqrt(image[cur++] / max) * chars.size());
            if (idx >= chars.size())
                idx = chars.size() - 1;
            if (idx < 0)
                idx = 0;
            std::cout << chars[idx] << " ";
        }
        std::cout << "\n";
    }
}

} // namespace tomo
