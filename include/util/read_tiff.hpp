#pragma once

#include <cstdint>
#include <stdexcept>

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgcodecs.hpp>

#include "../common.hpp"
#include "../image.hpp"

namespace tomo {

class invalid_tiff_file : public std::runtime_error {
    using runtime_error::runtime_error;
};

template <dimension D, typename T = default_scalar_type>
tomo::image<D, T> tiff_to_image(std::string file) {
    cv::Mat image = cv::imread(file, CV_LOAD_IMAGE_ANYDEPTH);

    if (!image.data)
        throw invalid_tiff_file("Could not open or find the image");

    int w = image.cols;
    int h = image.rows;

    auto v = tomo::volume<D, T>({w, h});
    auto result = tomo::image<D, T>(v);
    for (int i = 0; i < w; ++i) {
        for (int j = 0; j < h; ++j) {
            result({i, j}) = (T)image.at<uint8_t>(i, j);
        }
    }

    return result;
}

template <dimension D, typename T = default_scalar_type>
tomo::projections<D, T> tiff_stack_to_projections(const tomo::geometry::base<D, T>& g,
                                                  std::string file_pattern,
                                                  int count) {
    tomo::projections<D, T> result(g);
    // for each file, load as tomo image
    // set the appropriate value in result
    return result;
}

} // namespace tomo
