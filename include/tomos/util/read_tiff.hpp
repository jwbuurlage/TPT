#pragma once

#include <cstdint>
#include <regex>
#include <stdexcept>
#include <string>

#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

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
tomo::image<D, T> tiff_to_image(std::string filename) {
    cv::Mat image = cv::imread(filename, CV_LOAD_IMAGE_ANYDEPTH);

    if (!image.data)
        throw invalid_tiff_file("Could not open or find the image: " +
                                filename);

    int w = image.cols;
    int h = image.rows;

    auto v = tomo::volume<D, T>({w, h});
    auto result = tomo::image<D, T>(v);
    for (int i = 0; i < w; ++i) {
        for (int j = 0; j < h; ++j) {
            result({i, j}) = (T)image.at<int32_t>(i, j);
        }
    }

    return result;
}

template <dimension D, typename T = default_scalar_type>
tomo::projections<D, T>
tiff_stack_to_projections(const tomo::geometry::base<D, T>& g,
                          std::string filename_pattern, int projection_count,
                          fs::path root_directory = {}) {
    // TODO implement for other D's
    static_assert(D == 3_D, "can only read 3D projection data");

    tomo::projections<D, T> result(g);

    for (int i = 0; i < projection_count; ++i) {
        auto filename =
            std::regex_replace(root_directory.string() + filename_pattern,
                               std::regex("[*]"), std::to_string(i));
        auto image = tiff_to_image<D - 1, T>(filename);
        result.set_projection(i, image);
    }

    return result;
}

} // namespace tomo
