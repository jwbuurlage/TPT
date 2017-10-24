#pragma once

#include <cstdint>
#include <regex>
#include <stdexcept>
#include <string>
#include <limits>

#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

#include <opencv2/opencv.hpp>

#include "../common.hpp"
#include "../image.hpp"
#include "../utilities.hpp"

namespace tomo {

class invalid_tiff_file : public std::runtime_error {
    using runtime_error::runtime_error;
};

template <typename T = default_scalar_type>
tomo::image<2_D, T> tiff_to_image(std::string filename) {
//    cv::Mat image = cv::imread(filename, CV_LOAD_IMAGE_ANYDEPTH);
//
//    if (!image.data)
//        throw invalid_tiff_file("Could not open or find the image: " +
//                                filename);
//
//    int w = image.cols;
//    int h = image.rows;
//
//    auto v = tomo::volume<D, T>({w, h});
//    auto result = tomo::image<D, T>(v);
//    for (int i = 0; i < w; ++i) {
//        for (int j = 0; j < h; ++j) {
//            result({i, j}) = (T)image.at<int32_t>(i, j);
//        }
//    }
(void)filename;

    auto v = tomo::volume<2_D, T>({1, 1});
    return tomo::image<2_D, T>(v);
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
        auto image = tiff_to_image<T>(filename);
        result.set_projection(i, image);
    }

    return result;
}

template <typename T = default_scalar_type>
void write_png(tomo::image<2_D, T> x, std::string filename) {
    auto v = x.get_volume();

    T min = 0.0;
    for (auto y : x.data()) {
        if (y < min) {
            min = y;
        }
    }
    for (auto& y : x.mutable_data()) {
        y -= min;
    }
    T max = 0.0;
    for (auto y : x.data()) {
        if (y > max) {
            max = y;
        }
    }

    auto image = cv::Mat(v.voxels()[0], v.voxels()[1], CV_16UC1);
    for (int i = 0; i < v.voxels()[0]; ++i) {
        for (int j = 0; j < v.voxels()[1]; ++j) {
            image.at<unsigned short>(i, j) =
                (unsigned short)(x[x.index({i, j})] / max *
                                std::numeric_limits<unsigned short>::max());
        }
    }

    std::vector<int> compression_params;
    compression_params.push_back(CV_IMWRITE_PNG_COMPRESSION);
    compression_params.push_back(9);

    try {
	IplImage tmp = image;
        cvSaveImage((filename + ".png").c_str(), &tmp);
    } catch (std::runtime_error& ex) {
        fprintf(stderr, "Exception converting image to PNG format: %s\n",
                ex.what());
    }

    std::cout << "Saved: " << filename << ".png\n";
}


template <typename T = default_scalar_type>
void write_png(tomo::image<3_D, T> x, std::string filename) {
    write_png(tomo::slice(x, x.get_volume()[0] / 2, 0), filename + "_x");
    write_png(tomo::slice(x, x.get_volume()[1] / 2, 1), filename + "_y");
    write_png(tomo::slice(x, x.get_volume()[2] / 2, 2), filename + "_z"); 
}

} // namespace tomo
