#pragma once

#include <zmq.hpp>

namespace tomo {
namespace util {

struct slice {
    math::mat<4_D, T> orientation;
}

/*
 * the reconstructor should subscribe to slice change updates
 * and maybe update the plotter where necessary?
 *
 * - plotter
 * - geometry
 * - image
 */
class reconstructor {
  public:
    reconstructor();

    virtual tomo::image<2_D> get_slice_data(slice s) = 0;
    virtual tomo::image<3_D> get_volume_data(int resolution) = 0;
};

/**
 * a dummy reconstructor that caches the result of a reconstruction
 * and supplies the data where necessary
 */
class dummy_reconstructor : public reconstructor {
  public:
  private:
};

} // namespace util
} // namespace tomo
