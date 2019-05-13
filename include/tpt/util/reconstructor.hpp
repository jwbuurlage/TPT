#pragma once

#include <mutex>
#include <vector>

#include <zmq.hpp>

#include "../algorithms/sirt.hpp"
#include "../geometry.hpp"
#include "../image.hpp"
#include "../math.hpp"
#include "../operations.hpp"
#include "../phantoms.hpp"
#include "../projections.hpp"
#include "../projector.hpp"

namespace tpt {
namespace util {

template <typename T>
class on_demand_reconstructor;

template <typename T>
class data_update_observer {
  public:
    virtual void
    new_data_available(on_demand_reconstructor<T>& the_reconstructor) = 0;

  protected:
    ~data_update_observer() = default;
};

/*
 * the reconstructor should subscribe to slice change updates
 * and maybe update the plotter where necessary?
 *
 * - plotter
 * - geometry
 * - image
 */
template <typename T>
class on_demand_reconstructor {
  public:
    on_demand_reconstructor() = default;
    virtual ~on_demand_reconstructor() = default;

    virtual image<2_D, T> get_slice_data(math::slice<T> s) = 0;
    virtual image<3_D, T> get_volume_data(int resolution) = 0;

    void add_observer(data_update_observer<T>* observer) {
        observers_.push_back(observer);
    }

    void notify() {
        for (auto observer : observers_) {
            observer->new_data_available(*this);
        }
    }

  protected:
    std::vector<data_update_observer<T>*> observers_;
};

/**
 * a dummy reconstructor that caches the result of a reconstruction
 * and supplies the data where necessary
 */
template <typename T>
class dummy_reconstructor : public on_demand_reconstructor<T> {
  public:
    dummy_reconstructor(volume<3_D, T> volume, dim::base<3_D, T>& kernel,
                        geometry::base<3_D, T>& geometry,
                        image<3_D, T>& phantom_image,
                        projections<3_D, T> projection_stack)
        : volume_(volume), current_image_(volume), kernel_(kernel),
          geometry_(geometry), phantom_image_(phantom_image),
          projection_stack_(projection_stack) {}

    void reconstruct() {
        current_image_ = phantom_image_;
        this->notify();

        auto image = reconstruction::sirt(
            volume_, geometry_, kernel_, projection_stack_, 0.5, 10,
            {[&](auto& iteration_result, int) {
                {
                    std::lock_guard<std::mutex> guard(image_mutex_);
                    current_image_ = iteration_result;
                }
                this->notify();
            }});

        std::lock_guard<std::mutex> guard(image_mutex_);
        current_image_ = image;
    }

    image<2_D, T> get_slice_data(math::slice<T> s) override {
        std::lock_guard<std::mutex> guard(image_mutex_);
        return slice_of_image(current_image_, s);
    }

    image<3_D, T> get_volume_data(int resolution) override {
        std::lock_guard<std::mutex> guard(image_mutex_);
        return downscale<3_D, T>(current_image_, math::vec3<int>{resolution});
        ;
    }

  private:
    volume<3_D, T> volume_;

    image<3_D, T> current_image_;
    /* Note: we will have to return the images by value. */
    std::mutex image_mutex_; // protects current_image_

    dim::base<3_D, T>& kernel_;
    geometry::base<3_D, T>& geometry_;
    image<3_D, T>& phantom_image_;
    projections<3_D, T> projection_stack_;
};

} // namespace util
} // namespace tpt
