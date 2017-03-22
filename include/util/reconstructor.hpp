#pragma once

#include <mutex>
#include <vector>

#include <zmq.hpp>

#include "../algorithms/sirt.hpp"
#include "../geometries/parallel.hpp"
#include "../math.hpp"
#include "../operations.hpp"
#include "../phantoms.hpp"
#include "../projectors/joseph.hpp"

namespace tomo {
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

    virtual tomo::image<2_D, T> get_slice_data(math::slice<T> s) = 0;
    virtual tomo::image<3_D, T> get_volume_data(int resolution) = 0;

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
    dummy_reconstructor(tomo::volume<3_D, T> volume)
        : volume_(volume), current_image_(volume) {}

    void reconstruct(int size) {
        auto v = volume_;
        auto g = tomo::geometry::parallel<3_D, T>(v, size, size);
        auto f = tomo::modified_shepp_logan_phantom<T>(v);

        auto kernel = tomo::dim::joseph<3_D, T>(v);

        this->notify();

        auto sino = tomo::forward_projection<3_D, T>(f, g, kernel);
        auto image = tomo::reconstruction::sirt(
            v, g, kernel, sino, 0.5, 10,
            {[&](tomo::image<3_D, T>& iteration_result) {
                {
                    std::lock_guard<std::mutex> guard(image_mutex_);
                    current_image_ = iteration_result;
                }
                this->notify();
            }});

        std::lock_guard<std::mutex> guard(image_mutex_);
        current_image_ = image;
    }

    tomo::image<2_D, T> get_slice_data(math::slice<T> s) override {
        std::lock_guard<std::mutex> guard(image_mutex_);
        return slice_of_image(current_image_, s);
    }

    tomo::image<3_D, T> get_volume_data(int resolution) override {
        std::lock_guard<std::mutex> guard(image_mutex_);
        return downscale<3_D, T>(current_image_, math::vec3<int>{resolution});;
    }

  private:
    tomo::volume<3_D, T> volume_;

    image<3_D, T> current_image_;
    /* Note: we will have to return the images by value. */
    std::mutex image_mutex_; // protects current_image_
};

} // namespace util
} // namespace tomo
