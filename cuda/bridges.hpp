#pragma once

#include "device_data.hpp"
#include "tomo.hpp"

namespace tomo {
namespace cuda {
namespace device {

template <typename T, typename Geometry>
class geometry_bridge {
  public:
    geometry_bridge(Geometry& g) {
        std::vector<device::line<T>> lines(g.lines());
        int k = 0;
        for (auto line : g) {
            lines[k].origin.x = line.origin.x;
            lines[k].origin.y = line.origin.y;
            lines[k].delta.x = line.delta.x;
            lines[k].delta.y = line.delta.y;
            k++;
        }
        auto data_size = sizeof(device::line<T>) * lines.size();
        cudaMalloc(&device_data_, data_size);
        cudaMemcpy(device_data_, lines.data(), data_size,
                   cudaMemcpyHostToDevice);
    }

    ~geometry_bridge() {
        if (device_data_)
            cudaFree(device_data_);
    }

    geometry_bridge(geometry_bridge& other) = delete;
    geometry_bridge(geometry_bridge&& other) = delete;

    device::line<T>* get() { return device_data_; };

  private:
    device::line<T>* device_data_ = nullptr;
};

class volume_bridge {
  public:
    volume_bridge(const tomo::volume<2_D>& host_volume) {
        device_volume_.x = host_volume.x();
        device_volume_.y = host_volume.y();
    }

    device::volume get() { return device_volume_; }

  private:
    device::volume device_volume_;
};

template <typename T>
class image_bridge {
  public:
    image_bridge(const tomo::image<2_D, T>& f) {
        auto data_size = sizeof(T) * f.data().size();
        cudaMalloc(&device_data_, data_size);
        cudaMemcpy(device_data_, f.data().data(), data_size,
                   cudaMemcpyHostToDevice);
    }
    ~image_bridge() {
        if (device_data_)
            cudaFree(device_data_);
    }

    image_bridge(image_bridge& other) = delete;
    image_bridge(image_bridge&& other) = delete;

    T* get() { return device_data_; };

  private:
    T* device_data_ = nullptr;
};

template <typename T, typename Sinogram>
class sinogram_bridge {
  public:
    sinogram_bridge(const Sinogram& sino) {
        auto data_size = sino.geometry().lines() * sizeof(T);
        cudaMalloc(&device_data_, data_size);
        cudaMemcpy(device_data_, sino.data().data(), data_size,
                   cudaMemcpyHostToDevice);
    }
    ~sinogram_bridge() {
        if (device_data_)
            cudaFree(device_data_);
    }

    // disallow copying / moving
    sinogram_bridge(sinogram_bridge& other) = delete;
    sinogram_bridge(sinogram_bridge&& other) = delete;

    T* get() { return device_data_; };

  private:
    T* device_data_ = nullptr;
};

} // namespace device
} // namespace cuda
} // namespace tomo
