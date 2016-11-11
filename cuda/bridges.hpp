#pragma once

#include <cassert>

#include "device_data.hpp"
#include "tomo.hpp"

namespace tomo {
namespace cuda {
namespace device {

// TODO:
// lines are now rays
// should fix all this
// ...
template <typename T, typename Geometry>
class geometry_bridge {
  public:
    geometry_bridge(Geometry& g, tomo::volume<2_D> v) {
        std::vector<device::line<T>> lines(g.lines());
        int k = 0;
        for (auto line : g) {
            lines[k].origin.x = line.origin.x;
            lines[k].origin.y = line.origin.y;
            lines[k].delta.x = line.delta.x;
            lines[k].delta.y = line.delta.y;

            auto far_away = line.origin + math::sqrt2<T> * v.x() * line.delta;

            // compute number of steps
            auto end_point =
                math::box_intersection<T>(line.origin + line.delta, far_away,
                                          math::vec2<T>(v.x(), v.y()));
            lines[k].steps =
                (int)(math::distance<2_D, T>(line.origin, end_point)) - 1;
            k++;
        }
        auto line_size = sizeof(device::line<T>);

        auto groups = g.groups();
        // FIXME pitch not necessarily aligned
        cudaMallocPitch(&device_data_, &pitch_, groups[0] * line_size,
                        groups[1]);
        assert(pitch_ == groups[0] * line_size);

        cudaMemcpy2D(device_data_, pitch_, lines.data(), groups[0] * line_size,
                     groups[0] * line_size, groups[1], cudaMemcpyHostToDevice);
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
    size_t pitch_;
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
    image_bridge(const tomo::image<2_D, T>& f)
        : texture_object_{}, surface_object_{} {
        auto data_size = sizeof(T) * f.data().size();
        cudaMalloc(&device_data_, data_size);
        cudaMemcpy(device_data_, f.data().data(), data_size,
                   cudaMemcpyHostToDevice);

        cudaChannelFormatDesc desc =
            cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
        cudaMallocArray(&device_array_, &desc, f.get_volume().x(),
                        f.get_volume().y());
        cudaMemcpyToArray(device_array_, 0, 0, f.data().data(), data_size,
                          cudaMemcpyHostToDevice);

        create_objects_();
    }

    //    image_bridge(tomo::volume<2_D> v) {
    //        auto data_size = sizeof(T) * v.cells();
    //        cudaMalloc(&device_data_, data_size);
    //        cudaMemset(device_data_, 0, data_size);
    //
    //        cudaChannelFormatDesc desc =
    //            cudaCreateChannelDesc(32, 0, 0, 0,
    //            cudaChannelFormatKindFloat);
    //        cudaMallocArray(&device_array_, &desc, v.x(), v.y(),
    //        cudaArrayDefault);
    //        cudaMemcpyToArray(device_array_, 0, 0, device_data_, data_size,
    //                          cudaMemcpyDeviceToDevice);
    //
    //        create_objects_();
    //    };

    ~image_bridge() {
        if (device_data_)
            cudaFree(device_data_);
        if (device_array_) {
            cudaDestroySurfaceObject(surface_object_);
            cudaDestroyTextureObject(texture_object_);

            cudaFreeArray(device_array_);
        }
    }

    image_bridge(image_bridge& other) = delete;
    image_bridge(image_bridge&& other) = delete;

    T* get() { return device_data_; };
    cudaArray* get_array() { return device_array_; }

    auto get_surface() { return surface_object_; }
    auto get_texture() { return texture_object_; }

    void copy_to(tomo::image<2_D, T>& f) {
        auto data_size = sizeof(T) * f.data().size();
        cudaMemcpyFromArray(f.mutable_data().data(), device_array_, 0, 0,
                            data_size, cudaMemcpyDeviceToHost);
    }

  private:
    void create_objects_() {

        // specify resource type
        struct cudaResourceDesc resDesc = {};
        resDesc.resType = cudaResourceTypeArray;
        resDesc.res.array.array = device_array_;

        // specify texture object parameters
        struct cudaTextureDesc texDesc = {};
        texDesc.addressMode[0] = cudaAddressModeWrap;
        texDesc.addressMode[1] = cudaAddressModeWrap;
        texDesc.filterMode = cudaFilterModeLinear;
        texDesc.readMode = cudaReadModeElementType;
        texDesc.normalizedCoords = 0;

        // create texture object
        cudaCreateTextureObject(&texture_object_, &resDesc, &texDesc, NULL);

        // create surface object
        cudaCreateSurfaceObject(&surface_object_, &resDesc);
    }

    T* device_data_ = nullptr;
    cudaArray* device_array_ = nullptr;
    cudaTextureObject_t texture_object_;
    cudaSurfaceObject_t surface_object_;
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
