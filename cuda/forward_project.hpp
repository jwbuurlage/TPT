#include "tomo.hpp"

#include <cuda_runtime.h>
#include <thrust/device_vector.h>

#include "bridges.hpp"
#include "device_data.hpp"

namespace tomo {
namespace cuda {

extern void test();

template <typename T>
void run_forward_projection(T* device_image, device::line<T>* device_lines,
                            int lines, T* result, device::volume);

extern template void run_forward_projection<float>(float*, device::line<float>*,
                                                   int, float*, device::volume);
extern template void run_forward_projection<double>(double*,
                                                    device::line<double>*, int,
                                                    double*, device::volume);

struct external_cuda_projector {};

template <dimension D, typename T, class Geometry>
sinogram<D, T, Geometry, external_cuda_projector>
forward_projection(const image<D, T>& f, const Geometry& g) {
    sinogram<D, T, Geometry, external_cuda_projector> result(g);

    // 1. copy data
    device::geometry_bridge<T, decltype(g)> device_geometry(g);
    device::volume_bridge device_volume(f.get_volume());
    device::image_bridge<T> device_image(f);

    // 2. run alg
    run_forward_projection(device_image.get(), device_geometry.get(), g.lines(),
                           result.mutable_data().data(), device_volume.get());

    return result;
}

} // namespace cuda
} // namespace tomo
