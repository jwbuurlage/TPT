#include "tomo.hpp"

#include <cuda_runtime.h>

#include "bridges.hpp"
#include "device_data.hpp"

namespace tomo {
namespace cuda {

template <typename T>
void run_sart(device::volume v, device::line<T>* device_lines, int lines,
              T* device_sino, T* host_image, int group_count, T beta = 0.5,
              int iterations = 10);

extern template void run_sart(device::volume v,
                              device::line<float>* device_lines, int lines,
                              float* device_sino, float* host_image,
                              int group_count, float beta, int iterations);
extern template void run_sart(device::volume v,
                              device::line<double>* device_lines, int lines,
                              double* device_sino, double* host_image,
                              int group_count, double beta, int iterations);

template <dimension D, typename T, class Geometry, class Projector>
image<D, T> sart(const volume<D>& v, const Geometry& g,
                 const sinogram<D, T, Geometry, Projector>& p, T beta = (T)0.5,
                 int iterations = 10) {
    image<D, T> f(v);

    // 1. copy data
    device::geometry_bridge<T, decltype(g)> device_geometry(g);
    device::volume_bridge device_volume(v);
    device::sinogram_bridge<T, decltype(p)> device_sino(p);

    // 2. run alg
    run_sart(device_volume.get(), device_geometry.get(), g.lines(),
             device_sino.get(), f.mutable_data().data(), g.groups()[0], beta,
             iterations);

    return f;
}

} // namespace cuda
} // namespace tomo
