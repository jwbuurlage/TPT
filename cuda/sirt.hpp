
#include "tomo.hpp"

#include <cuda_runtime.h>

#include "bridges.hpp"
#include "device_data.hpp"

namespace tomo {
namespace cuda {

template <typename T>
void run_sirt(cudaTextureObject_t, cudaSurfaceObject_t, cudaTextureObject_t,
              cudaSurfaceObject_t, device::volume v,
              device::line<T>* device_lines, int lines, T* device_sino,
              int group_count, T beta = 0.5,
              int iterations = 10);

extern template void run_sirt(cudaTextureObject_t, cudaSurfaceObject_t,
                              cudaTextureObject_t, cudaSurfaceObject_t,
                              device::volume, device::line<float>*, int,
                              float* device_sino, int, float,
                              int);

template <dimension D, typename T, class Geometry, class Projector>
image<D, T> sirt(const volume<D>& v, const Geometry& g,
                 const sinogram<D, T, Geometry, Projector>& p, T beta = (T)0.5,
                 int iterations = 10) {
    image<D, T> f(v);

    // 1. copy data
    device::geometry_bridge<T, decltype(g)> device_geometry(g, v);
    device::volume_bridge device_volume(v);
    device::sinogram_bridge<T, decltype(p)> device_sino(p);

    // 1b. allocate buffers
    device::image_bridge<T> image(v);
    device::image_bridge<T> image_buffer(v);

    // 2. run alg
    run_sirt(image.get_texture(), image.get_surface(),
             image_buffer.get_texture(), image_buffer.get_surface(),
             device_volume.get(), device_geometry.get(), g.lines(),
             device_sino.get(), g.groups()[0], beta,
             iterations);

    // 3. copy image from device
    image.copy_to(f);

    return f;
}

} // namespace cuda
} // namespace tomo
