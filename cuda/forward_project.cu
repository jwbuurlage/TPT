#include "fmt/format.h"
#include "hemi/parallel_for.h"
#include <thrust/device_vector.h>

#include "device_data.hpp"
#include "projectors/closest.cu.hpp"

namespace tomo {
namespace cuda {

template <typename T>
__global__ void forward_project_kernel(cudaTextureObject_t image_texture,
                                       device::line<T>* device_lines,
                                       T* device_sino, device::volume v) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    T result = 0;

    project_closest_texture(device_lines[i], v,
                            [&result, image_texture](T x, T y) {
                                result += tex2D<float>(image_texture, x, y);
                            });

    device_sino[i] = result;
}

template <typename T>
void run_forward_projection(cudaTextureObject_t image_texture,
                            device::line<T>* device_lines, int lines,
                            T* host_sino, device::volume v) {
    thrust::device_vector<T> device_sino(lines);

    int threads = 256;
    forward_project_kernel<<<lines / threads, threads>>>(
        image_texture, device_lines, device_sino.data().get(), v);

    cudaMemcpy(host_sino, device_sino.data().get(),
               device_sino.size() * sizeof(T), cudaMemcpyDeviceToHost);
}

template void run_forward_projection<float>(cudaTextureObject_t,
                                            device::line<float>*, int, float*,
                                            device::volume);

} // namespace cuda
} // namespace tomo
