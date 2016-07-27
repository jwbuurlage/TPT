#include "hemi/parallel_for.h"

#include "device_data.hpp"

#include "fmt/format.h"

namespace tomo {
namespace cuda {

#define EPSILON 1e-6

template <typename T>
__global__ void forward_project_kernel(T* device_image,
                                       device_line<T>* device_lines,
                                       T* device_sino, device_volume v) {
    // gpu, first just 'closest' kernel (value 1)
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    T result = 0;
    vec2<T> current = device_lines[i].origin;
    while (current.x <= v.x + EPSILON && current.y <= v.y + EPSILON && current.x + EPSILON >= 0 &&
           current.y + EPSILON >= 0) {
        int index = (int)(current.x + (T)0.5) + v.x * (int)(current.y + (T)0.5);
        result += device_image[index];
        current.x += device_lines[i].delta.x;
        current.y += device_lines[i].delta.y;
    }
    device_sino[i] = result;
}

template <typename T>
void run_forward_projection(T* device_image, device_line<T>* device_lines,
                            int lines, T* host_sino, device_volume v) {
    T* device_sino = nullptr;
    auto sino_bytes = lines * sizeof(T);
    cudaMalloc(&device_sino, sino_bytes);
    cudaMemset(device_sino, 0, sino_bytes);

    int threads = 256;
    forward_project_kernel<<<lines / threads, threads>>>(
        device_image, device_lines, device_sino, v);

    cudaMemcpy(host_sino, device_sino, sino_bytes, cudaMemcpyDeviceToHost);

    cudaFree(device_sino);
}

template void run_forward_projection<float>(float*, device_line<float>*, int,
                                            float*, device_volume);

template void run_forward_projection<double>(double*, device_line<double>*, int,
                                             double*, device_volume);

} // namespace cuda
} // namespace tomo
