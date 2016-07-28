#include "fmt/format.h"
#include "hemi/parallel_for.h"

#include "device_data.hpp"
#include "projectors/closest.cu.hpp"

namespace tomo {
namespace cuda {

template <typename T>
__global__ void forward_project_kernel(T* device_image,
                                       device::line<T>* device_lines,
                                       T* device_sino, device::volume v) {
    // gpu, first just 'closest' kernel (value 1)
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    T result = 0;
    project_closest(device_lines[i], v, [&result, &device_image](int index) {
        result += device_image[index];
    });
    device_sino[i] = result;
}

template <typename T>
void run_forward_projection(T* device_image, device::line<T>* device_lines,
                            int lines, T* host_sino, device::volume v) {
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

template void run_forward_projection<float>(float*, device::line<float>*, int,
                                            float*, device::volume);

template void run_forward_projection<double>(double*, device::line<double>*, int,
                                             double*, device::volume);

} // namespace cuda
} // namespace tomo
