#include <cmath>
#include <thrust/device_vector.h>

#include "fmt/format.h"
#include "hemi/parallel_for.h"

#include "device_data.hpp"
#include "projectors/closest.cu.hpp"

namespace tomo {
namespace cuda {

#define EPSILON 1e-6

template <typename T>
__global__ void w_norms_kernel(const device::line<T>* device_lines,
                               device::volume v, T* w_norms) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    T result = 0;
    project_closest(device_lines[i], v, [&result](int index) { result += 1; });
    w_norms[i] = result;
}

template <typename T>
__global__ void
sart_kernel(T* device_image, const device::line<T>* device_lines,
            const T* device_sino, device::volume v, T beta, const T* w_norms) {
    // gpu, first just 'closest' kernel (value 1)
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    T alpha = 0;
    project_closest(device_lines[i], v, [&alpha, &device_image](int index) {
        alpha += device_image[index];
    });

    auto factor = beta * (device_sino[i] - alpha) / w_norms[i];

    project_closest(device_lines[i], v, [&factor, &device_image](int index) {
        device_image[index] += factor;
    });
}

template <typename T>
void run_sart(device::volume v, device::line<T>* device_lines, int lines,
              T* device_sino, T* host_image, int group_count, T beta = 0.5,
              int iterations = 10) {
    thrust::device_vector<T> device_image(v.x * v.y);

    int group_size = lines / group_count;
    int threads = 256;

    thrust::device_vector<T> w_norms(lines);
    w_norms_kernel<<<lines / threads, threads>>>(device_lines, v,
                                                 w_norms.data().get());

    for (int i = 0; i < iterations; ++i) {
        for (int k = 0; k < group_count; ++k) {
            sart_kernel<<<group_size / threads, threads>>>(
                device_image.data().get(), &device_lines[k * group_size],
                &device_sino[k * group_size], v, beta,
                &w_norms.data().get()[k * group_size]);
        }
    }

    cudaMemcpy(host_image, device_image.data().get(),
               device_image.size() * sizeof(T), cudaMemcpyDeviceToHost);
}

template void run_sart(device::volume v, device::line<float>* device_lines,
                       int lines, float* device_sino, float* host_image, int,
                       float beta, int iterations);

template void run_sart(device::volume v, device::line<double>* device_lines,
                       int lines, double* device_sino, double* host_image,
                       int group_count, double beta, int iterations);

} // namespace cuda
} // namespace tomo
