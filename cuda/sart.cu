#include "hemi/parallel_for.h"
#include <cmath>

#include "device_data.hpp"

#include "fmt/format.h"

namespace tomo {
namespace cuda {

#define EPSILON 1e-6

template <typename T>
__global__ void w_norms_kernel(const device_line<T>* device_lines,
                               device_volume v, T* w_norms) {
    // gpu, first just 'closest' kernel (value 1)
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    T result = 0;
    vec2<T> current = device_lines[i].origin;
    while (current.x <= v.x + EPSILON && current.y <= v.y + EPSILON &&
           current.x + EPSILON >= 0 && current.y + EPSILON >= 0) {
        int index = (int)(current.x) + v.x * (int)(current.y);
        if (index < v.x * v.y && index >= 0) {
            result += (T)1.0;
        }

        current.x += device_lines[i].delta.x;
        current.y += device_lines[i].delta.y;
    }

    w_norms[i] = result;
}

template <typename T>
__global__ void sart_kernel(T* device_image, const device_line<T>* device_lines,
                            const T* device_sino, device_volume v, T beta,
                            const T* w_norms) {
    // gpu, first just 'closest' kernel (value 1)
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (w_norms[i] < EPSILON)
        return;

    T alpha = 0;
    vec2<T> current = device_lines[i].origin;
    while (current.x <= v.x + EPSILON && current.y <= v.y + EPSILON &&
           current.x + EPSILON >= 0 && current.y + EPSILON >= 0) {
        int index = (int)(current.x) + v.x * (int)(current.y);
        if (index < v.x * v.y && index >= 0) {
            alpha += device_image[index];
        }

        current.x += device_lines[i].delta.x;
        current.y += device_lines[i].delta.y;
    }

    auto factor = beta * (device_sino[i] - alpha) / w_norms[i];

    current = device_lines[i].origin;
    while (current.x <= v.x + EPSILON && current.y <= v.y + EPSILON &&
           current.x + EPSILON >= 0 && current.y + EPSILON >= 0) {
        int index = (int)(current.x) + v.x * (int)(current.y);
        if (index < v.x * v.y && index >= 0) {
            device_image[index] += factor;
        }

        current.x += device_lines[i].delta.x;
        current.y += device_lines[i].delta.y;
    }
}

template <typename T>
void run_sart(device_volume v, device_line<T>* device_lines, int lines,
              T* device_sino, T* host_image, int group_count, T beta = 0.5,
              int iterations = 10) {
    T* device_image = nullptr;
    auto image_bytes = v.x * v.y * sizeof(T);

    cudaMalloc(&device_image, image_bytes);
    cudaMemset(device_image, 0, image_bytes);

    int group_size = lines / group_count;
    int threads = 256;

    T* w_norms = nullptr;
    auto w_norms_bytes = lines * sizeof(T);
    cudaMalloc(&w_norms, w_norms_bytes);
    cudaMemset(w_norms, 0, w_norms_bytes);
    w_norms_kernel<<<lines / threads, threads>>>(device_lines, v, w_norms);

    cudaDeviceSynchronize();
    std::vector<T> w(lines);
    cudaMemcpy(w.data(), w_norms, w_norms_bytes, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    for (int i = 0; i < iterations; ++i) {
        for (int k = 0; k < group_count; ++k) {
            sart_kernel<<<group_size / threads, threads>>>(
                device_image, &device_lines[k * group_size],
                &device_sino[k * group_size], v, beta,
                &w_norms[k * group_size]);
        }
    }

    cudaMemcpy(host_image, device_image, image_bytes, cudaMemcpyDeviceToHost);

    cudaFree(w_norms);
    cudaFree(device_image);
}

template void run_sart(device_volume v, device_line<float>* device_lines,
                       int lines, float* device_sino, float* host_image, int,
                       float beta, int iterations);

template void run_sart(device_volume v, device_line<double>* device_lines,
                       int lines, double* device_sino, double* host_image,
                       int group_count, double beta, int iterations);

} // namespace cuda
} // namespace tomo
