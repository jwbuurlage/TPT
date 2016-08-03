#include <cmath>
#include <thrust/device_vector.h>

#include "fmt/format.h"
#include "hemi/parallel_for.h"

#include "device_data.hpp"
#include "projectors/closest.cu.hpp"

namespace tomo {
namespace cuda {

template <typename T>
__global__ void compute_sums(const device::line<T>* device_lines,
                             device::volume v, T* row, T* col) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    T partial_row = (T)0.0;
    project_closest(device_lines[i], v, [&partial_row, &col](int index) {
        partial_row += (T)1.0;
        col[index] += (T)1.0;
    });

    row[i] = partial_row;
}

template <typename T>
__global__ void invert_sums(T* row, T* col, T beta) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    col[i] = beta / col[i];
    row[i] = (T)1.0 / row[i];
}

template <typename T>
__global__ void sirt_fp(const device::line<T>* device_lines, T* sino_buffer,
                        const T* device_sino, device::volume v, const T* r,
                        cudaTextureObject_t image_texture) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    T result = 0;
    project_closest_texture(device_lines[i], v,
                            [&image_texture, &result](T x, T y) {
                                result += tex2D<T>(image_texture, x, y);
                            });

    sino_buffer[i] = (device_sino[i] - result) * r[i];
}

template <typename T>
__global__ void sirt_bp(cudaSurfaceObject_t image_buffer_surface,
                        const device::line<T>* device_lines, T* sino_buffer,
                        device::volume v) {
    int idx = blockDim.x * blockIdx.x + threadIdx.x;

    auto line_sum = sino_buffer[idx];
    project_closest_surface(
        device_lines[idx], v,
        [&line_sum, &image_buffer_surface, &idx](int i, int j) {
            T old = 0;
            surf2Dread(&old, image_buffer_surface, sizeof(T) * i, j);
            surf2Dwrite(old + line_sum, image_buffer_surface, sizeof(T) * i, j);
        });
}

template <typename T>
__global__ void sirt_scale(cudaSurfaceObject_t image_surface,
                           cudaTextureObject_t image_buffer_texture,
                           const T* col_sum, device::volume v) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int j = blockDim.y * blockIdx.y + threadIdx.y;
    T c = col_sum[i + j * v.x];
    surf2Dwrite(c * tex2D<T>(image_buffer_texture, i, j), image_surface,
                sizeof(T) * i, j);
}

template <typename T>
__global__ void sirt_clear(cudaSurfaceObject_t surface) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int j = blockDim.y * blockIdx.y + threadIdx.y;

    surf2Dwrite((T)0, surface, sizeof(T) * i, j);
}

template <typename T>
void run_sirt(cudaTextureObject_t image_texture,
              cudaSurfaceObject_t image_surface,
              cudaTextureObject_t image_buffer_texture,
              cudaSurfaceObject_t image_buffer_surface, device::volume v,
              device::line<T>* device_lines, int lines, T* device_sino,
              int group_count, T beta = 0.5, int iterations = 10) {
    int cells = v.x * v.y;

    constexpr int threads_x = 16;
    constexpr int threads_y = 16;
    constexpr int threads = threads_x * threads_y;

    // precompute R and beta * C
    thrust::device_vector<T> row_sums(lines);
    thrust::device_vector<T> relaxed_col_sums(cells);

    compute_sums<<<lines / threads, threads>>>(
        device_lines, v, row_sums.data().get(), relaxed_col_sums.data().get());

    invert_sums<<<lines / threads, threads>>>(
        row_sums.data().get(), relaxed_col_sums.data().get(), beta);

    thrust::device_vector<T> sino_buffer(lines);
    for (int i = 0; i < iterations; ++i) {
        // forward project
        sirt_fp<<<lines / threads, threads>>>(
            device_lines,             // device geometry
            sino_buffer.data().get(), // a buffer for the fp operation
            device_sino,              // the actual projection data
            v,                        // volume information
            row_sums.data().get(),    // the row sums
            image_texture);

        sirt_clear<T><<<dim3(v.x / threads_x, v.y / threads_y, 1),
                        dim3(threads_x, threads_y, 1)>>>(image_buffer_surface);

        // back project
        sirt_bp<<<lines / threads, threads>>>(
            image_buffer_surface,     // device memory for holding the image
            device_lines,             // device geometry
            sino_buffer.data().get(), // a buffer for the fp operation
            v                         // volume information
            );

        // scale, now with 2d blocks
        sirt_scale<<<dim3(v.x / threads_x, v.y / threads_y, 1),
                     dim3(threads_x, threads_y, 1)>>>(
            image_surface,                 // the partial image
            image_buffer_texture,          // the buffer
            relaxed_col_sums.data().get(), // the row sums
            v);
    }
}

template void run_sirt(cudaTextureObject_t, cudaSurfaceObject_t,
                       cudaTextureObject_t, cudaSurfaceObject_t, device::volume,
                       device::line<float>*, int, float* device_sino, int,
                       float, int);

} // namespace cuda
} // namespace tomo
