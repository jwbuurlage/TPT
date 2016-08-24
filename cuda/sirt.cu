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
__global__ void sirt_forward_project(const device::line<T>* device_lines,
                                     T* sino_buffer, const T* device_sino,
                                     device::volume v, const T* r,
                                     cudaTextureObject_t image_texture) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    T result = 0;
    project_closest_texture(device_lines[i], v,
                            [&image_texture, &result](T x, T y) {
                                result += tex2D<T>(image_texture, x, y);
                            });

    sino_buffer[i] = (device_sino[i] - result) * r[i];
}

// FIXME rewrite this to be 'voxel-centric'
// - we have 48 kb per block, if we want to assign subvolume to block
// 12 000 floats, +- 100 x 100 size subvolumes, could write atomically to
// shared memory, and then once atomically to global memory
// - could also do the step idea
// - slow solution is to have another buffer we can write atomically to (current
// impl)
template <typename T>
__global__ void sirt_back_project(const device::line<T>* device_lines,
                                  T* sino_buffer, device::volume v,
                                  T* image_second_buffer) {
    // one thread per line..
    int idx = blockDim.x * blockIdx.x + threadIdx.x;

    T line_sum = sino_buffer[idx];
    project_closest(device_lines[idx], v,
                    [&line_sum, &image_second_buffer](int i) {
                        atomicAdd(&image_second_buffer[i], line_sum);
                    });
}

template <typename T>
__global__ void sirt_scale(cudaSurfaceObject_t image_surface,
                           const T* image_buffer, const T* col_sum,
                           device::volume v) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int j = blockDim.y * blockIdx.y + threadIdx.y;
    int k = i + j * v.x;

    T old = 0;
    surf2Dread(&old, image_surface, sizeof(T) * i, j);
    surf2Dwrite(old + col_sum[k] * image_buffer[k], image_surface,
                sizeof(T) * i, j);
}

template <typename T>
void run_sirt(cudaTextureObject_t image_texture,
              cudaSurfaceObject_t image_surface, device::volume v,
              device::line<T>* device_lines, int lines, T* device_sino,
              int group_count, T beta = 0.5, int iterations = 10) {
    int cells = v.x * v.y;

    constexpr int threads_x = 16;
    constexpr int threads_y = 16;
    constexpr int threads = threads_x * threads_y;

    // precompute R and beta * C
    thrust::device_vector<T> row_sums(lines);
    thrust::device_vector<T> relaxed_col_sums(cells);

    thrust::device_vector<T> image_second_buffer(cells);

    compute_sums<<<lines / threads, threads>>>(
        device_lines, v, row_sums.data().get(), relaxed_col_sums.data().get());

    invert_sums<<<lines / threads, threads>>>(
        row_sums.data().get(), relaxed_col_sums.data().get(), beta);

    thrust::device_vector<T> sino_buffer(lines);
    for (int i = 0; i < iterations; ++i) {
        // forward project
        sirt_forward_project<<<lines / threads, threads>>>(
            device_lines,             // device geometry
            sino_buffer.data().get(), // a buffer for the fp operation
            device_sino,              // the actual projection data
            v,                        // volume information
            row_sums.data().get(),    // the row sums
            image_texture);

        thrust::fill(image_second_buffer.begin(), image_second_buffer.end(), 0);

        // back project
        sirt_back_project<<<lines / threads, threads>>>(
            device_lines,             // device geometry
            sino_buffer.data().get(), // data for the bp operation
            v,                        // volume info
            image_second_buffer.data().get());

        // scale, now with 2d blocks
        sirt_scale<<<dim3(v.x / threads_x, v.y / threads_y, 1),
                     dim3(threads_x, threads_y, 1)>>>(
            image_surface,                    // the partial image
            image_second_buffer.data().get(), // the buffer
            relaxed_col_sums.data().get(),    // the row sums
            v);
    }
}

template void run_sirt(cudaTextureObject_t, cudaSurfaceObject_t, device::volume,
                       device::line<float>*, int, float* device_sino, int,
                       float, int);

} // namespace cuda
} // namespace tomo
