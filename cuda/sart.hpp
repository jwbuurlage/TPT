#include "tomo.hpp"

#include <cuda_runtime.h>

#include "device_data.hpp"

namespace tomo {
namespace cuda {

template <typename T>
void run_sart(device_volume v, device_line<T>* device_lines, int lines,
              T* device_sino, T* host_image, int group_count, T beta = 0.5,
              int iterations = 10);

extern template void run_sart(device_volume v, device_line<float>* device_lines,
                              int lines, float* device_sino, float* host_image,
                              int group_count, float beta, int iterations);
extern template void run_sart(device_volume v,
                              device_line<double>* device_lines, int lines,
                              double* device_sino, double* host_image,
                              int group_count, double beta, int iterations);

template <dimension D, typename T, class Geometry, class Projector>
image<D, T> sart(const volume<D>& v, const Geometry& g,
                 const sinogram<D, T, Geometry, Projector>& p,
                 T beta = (T)0.5, int iterations = 10) {
    image<D, T> f(v);

    // FIXME: this duplicates FP code, would be really cool to make abstractions
    // for this. See Thrust as well.

    // 1. copy data
    std::vector<device_line<T>> lines(g.lines());
    int k = 0;
    for (auto line : g) {
        lines[k].origin.x = line.origin.x;
        lines[k].origin.y = line.origin.y;
        lines[k].delta.x = line.delta.x;
        lines[k].delta.y = line.delta.y;
        k++;
    }

    device_volume device_v = {f.get_volume().x(), f.get_volume().y()};

    device_line<T>* device_lines = nullptr;
    T* device_sino = nullptr;
    auto line_bytes = g.lines() * sizeof(device_line<T>);
    auto sino_bytes = g.lines() * sizeof(T);
    cudaMalloc(&device_lines, line_bytes);
    cudaMalloc(&device_sino, sino_bytes);

    cudaMemcpy(device_sino, p.data().data(), sino_bytes,
               cudaMemcpyHostToDevice);
    cudaMemcpy(device_lines, lines.data(), line_bytes, cudaMemcpyHostToDevice);

    // 2. run alg
    run_sart(device_v, device_lines, g.lines(), device_sino,
             f.mutable_data().data(), g.groups()[0], beta, iterations);

    // 3. free memory
    cudaFree(device_lines);
    cudaFree(device_sino);

    return f;
}

} // namespace cuda
} // namespace tomo
