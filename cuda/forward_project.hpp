#include "tomo.hpp"

#include <cuda_runtime.h>

#include "device_data.hpp"

namespace tomo {
namespace cuda {

extern void test();

template <typename T>
void run_forward_projection(T* device_image, device_line<T>* device_lines,
                            int lines, T* result, device_volume);

extern template void run_forward_projection<float>(float*, device_line<float>*,
                                                   int, float*, device_volume);
extern template void run_forward_projection<double>(double*,
                                                    device_line<double>*, int,
                                                    double*, device_volume);

struct external_cuda_projector {};

template <dimension D, typename T, class Geometry>
sinogram<D, T, Geometry, external_cuda_projector>
forward_projection(const image<D, T>& f, const Geometry& g) {
    sinogram<D, T, Geometry, external_cuda_projector> result(g);
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

    device_volume v = {f.get_volume().x(), f.get_volume().y()};
    T* device_image = nullptr;
    device_line<T>* device_lines = nullptr;
    auto image_bytes = f.data().size() * sizeof(T);
    auto line_bytes = g.lines() * sizeof(device_line<T>);
    cudaMalloc(&device_image, image_bytes);
    cudaMalloc(&device_lines, line_bytes);

    cudaMemcpy(device_image, f.data().data(), image_bytes,
               cudaMemcpyHostToDevice);
    cudaMemcpy(device_lines, lines.data(), line_bytes, cudaMemcpyHostToDevice);

    // 2. run alg
    run_forward_projection(device_image, device_lines, g.lines(),
                           result.mutable_data().data(), v);



    // 3. free data
    cudaFree(device_image);
    cudaFree(device_lines);

    return result;
}

} // namespace cuda
} // namespace tomo
