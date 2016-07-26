#include "tomo.hpp"

namespace tomo {
namespace cuda {

extern void internal_cuda_stuff(...);

// front end for cuda, calling 'internal stuff' for each geometry etc.

struct external_cuda_projector {};

template <dimension D, typename T, class Geometry>
sinogram<D, T, Geometry, external_cuda_projector>
forward_projection(const image<D, T>& f, const Geometry& g, Projector& proj) {
    sinogram<D, T, Geometry, external_cuda_projector> result(g);
    // 1. copy data
    // 2. run alg
    // 3. copy data back
    // 4. return sino
    return result;
}

} // namespace cuda
} // namespace tomo
