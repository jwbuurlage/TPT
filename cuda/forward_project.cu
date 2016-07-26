#include "hemi/parallel_for.h"

#include <iostream>

namespace tomo {
namespace cuda {

template <typename T>
void internal_cuda_stuff(T* image_data, T* sinogram_data /* ... */) {

}

void test() {
    hemi::launch([=] HEMI_LAMBDA() {
        printf("Hello World from Lambda in thread %d of %d\n",
               hemi::globalThreadIndex(), hemi::globalThreadCount());
    });

    hemi::deviceSynchronize();
}

} // namespace cuda
} // namespace tomo
