#include "tomo.hpp"

/**
 * TODO:
 * - Make a forward projector that works with distributed image instead of image
 *   > First generalize projectors to work with 'non centered volumes'
 *   > This is also in geometry
 *   > distributed image and distributed sinogram, or sinogram as coarray maybe?
 * - Add functions that manipulate in place
 * - Add functions for gradient (norm) computations
 * - Add functions for blurs
 * - Add GPU storage for partitioned array
 * - Generalize SIRT to use distributed image
 * - Optimize CUDA kernels for trajectory based geometries.
 * - How to 'check' and visualize result?
 * - Compare reconstruction with distributed ASTRA code.
 */

// // for reference, this is sequential code
// auto v = tomo::volume<2_D>(k, k);
// auto f = tomo::modified_shepp_logan_phantom<T>(v);
// auto g = tomo::geometry::parallel<2_D, T>(180, 250, v);
//
// auto proj = tomo::dim::closest<2_D, T>(v);
// auto sino = tomo::forward_projection<2_D, T>(f, g, proj);
//
// auto x = tomo::reconstruction::sirt(v, g, sino, beta, iterations);
// tomo::ascii_plot(x);


#include <array>
#include <cassert>
#include <cmath>

#include "bulk/bulk.hpp"
#include "set_backend.hpp"

template <typename T, int D, typename World>
class distributed_image {
   public:
    distributed_image(World& world, bulk::partitioning<D>& part)
        : data_(world, part) {}

   private:
    bulk::partitioned_array<T, D, World> data_;
};

template <typename Func, typename T, int D, typename World>
void transform_in_place(distributed_image<T, D, World>& img, Func func) {
    for (auto& elem : img) {
        elem = func(elem);
    }
}

// TODO border operations with ghost voxels, boundary on cpu rest on gpu?

int main() {
    bulk::environment<provider> env;

    env.spawn(env.available_processors(), [](auto world, int s, int p) {
        int N = sqrt(p);
        assert(N * N == p);

        std::array<int, 2> size{100, 100};
        auto block = bulk::block_partitioning<2, 2>({N, N}, size);

        // on the master node, distribute the geometry
        // 1. initialize the geometry and the imaged volume
        if (s == 0) {
            // 2. partition the volume using the geometry
            // 3. distribute the partitioned volume
        }

        // this should work for 'hyperrect partitionings'
        auto local_v = part_vol.local_bounds(s);

        // restrict the geometry to the local volume
        auto local_g = tomo::distributed::geometry(g, local_v);

        // now we should be ready to reconstruct
        auto f = tomo::distributed::localized_phantom(...);

        // now we should be ready to reconstruct
        auto proj = tomo::dim::closest<..>;

        // dont/do want distributed here..? how to sync the geometries and so on..
        auto local_sino = tomo::distributed::forward_projection<2_D, T>(f, g, proj);

        // reconstruction
        tomo::distributed::psirt(local_v, local_g, local_sino, beta, iterations);
    });
}
