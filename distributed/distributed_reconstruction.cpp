#include "tomo.hpp"

/**
 * TODO:
 * [ ] How to 'check' and visualize result?
 * [ ] Compare reconstruction with distributed ASTRA code.
 */
// TODO border operations with ghost voxels, boundary on cpu rest on gpu?

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
#include <limits>

#include "bulk/backends/mpi/mpi.hpp"
#include "bulk/bulk.hpp"
using provider = bulk::mpi::provider;

int main() {
    using T = float;
    int k = 64;

    bulk::environment<provider> env;

    env.spawn(env.available_processors(), [=](auto world, int /* s */, int p) {
        int N = sqrt(p);
        assert(N * N == p);

        std::array<int, 2> size{k, k};
       auto block = bulk::block_partitioning<2, 2>({N, N}, size);

        // [x] 1A. construct distributed image
        auto img = partitioned_image<T, 2, decltype(world), decltype(block)>(
            world, block);

        auto global_volume =
            tomo::volume<2>(img.global_size()[0], img.global_size()[1]);

        // [x] 1B. non-centered (local) volumes can be initialized here
        // see img.local_volume().origin();
        auto v = img.local_volume();

        // [x] 1C. support for distributed phantom construction
        partitioned_phantom(img);
        plot(img);

        // [x] 1D. projectors and geometries are modified so that they are
        // intersected with volumes at proper location
        auto geom = tomo::geometry::parallel<2_D, T>(k, k, global_volume);
        auto proj = tomo::dim::closest<2_D, T>(v);

        // [x] 1E. The forward projection is modified so that we can perform it
        // in
        // parallel on a distributed image, obtaining a 'distributed sinogram'

        // [ ] 2. D-SIRT

        // [ ] 3. Support for boxing and so on

        // [ ] 4. Non-local (boundary) operations

        // [ ] 5. CUDA storage backend here

        /*
        // this should work for 'hyperrect partitionings'
         auto local_v = part_vol.local_bounds(s);

         // restrict the geometry to the local volume
         auto local_g = tomo::distributed::geometry(g, local_v);

         // now we should be ready to reconstruct
         auto f = tomo::distributed::localized_phantom(...);

         // now we should be ready to reconstruct
         auto proj = tomo::dim::closest<..>;

         // dont/do want distributed here..? how to sync the geometries and so
        on..
         auto local_sino = tomo::distributed::forward_projection<2_D, T>(f, g,
        proj);

         // reconstruction
         tomo::distributed::psirt(local_v, local_g, local_sino, beta,
        iterations); */
    });
}
