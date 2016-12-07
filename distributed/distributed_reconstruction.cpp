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

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>

#include "bulk/backends/mpi/mpi.hpp"
#include "bulk/bulk.hpp"
using provider = bulk::mpi::provider;

int main() {
    using T = float;
    int k = 256;
    constexpr tomo::dimension D = 2;

    namespace td = tomo::distributed;

    bulk::environment<provider> env;

    env.spawn(env.available_processors(), [=](auto world, int /* s */, int p) {
        std::array<int, D> size{};
        std::fill(size.begin(), size.end(), k);

        auto block = bulk::block_partitioning<D, 1>(size, {p});

        // construct distributed image
        auto img = td::partitioned_image<T, D, decltype(world)>(world, block);

        // global_size why doesnt it match?
        auto global_volume = tomo::volume<D>(img.global_size());

        // non-centered (local) volumes can be initialized here
        // see img.local_volume().origin();
        auto v = img.local_volume();

        // support for distributed phantom construction
        td::partitioned_phantom(img);
        td::plot(img);

        // projectors and geometries are modified so that they are
        // intersected with volumes at proper location
        auto geom = tomo::geometry::parallel<D, T>(k, k, global_volume);
        auto proj = tomo::dim::closest<D, T>(v);

        // TODO restrict the geometry to the local volume
        //auto local_g = tomo::distributed::geometry(g, local_v);

        // The forward projection is modified so that we can perform it
        // in parallel on a distributed image, obtaining a 'distributed
        // sinogram'
        auto partitioned_sino =
            td::partitioned_sinogram<D, T, decltype(geom), decltype(world)>(
                world, block, geom);
        partitioned_sino.compute_overlap(proj);

        td::forward_project(img, geom, proj, partitioned_sino);

        world.sync();

        // [ ] 2. D-SIRT
        // TODO compute r and c
        std::vector<T> rs(geom.lines());
        std::vector<T> cs(geom.lines());

        int row = 0;
        for (auto l : geom) {
            for (auto elem : proj(l)) {
                rs[row] += elem.value;
                cs[elem.index] += elem.value;
            }
            row++;
        }

        // communicate row elements
        // note that this is exactly exchange in sino
        auto row_queue = bulk::create_queue<int, T>(world);
        for (auto exchange : partitioned_sino.exchanges()) {
            row_queue(exchange.target).send(exchange.line, rs[exchange.line]);
        }
        world.sync();

        for (auto msg : row_queue) {
            rs[msg.tag] += msg.content;
        }

        T beta = (T)0.5;

        // invert the sums
        for (auto& r : rs)
            r = 1.0 / r;
        for (auto& c : cs)
            c = beta / c;

        // temporary sino
        auto tmp_sino = td::partitioned_sinogram<D, T, decltype(geom), decltype(world)>(
                world, block, geom);
        // TODO construct using already computed exchanges
        tmp_sino.compute_overlap(proj);

        auto x = td::partitioned_image<T, D, decltype(world)>(world, block);
        auto tmp_image = td::partitioned_image<T, D, decltype(world)>(world, block);

        for (int iter = 0; iter < 100; ++iter) {
            tmp_sino.clear();

            // forward project x
            td::forward_project(x, geom, proj, tmp_sino);

            for (int j = 0; j < geom.lines(); ++j) {
                tmp_sino[j] = rs[j] * (partitioned_sino[j] - tmp_sino[j]);
            }

            tmp_image.clear();
            td::back_project(tmp_image, geom, proj, tmp_sino);

            for (int j = 0; j < v.cells(); ++j) {
                x[j] += cs[j] * tmp_image[j];
            }
        }

        td::plot(x);

        // [ ] 3. Support for boxing and so on

        // [ ] 4. Non-local (boundary) operations

        // [ ] 5. CUDA storage backend here
    });
}
