#include "distributed/plotter.hpp"
#include "tomo.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>

#include "bulk/backends/mpi/mpi.hpp"
#include "bulk/bulk.hpp"

namespace td = tomo::distributed;

int main() {
    using T = float;
    int k = 512;
    constexpr tomo::dimension D = 2;

    bulk::mpi::environment env;
    env.spawn(env.available_processors(), [=](auto& world, int s, int p) {
        td::ext_plotter<D, T> plotter("tcp://localhost:5555",
                                      "Distributed reconstruction", world);

        auto bench = tomo::benchmark("reconstruction");
        if (s != 0)
            bench.silence();

        bench.phase("initialize");

        // set up the partitioning
        std::array<int, D> size{};
        std::fill(size.begin(), size.end(), k);
        auto block = bulk::block_partitioning<D, 2>(world, size, {2, p / 2});

        // construct distributed image and the associated volumes
        auto img = td::partitioned_image<D, 2, T>(world, block);
        auto global_volume = tomo::volume<D>(img.global_size());
        auto local_volume = img.local_volume();

        // we initialize a Shepp-Logan phantom inside the image
        td::partitioned_phantom(img);
        plotter.plot(img);

        // projectors and geometries are modified so that they are
        // intersected with volumes at proper location
        auto geom = tomo::geometry::parallel<D, T>(k, k, global_volume);
        auto proj = tomo::dim::closest<D, T>(local_volume);

        // the forward projection is modified so that we can perform it
        // in parallel on a distributed image, obtaining a 'distributed
        // projection_stack'
        auto partitioned_ps =
            td::partitioned_projection_stack<D, T, decltype(geom)>(world, block,
                                                                   geom);

        // we compute the overlap on the projection stack
        bench.phase("compute overlap");
        partitioned_ps.compute_overlap(proj);
        world.log("com vol = %i\n", partitioned_ps.communication_volume());

        // we perform one forward projection
        bench.phase("one fp (experiment)");
        td::forward_project(img, geom, proj, partitioned_ps);

        bench.phase("compute rs, cs");
        std::vector<T> rs(geom.lines());
        std::vector<T> cs(local_volume.cells());

        int row = 0;
        for (auto l : geom) {
            for (auto elem : proj(l)) {
                rs[row] += elem.value;
                cs[elem.index] += elem.value;
            }
            row++;
        }

        // communicate row elements, note that this corresponds to the computed
        // overlap
        auto row_queue = bulk::queue<int, T>(world);
        for (auto exchange : partitioned_ps.exchanges()) {
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

        bench.phase("initialize sirt");
        // temporary ps
        auto buffer_ps =
            td::partitioned_projection_stack<D, T, decltype(geom)>(world, block,
                                                                   geom);
        // TODO construct using already computed exchanges
        buffer_ps.compute_overlap(proj);

        auto x = td::partitioned_image<D, D, T>(world, block);
        auto buffer_image = td::partitioned_image<D, D, T>(world, block);

        bench.phase("10 times sirt");
        for (int iter = 0; iter < 10; ++iter) {
            buffer_ps.clear();

            // forward project x
            td::forward_project(x, geom, proj, buffer_ps);

            for (int j = 0; j < geom.lines(); ++j) {
                buffer_ps[j] = rs[j] * (partitioned_ps[j] - buffer_ps[j]);
            }

            buffer_image.clear();
            td::back_project(buffer_image, geom, proj, buffer_ps);

            for (int j = 0; j < local_volume.cells(); ++j) {
                x[j] += cs[j] * buffer_image[j];
            }

            plotter.plot(x);
        }
    });
}
