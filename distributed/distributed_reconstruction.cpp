#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>

#include "bulk/backends/mpi/mpi.hpp"
#include "bulk/bulk.hpp"

#include "tomos/distributed/plotter.hpp"
#include "tomos/tomos.hpp"
#include "tomos/util/trees.hpp"

namespace td = tomo::distributed;

using T = float;
constexpr tomo::dimension G = 1;

template <tomo::dimension D>
void run(tomo::util::args opt) {
    bulk::mpi::environment env;
    env.spawn(env.available_processors(), [=](auto& world) {
        auto s = world.processor_id();
        auto p = world.active_processors();

        td::ext_plotter<D, T> plotter("tcp://localhost:5555",
                                      "Distributed reconstruction", world);
        auto bench = tomo::benchmark("reconstruction");
        if (s != 0)
            bench.silence();

        bench.phase("initialize");

        auto global_volume = tomo::volume<D, T>(opt.k);

        // projectors and geometries are modified so that they are
        // intersected with volumes at proper location
        auto geom = tomo::geometry::parallel<D, T>(global_volume, opt.k, opt.k);

        // set up the partitioning
        std::array<int, D> size{};
        std::fill(size.begin(), size.end(), opt.k);
        // auto partitioning = td::partition_trivial(geom, global_volume, p);
        // auto partitioning = bulk::block_partitioning<D, G>(size, {p}, {2});
        auto partitioning_ptr = tomo::load_partitioning<T>(
            "data/partitionings/cone_test.toml", global_volume);
        auto& partitioning = *partitioning_ptr;

        // construct distributed image and the associated volumes
        auto img =
            td::partitioned_image<D, G, T>(world, partitioning, global_volume);
        auto local_volume = img.local_volume();

        // we initialize a Shepp-Logan phantom inside the image
        td::partitioned_phantom(img);
        plotter.plot(img);

        // auto geom = tomo::geometry::cone_beam<T>(global_volume, opt.k, 1.5,
        // {opt.k, opt.k});
        auto proj = tomo::dim::linear<D, T>(local_volume);

        // the forward projection is modified so that we can perform it
        // in parallel on a distributed image, obtaining a 'distributed
        // projection_stack'
        auto partitioned_ps =
            td::partitioned_projection_stack<D, T>(world, partitioning, geom);

        // we compute the overlap on the projection stack
        bench.phase("compute overlap");
        partitioned_ps.compute_overlap(proj);
        world.log("com vol = %i / %i\n", partitioned_ps.communication_volume(),
                  geom.lines() * (p - 1));

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

        for (auto [tag, content] : row_queue) {
            rs[tag] += content;
        }

        T beta = (T)1.0;

        // invert the sums
        for (auto& r : rs)
            r = 1.0 / r;
        for (auto& c : cs)
            c = (tomo::math::abs(beta) > tomo::math::epsilon<T>) ? (beta / c)
                                                                 : (T)0;

        bench.phase("initialize sirt");
        // temporary ps
        auto buffer_ps =
            td::partitioned_projection_stack<D, T>(world, partitioning, geom);
        // TODO construct using already computed exchanges
        buffer_ps.compute_overlap(proj);

        auto x =
            td::partitioned_image<D, G, T>(world, partitioning, global_volume);
        auto buffer_image =
            td::partitioned_image<D, G, T>(world, partitioning, global_volume);

        bench.phase("iterating sirt");
        for (int iter = 0; iter < opt.iterations; ++iter) {
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

            // plotter.plot(x);
        }
        plotter.plot(x);
    });
}

int main(int argc, char* argv[]) {
    auto opt = tomo::util::args(argc, argv);
    run<3_D>(opt);

    return 0;
}
