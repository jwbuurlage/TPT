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

template <tomo::dimension D>
auto to_vec(std::array<int, D> in) {
    tomo::math::vec<D, int> out;
    for (int d = 0; d < D; ++d) {
        out[d] = in[d];
    }
    return out;
}

template <typename T, tomo::dimension D, typename World, typename Part>
class partitioned_image {
  public:
    partitioned_image(World& world, Part& part)
        : data_(world, part), world_(world), part_(part),
          multi_index_(bulk::unflatten<D>(part_.grid(), world.processor_id())),
          local_volume_(to_vec<D>(part_.origin(multi_index_)),
                        to_vec<D>(part_.local_extent(multi_index_))) {}

    auto world() const { return world_; }

    tomo::math::vec<D, int> local_size() const {
        auto size = part_.local_extent(multi_index_);
        tomo::math::vec<D, int> result;
        for (int d = 0; d < D; ++d) {
            result[d] = size[d];
        }
        return result;
    }

    auto global_size() const { return to_vec<D>(part_.global_size()); }
    auto local_volume() const { return local_volume_; }

    T& operator[](int i) {
        auto idx = bulk::unflatten<D>(part_.global_size(), i);
        return data_.local(idx[0], idx[1]);
    }

    T& operator()(std::array<int, D> xs) {
        // FIXME this is unnecessarily involved because of the
        // std::array/tomo::math::vec discrepency. Come up with a solution.
        return (*this)[bulk::flatten<D>(part_.local_extent(multi_index_), xs)];
    }

  private:
    bulk::partitioned_array<T, D, World> data_;
    World& world_;
    Part& part_;

    std::array<int, D> multi_index_;
    tomo::volume<D> local_volume_;
};

template <typename Func, typename T, int D, typename World, typename Part>
void transform_in_place(partitioned_image<T, D, World, Part>& img, Func func) {
    for (auto& elem : img) {
        elem = func(elem);
    }
}

template <typename T, typename World, typename Part>
void partitioned_phantom(partitioned_image<T, 2, World, Part>& img) {
    auto global_volume =
        tomo::volume<2>(img.global_size()[0], img.global_size()[1]);
    tomo::fill_ellipses_(img, tomo::mshl_ellipses_<T>(), img.local_volume(),
                         global_volume);
}

template <typename World, typename Func>
void output_in_turn(World& world, Func func) {
    // easiest is to copy the local elements into a tomo::image, then output the
    // pid, origin, extent, and finally plot it.
    // we ensure correct oder here
    for (int t = 0; t < world.active_processors(); ++t) {
        if (world.processor_id() == t) {
            func();
        }
        sleep(1);
        world.sync();
    }
}

template <typename T, typename World, typename Part>
void plot(partitioned_image<T, 2_D, World, Part>& img) {
    auto world = img.world();
    auto volume = img.local_volume();
    auto cells = volume.cells();

    T max = std::numeric_limits<T>::min();

    auto local_image = tomo::image<2_D, T>(img.local_volume());
    for (int i = 0; i < cells; ++i) {
        local_image[i] = img[i];
        if (max < img[i]) {
            max = img[i];
        }
    }

    // exchange global maximum
    auto maxs = bulk::gather_all(world, max);
    max = *std::max_element(maxs.begin(), maxs.end());

    output_in_turn(world, [&]() {
        world.log("(%i x %i) + (%i, %i) \n", img.local_size()[0],
                  img.local_size()[1], volume.origin()[0], volume.origin()[1]);
        tomo::ascii_plot<T>(local_image, max);
    });
}

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
        auto partitioned_sino =
            tomo::distributed::partitioned_sinogram<2_D, T, decltype(geom),
                                                    decltype(proj)>(geom);
        tomo::distributed::forward_project<2_D, T>(img, geom, proj, partitioned_sino);

        // TODO define distributed sinogram
        // restricting geometry
        // restricting image
        // TODO: exchange all sinograms, we shall call this 'harmonizing'

        output_in_turn(world, [&]() { tomo::ascii_plot(sino); });
        sino.harmonize();

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
