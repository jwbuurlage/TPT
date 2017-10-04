#include "tomos/tomos.hpp"
using namespace tomo;

int main() {
    using T = float;

    constexpr dimension D = 3_D;

    int size = 32;
    int sweeps = 10;
    auto volume = tomo::volume<D, T>(size);
    auto phantom = tomo::modified_shepp_logan_phantom<T>(volume);
    tomo::ascii_plot(phantom);

    auto geometry = tomo::geometry::cone_beam<T>(volume, size, {(T)2, (T)2},
                                                 {size, size}, 2.0, 2.0);
    auto kernel = tomo::dim::joseph<D, T>(volume);
    auto b = tomo::forward_projection(phantom, geometry, kernel);

    auto cs = tomo::image<D, T>(volume);
    for (auto j = 0u; j < volume.cells(); ++j) {
        auto voxel_index = volume.unroll(j);
        auto voxel_volume = tomo::volume<D, T>(
            {1, 1, 1},
            math::vec<D, T>(voxel_index) / math::vec<D, T>(volume.voxels()),
            math::vec<D, T>(1.0) / math::vec<D, T>(volume.voxels()));
        auto voxel_rays =
            tomo::distributed::restricted_geometry(geometry, voxel_volume);
        for (auto[proj, local_i, ray] : voxel_rays) {
            (void)proj;
            (void)local_i;
            auto value = kernel.matrix_value(ray, voxel_index);
            cs[j] += value * value;
        }
    }

    auto x = tomo::image<D, T>(volume);
    auto r = b;

    // sequential cyclic
    for (auto k = 0; k < sweeps; ++k) {
        for (auto j = 0u; j < volume.cells(); ++j) {
            if (cs[j] < tomo::math::epsilon<T>) {
                continue;
            }
            auto voxel_index = volume.unroll(j);
            auto voxel_volume = tomo::volume<D, T>(
                {1, 1, 1},
                math::vec<D, T>(voxel_index) / math::vec<D, T>(volume.voxels()),
                math::vec<D, T>(1.0) / math::vec<D, T>(volume.voxels()));
            auto voxel_rays =
                tomo::distributed::restricted_geometry(geometry, voxel_volume);
            auto delta = (T)0;
            for (auto[proj, local_i, ray] : voxel_rays) {
                auto i = voxel_rays.global_index(proj, local_i);
                auto value = kernel.matrix_value(ray, voxel_index);
                delta += value * r[i];
            }
            delta /= cs[j];
            for (auto[proj, local_i, ray] : voxel_rays) {
                auto i = voxel_rays.global_index(proj, local_i);
                auto value = kernel.matrix_value(ray, voxel_index);
                r[i] -= value * delta;
            }
            x[j] += delta;
        }
    }

    // print result
    tomo::ascii_plot(x);
    tomo::write_png(tomo::slice<T>(x, size / 2, 2), "central_slice");
}
