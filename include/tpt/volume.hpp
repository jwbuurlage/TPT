#pragma once

#include <cstdint>
#include <numeric>

#include "common.hpp"
#include "math/vector.hpp"
#include "math/vector_operations.hpp"

namespace tpt {

/**
 * The region which is being imaged.
 *
 * The volume describes the resolution of the reconstructed image. This is used
 * to properly construct the geometry of the experiment, but it is also used in
 * the descritization and reconstruction methods themselves.
 *
 * \tparam D the dimension of the volume (and thus the reconstruction problem).
 * \tparam T the scalar type to use
*/
template <dimension D, typename T>
class volume {
  public:
    /**
     * Construct a (hyper)rectangular volume that has an offset.
     *
     * \param origin the origin of the rectangular volume
     * \param lengths an array of which the i-th element corresponds to
     * the spanning width of the volume on the i-th axis.
     */
    volume(math::vec<D, int> voxels, math::vec<D, T> origin,
           math::vec<D, T> lengths)
        : voxels_(voxels), lengths_(lengths), origin_(origin) {}

    /**
     *  Construct a unit cubic volume spanning `voxels` voxels on each axis,
     * that
     *  sits on the (mathematical) origin.
     */
    volume(math::vec<D, int> voxels)
        : volume(voxels, math::vec<D, T>((T)0.0), math::vec<D, T>((T)1.0)) {}

    /// Construct a unit cubic volume spanning `k` voxels on each axis, that
    /// sits on the (mathematical) origin.
    volume(int k) : volume(math::vec<D, int>(k)) {}

    /**
     * Obtain the physical size in the given axis
     *
     * \param i the index of the dimension
     * \returns size of the volume in the axis
     */
    T operator[](size_t i) const { return lengths_[i]; }

    /**
     * Obtain the index corresponding to a voxel
     *
     * \param xs a vector describing the voxel
     * \returns the index of the voxel
     */
    template <typename Vector>
    int index_by_vector(Vector xs) const {
        return index_by_vector_(xs);
    }

    /**
     * Obtain the index corresponding to a voxel
     *
     * \param xs an array describing the voxel
     * \returns the index of the voxel
     */
    int index(math::vec<D, int> xs) const { return index_by_vector_(xs); }

    /**
     * Obtain the index corresponding to a voxel
     *
     * The i-th parameter corresponds to the index of the voxel in the i-th
     * dimension.
     *
     * \returns the index of the voxel
     */
    template <typename... Ss, typename = check_dim<D, Ss...>>
    int index(Ss... xs) const {
        return index_(0, 1, xs...);
    }

    /**
     * Obtain the voxels of the volume.
     *
     * \returns an array whose i-th element corresponds to the number of voxels
     * of the volume in the i-th axis.
     */
    math::vec<D, int> voxels() const { return voxels_; }

    /**
     * Obtain the origin of the volume.
     *
     * \returns an array whose i-th element corresponds to the spanning width
     * of the volume in the i-th axis.
     */
    math::vec<D, T> origin() const { return origin_; }

    /**
     * Obtain the lengths of the sides of the volume.
     *
     * \returns a vector whose i-th element corresponds to the physical spanning
     * width
     * of the volume in the i-th axis.
     */
    math::vec<D, T> physical_lengths() const { return lengths_; }

    /**
     * Obtain the total number of cells (voxels) in the volume.
     *
     * \returns number of voxels in the volume.
     */
    uint64_t cells() const { return math::product<D, int>(voxels_); }

    /** Unroll on index, i.e. obtain the multi-index. */
    math::vec<D, int> unroll(int idx) {
        math::vec<D, int> cell;
        for (int d = 0; d < D; ++d) {
            cell[d] = idx % voxels_[d];
            idx /= voxels_[d];
        }
        return cell;
    }

    void set_voxels(math::vec<D, int> new_voxels) { voxels_ = new_voxels; }

    /**
     * Obtain a std::vector of D-dimensional vectors containing the corners of
     * the volume.
     *
     * Note: only implemented for D = 3.
     */
    std::vector<math::vec<D, T>> corners() {
        static_assert(D == 3_D, "volume corners only implemented for D = 3");
        auto axes = std::array<math::vec3<T>, 3>{
            lengths_[0] * math::standard_basis<3_D, T>(0),
            lengths_[1] * math::standard_basis<3_D, T>(1),
            lengths_[2] * math::standard_basis<3_D, T>(2)};
        return {origin_,
                origin_ + axes[0],
                origin_ + axes[1],
                origin_ + axes[0] + axes[1],
                origin_ + axes[2],
                origin_ + axes[0] + axes[2],
                origin_ + axes[1] + axes[2],
                origin_ + axes[0] + axes[1] + axes[2]};
    }

  private:
    template <typename S, typename... Ss>
    int index_(int current, int offset, S x, Ss... xs) const {
        current += offset * x;
        offset *= voxels_[D - 1 - sizeof...(xs)];
        return index_(current, offset, xs...);
    }

    template <typename Vector>
    inline int index_by_vector_(Vector xs) const {
        int result = xs[0];
        int offset = voxels_[0];
        for (int i = 1; i < D; ++i) {
            result += offset * xs[i];
            offset *= voxels_[i];
        }
        return result;
    }

    int index_(int current, int /* offset */) const { return current; }

    math::vec<D, int> voxels_;
    math::vec<D, T> lengths_;
    math::vec<D, T> origin_;
};

} // namespace tpt
