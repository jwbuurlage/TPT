#pragma once

#include <limits>

#include "../volume.hpp"
#include "../projectors/closest.hpp"

namespace tomo {
namespace distributed {

/* what distributions do we want to support at first.. and how do we want to
 * measure */
template <dimension D>
class partitioned_volume {
  public:
    partitioned_volume(tomo::volume<D> v, int processors)
        : v_(v), processors_(processors) {}

    virtual int owner(int index) = 0;
    int processors() const { return processors_; }

    const tomo::volume<D>& volume() const { return v_; };

  protected:
    virtual ~partitioned_volume() = default;

    tomo::volume<D> v_;
    int processors_ = 0;
};

template <dimension D>
class slabbed_volume : public partitioned_volume<D> {
  public:
    slabbed_volume(tomo::volume<D> v, int processors, int dir)
        : partitioned_volume<D>(v, processors), dir_(dir) {
        int side = v[dir_];
        per_proc_ = (side - 1) / this->processors_ + 1;
    }

    int owner(int index) final override {
        auto voxel = this->v_.unroll(index);
        return voxel[dir_] / per_proc_;
    }

  private:
    int dir_ = 0;
    int per_proc_ = 0;
};

template <dimension D>
class erjinzhi_volume : public partitioned_volume<D> {
  public:
    erjinzhi_volume(tomo::volume<D> v, int processors)
        : partitioned_volume<D>(v, processors) {}

    int owner(int index) final override { return index % this->processors_; }

  private:
    // tree...
};

template <typename Geometry>
int overlap_count(Geometry& g, partitioned_volume<3_D>& v) {
    using T = typename Geometry::value_type;
    tomo::dim::closest<3_D, T> proj(v.volume());
    int overlap = 0;

    std::vector<int> seen(v.processors());
    for (auto line : g) {
        std::fill(seen.begin(), seen.end(), 0);
        for (auto& elem : proj(line)) {
            seen[v.owner(elem.index)] = 1;
        }
        int intersections =
            std::accumulate(seen.begin(), seen.end(), 0, std::plus<int>());
        if (intersections > 0)
            overlap += intersections - 1;
    }

    return overlap;
}

template <dimension D, typename Geometry>
slabbed_volume<D> partition_trivial(Geometry& g, tomo::volume<D> v,
                                    int processors) {

    auto best_overlap = std::numeric_limits<int>::max();
    auto vol = slabbed_volume<D>{v, processors, 0};

    for (int d = 0; d < D; ++d) {
        slabbed_volume<D> slabbed{v, processors, d};
        auto overlap = overlap_count(g, slabbed);
        if (overlap < best_overlap) {
            vol = slabbed;
            best_overlap = overlap;
        }
    }

    return vol;
}

template <dimension D, typename Geometry>
erjinzhi_volume<D> partition_smart(Geometry& g, tomo::volume<D> v,
                                   int processors) {
    (void)g;
    (void)v;
    (void)processors;
}

} // namespace distributed
} // namespace tomo
