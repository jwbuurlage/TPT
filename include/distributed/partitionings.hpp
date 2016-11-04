#pragma once

#include <cmath>
#include <cstddef>
#include <experimental/optional>
#include <limits>
#include <queue>
#include <set>
#include <stack>
#include <utility>
#include <vector>
#include <sstream>

template <typename T>
using optional = std::experimental::optional<T>;

#include "../common.hpp"
#include "../math.hpp"
#include "../projectors/closest.hpp"
#include "../volume.hpp"

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

    virtual ~partitioned_volume() = default;

  protected:
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

    int axis() const { return dir_; }
    int slab_height() const { return per_proc_; }

  private:
    int dir_ = 0;
    int per_proc_ = 0;
};

template <typename Dim, typename Geometry>
int overlap_count(Geometry& g, partitioned_volume<3_D>& v) {
    Dim proj(v.volume());

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

template <typename Dim, dimension D, typename Geometry>
slabbed_volume<D> partition_trivial(Geometry& g, tomo::volume<D> v,
                                    int processors) {

    auto best_overlap = std::numeric_limits<int>::max();
    auto vol = slabbed_volume<D>{v, processors, 0};

    for (int d = 0; d < D; ++d) {
        slabbed_volume<D> slabbed{v, processors, d};
        auto overlap = overlap_count<Dim>(g, slabbed);
        if (overlap < best_overlap) {
            vol = slabbed;
            best_overlap = overlap;
        }
    }

    return vol;
}

} // namespace distributed
} // namespace tomo
