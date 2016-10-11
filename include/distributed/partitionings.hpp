#pragma once

#include <cmath>
#include <limits>
#include <queue>
#include <stack>
#include <utility>
#include <vector>

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
class bisected_volume : public partitioned_volume<D> {
  public:
    struct split {
        // dimension of the split
        int d;

        // voxel index at which the *second* subvolume begins
        int a;
    };

    bisected_volume(tomo::volume<D> v, int processors)
        : partitioned_volume<D>(v, processors) {}

    int owner(int index) final override {
        auto voxel = this->v_.unroll(index);

        // we encode the path to the final volume as a bit pattern, that will be
        // the processor id
        auto node = splits_.root.get();
        int proc = 0;
        int depth = 0;
        while (node) {
            if (voxel[node->value.d] < node->value.a) {
                node = node->left.get();
            } else {
                proc += 1 << depth;
                node = node->right.get();
            }
            ++depth;
        }
        return proc;
    }

    core::binary_tree<split>& splits() { return splits_; }
    const core::binary_tree<split>& splits() const { return splits_; }

  private:
    core::binary_tree<split> splits_;
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

template <dimension D, typename Geometry,
          typename T = typename Geometry::value_type>
bisected_volume<D> partition_bisection(const Geometry& g, tomo::volume<D> v,
                                       int processors, T epsilon = 0.2) {
    auto result = bisected_volume<D>(v, processors);

    // assert that p is a power of two, see Hacker's Delight page 11
    assert((processors & (processors - 1)) == 0);

    // alias the split type of the form (d, a)
    using split = typename bisected_volume<D>::split;
    using box = std::array<math::vec2<int>, D>;

    // compute the 'depth' of the tree (i.e. log(p))
    auto depth = (int)log2(processors);
    assert(1 << depth == processors);

    auto find_split = [&](const std::vector<math::line<D, T>>& lines,
                          std::vector<math::line<D, T>>& lines_left,
                          std::vector<math::line<D, T>>& lines_right,
                          box bounds) -> split {
        (void)epsilon;
        (void)lines;
        (void)lines_left;
        (void)lines_right;

        // maintain two queues, or a single one?
        // auto q_enter = std::priority_queue<math::line<D, T>>();
        // auto q_exit = std::priority_queue<math::line<D, T>>();

        return split{0, bounds[0][0] + (bounds[0][1] - bounds[0][0]) / 2};
    };

    // containers for the current left and right lines
    // TODO: needed here?
    std::vector<math::line<D, T>> all_lines;
    for (auto l : g) {
        all_lines.emplace_back(l);
    }

    // we store the splits as a binary tree
    core::binary_tree<split> splits;

    using node = typename core::binary_tree<split>::node;
    using dir = typename core::binary_tree<split>::dir;

    // all the information we need to split a subvolume and save it
    struct subvolume {
        box bounds;
        node* parent;
        dir direction;
        std::vector<math::line<D, T>>& lines;
        int depth;
    };

    // state should be the correct lines and bounds
    std::stack<subvolume> split_stack;

    box bounds;
    for (int d = 0; d < D; ++d) {
        bounds[d][1] = v[d];
    }

    auto complete_volume = subvolume{bounds, nullptr, dir::left, all_lines, 0};
    split_stack.push(complete_volume);

    while (!split_stack.empty()) {
        const auto& sub = split_stack.top();

        int sub_depth = sub.depth;
        auto sub_bounds = sub.bounds;

        // containers to hold the left and right lines
        std::vector<math::line<D, T>> left;
        std::vector<math::line<D, T>> right;

        // we now have the subvolume, and want to find the best split
        auto best_split = find_split(sub.lines, left, right, sub.bounds);

        // add the split to the tree
        node* current_node =
            result.splits().add(sub.parent, sub.direction, best_split);

        split_stack.pop();

        (void)current_node;

        // now right and left contain lines for recursion
        // we move these into the 'lower splits'
        if (sub.depth + 1 < depth) {
            auto bounds_left = sub_bounds;
            bounds_left[best_split.d][1] = best_split.a;
            split_stack.push(subvolume{bounds_left, current_node, dir::left,
                                       left, sub_depth + 1});

            auto bounds_right = sub_bounds;
            bounds_right[best_split.d][0] = best_split.a;
            split_stack.push(subvolume{bounds_right, current_node, dir::right,
                                       right, sub_depth + 1});
        }
    }

    return result;
}

} // namespace distributed
} // namespace tomo
