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

template <typename T>
using optional = std::experimental::optional<T>;

#include "../common.hpp"
#include "../math.hpp"
#include "../math/stringify.hpp"
#include "../projectors/closest.hpp"
#include "../volume.hpp"

#include <bulk/bulk.hpp>

namespace tomo {
namespace distributed {

// alias the split type of the form (d, a)
template <dimension D>
using box_t = std::array<math::vec2<int>, D>;

using split_t = bulk::util::split;

template <dimension D, typename T>
image<D, T> partial_sums(const image<D, T>& w) {
    auto w_sums = w;
    auto v = w.get_volume();

    auto idx = [&](auto x, auto y, auto z) { return w_sums.index({x, y, z}); };

    // 2. Compute 3D partial sums everywhere
    for (int k = 0; k < v.voxels()[2]; ++k) {
        for (int j = 0; j < v.voxels()[1]; ++j) {
            for (int i = 1; i < v.voxels()[0]; ++i) {
                w_sums[idx(i, j, k)] += w_sums[idx(i - 1, j, k)];
            }
        }
    }

    for (int k = 0; k < v.voxels()[2]; ++k) {
        for (int j = 1; j < v.voxels()[1]; ++j) {
            for (int i = 0; i < v.voxels()[0]; ++i) {
                w_sums[idx(i, j, k)] += w_sums[idx(i, j - 1, k)];
            }
        }
    }

    for (int k = 1; k < v.voxels()[2]; ++k) {
        for (int j = 0; j < v.voxels()[1]; ++j) {
            for (int i = 0; i < v.voxels()[0]; ++i) {
                w_sums[idx(i, j, k)] += w_sums[idx(i, j, k - 1)];
            }
        }
    }

    return w_sums;
}

template <dimension D, typename T>
image<D, T> voxel_weights(const geometry::base<D, T>& geometry,
                          const volume<D, T>& volume) {
    auto w = image<D, T>(volume, (T)0.0);
    auto proj = dim::closest<D, T>(volume);
    for (auto [idx, line] : geometry) {
        (void)idx;
        for (auto melm : proj(line)) {
            w[melm.index] += 1.0;
        }
    }

    return partial_sums(w);
}

template <dimension D, typename T>
T weight(math::vec3<int> base, math::vec3<int> end, const image<D, T>& ws) {
    auto x = end.x;
    auto y = end.y;
    auto z = end.z;
    auto x0 = base.x;
    auto y0 = base.y;
    auto z0 = base.z;
    if (!(x >= x0 && y >= y0 && z >= z0)) {
        std::cout << math::to_string<3_D, int>(base) << " .. "
                  << math::to_string<3_D, int>(end) << "\n";
        assert(false);
    }

    auto sum = [&](auto i, auto j, auto k) -> T {
        if (i < 0 || j < 0 || k < 0) {
            return (T)0.0;
        }
        return ws[ws.index({i, j, k})];
    };

    auto a = sum(x, y, z);
    auto b = sum(x0 - 1, y, z);
    auto c = sum(x, y0 - 1, z);
    auto d = sum(x0 - 1, y0 - 1, z);
    auto e = sum(x, y, z0 - 1);
    auto f = sum(x0 - 1, y, z0 - 1);
    auto g = sum(x, y0 - 1, z0 - 1);
    auto h = sum(x0 - 1, y0 - 1, z0 - 1);

    auto v1 = a - b - c + d;
    auto v2 = e - f - g + h;

    return v1 - v2;
}

template <dimension D, typename T>
auto bounds_weight(box_t<D> bound, const image<D, T>& ws) {
    auto base = math::vec3<int>{bound[0][0], bound[1][0], bound[2][0]};
    auto end = math::vec3<int>{bound[0][1], bound[1][1], bound[2][1]};
    return weight(base, end, ws);
}

template <dimension D, typename T>
std::pair<split_t, int> find_split(const std::vector<math::line<D, T>>& lines,
                                   std::vector<math::line<D, T>>& lines_left,
                                   std::vector<math::line<D, T>>& lines_right,
                                   box_t<D> bounds, const image<D, T>& ws,
                                   T max_epsilon = 0.2) {
    std::cout << "Finding split with epsilon: " << max_epsilon << "\n"
              << "for bounds: "
              << "[" << bounds[0][0] << ", " << bounds[0][1] << "], "
              << "[" << bounds[1][0] << ", " << bounds[1][1] << "], "
              << "[" << bounds[2][0] << ", " << bounds[2][1] << "] \n";

    struct crossing_event {
        math::vec<D, T> point;
        std::size_t line_index;
        math::vec<D, int> direction;
    };

    auto avg_weight = (T)0.5 * bounds_weight<D, T>(bounds, ws);

    std::vector<crossing_event> crossings;

    // we want to tag each line with an entry point and exit point
    // the priority queue depends on the dimension!
    // first intersect each line with the bounds
    size_t idx = 0;
    for (auto line : lines) {
        auto intersections = math::intersect_bounds<D, T>(line, bounds);
        if (intersections) {
            auto a = intersections.value().first;
            auto b = intersections.value().second;
            crossings.push_back({a, idx, math::sign<D, T>(b - a)});
            crossings.push_back({b, idx, math::sign<D, T>(a - b)});
        } else {
            std::cout << math::line_to_string(line) << " does not intersect\n"
                      << "for bounds: "
                      << "[" << bounds[0][0] << ", " << bounds[0][1] << "], "
                      << "[" << bounds[1][0] << ", " << bounds[1][1] << "], "
                      << "[" << bounds[2][0] << ", " << bounds[2][1] << "]\n";
        }
        ++idx;
    }

    T best_imbalance = max_epsilon;
    int best_overlap = std::numeric_limits<int>::max();
    int best_split = 0;
    int best_d = -1;

    for (int d = 0; d < D; ++d) {
        auto compare = [d](crossing_event lhs, crossing_event rhs) {
            return lhs.point[d] < rhs.point[d];
        };

        // everything is in voxel coordinates
        auto imbalance = [&](int i) -> T {
            auto base =
                math::vec3<int>{bounds[0][0], bounds[1][0], bounds[2][0]};
            auto end =
                math::vec3<int>{bounds[0][1], bounds[1][1], bounds[2][1]};
            auto split_before = end;
            auto split_after = base;
            split_before[d] = i;
            split_after[d] = i + 1;

            auto weight_left = weight(base, split_before, ws);
            auto weight_right = weight(split_after, end, ws);

            return math::max(weight_left / (T)avg_weight,
                             weight_right / (T)avg_weight) -
                   (T)1.0;
        };

        std::sort(crossings.begin(), crossings.end(), compare);

        int overlap = 0;

        // only if point changes
        int last_split = bounds[d][0];
        for (auto current_crossing = 0u; current_crossing < crossings.size(); ++current_crossing) {
            auto crossing = crossings[current_crossing];

            int split = (int)crossings[current_crossing].point[d];
            if (split != last_split && split > bounds[d][0]) {
                int half_split = ((split - 1) + last_split) / 2;
                // FIXME average between the two again
                last_split = split;
                if (half_split >= bounds[d][1]) {
                    break;
                }

                // FIXME we dont get here *after* the final crossing
                auto epsilon = imbalance(half_split);
                if ((overlap < best_overlap && epsilon < max_epsilon) ||
                    (overlap == best_overlap && epsilon < best_imbalance)) {
                    best_overlap = overlap;
                    best_imbalance = epsilon;
                    best_split = half_split;
                    best_d = d;
                }
            }

            overlap += crossing.direction[d];
        }
    }

    if (best_d < 0) {
        // FIXME output debug info
        // using box_t = std::array<math::vec2<int>, D>;
        std::cout << "Can't find split for bounds: "
                  << "[" << bounds[0][0] << ", " << bounds[0][1] << "], "
                  << "[" << bounds[1][0] << ", " << bounds[1][1] << "], "
                  << "[" << bounds[2][0] << ", " << bounds[2][1] << "] "
                  << "crossings.size(): " << crossings.size() << ", "
                  << "lines.size(): " << lines.size() << ", "
                  << "\n";
        best_d = 0;
        best_split = (T)0.5 * (bounds[0][0] + bounds[0][1]);
    }
    auto best_compare = [best_d](crossing_event lhs, crossing_event rhs) {
        return lhs.point[best_d] < rhs.point[best_d];
    };
    std::sort(crossings.begin(), crossings.end(), best_compare);

    std::set<std::size_t> indices_left;
    std::set<std::size_t> indices_right;

    for (auto crossing : crossings) {
        if ((int)crossing.point[best_d] <= best_split) {
            // event is happening left
            if (crossing.direction[best_d] > 0) {
                indices_left.insert(crossing.line_index);
                indices_right.insert(crossing.line_index);
            } else if (crossing.direction[best_d] < 0) {
                assert(indices_right.find(crossing.line_index) != indices_right.end());
                indices_right.erase(crossing.line_index);
            } else {
                indices_left.insert(crossing.line_index);
            }
        } else {
            // event is happening right
            if (crossing.direction[best_d] >= 0) {
                indices_right.insert(crossing.line_index);
            }
        }
    }

    for (auto left_idx : indices_left) {
        lines_left.push_back(lines[left_idx]);
    }

    for (auto right_idx : indices_right) {
        lines_right.push_back(lines[right_idx]);
    }

    std::cout << "Split found " << best_d << " " << best_split << "\n";
    return {split_t{best_d, (int)(best_split)}, best_overlap};
}

template <dimension D, typename T>
bulk::util::binary_tree<bulk::util::split>
partition_bisection(const tomo::geometry::base<D, T>& geometry,
                    tomo::volume<D, T> object_volume, int processors,
                    T max_epsilon = 0.2) {
    bulk::util::binary_tree<bulk::util::split> result;

    // We want the resulting partitioning to be portable for different detector
    // configurations. We do this in two steps:
    // - Do the splitting in physical coordinates here
    // - Make a conversion function in voxel-based coordinates

    // assert that the number of processors is a power of two
    // TODO: support for non-2^k-way partitionings
    assert(math::is_power_of_two(processors));

    // compute the 'depth' of the tree (i.e. log(p))
    auto depth = (int)log2(processors);
    assert(1 << depth == processors);

    box_t<D> bounds = {};
    for (int d = 0; d < D; ++d) {
        bounds[d][0] = 0;
        bounds[d][1] = object_volume.voxels()[d] - 1;
    }

    // container for all the lines intersecting the volume
    std::vector<math::line<D, T>> all_lines;
    for (auto [idx, l] : geometry) {
        (void)idx;
        auto src = math::to_voxel<D, T>(l.source, object_volume);
        auto det = math::to_voxel<D, T>(l.detector, object_volume);
        auto line = math::line<D, T>{src, math::normalize(det - src)};
        if (math::intersect_bounds<D, T>(line, bounds)) {
            all_lines.push_back(line);
        }
        // ... else it misses the volume entirely, so we ignore it
    }

    auto ws = voxel_weights<D, T>(geometry, object_volume);

    using node = typename bulk::util::binary_tree<split_t>::node;
    using dir = typename bulk::util::binary_tree<split_t>::dir;

    // all the information we need to split a subvolume and save it
    struct subvolume {
        box_t<D> bounds;
        node* parent;
        dir direction;
        std::vector<math::line<D, T>> lines;
        int depth;
        int p;
        T epsilon;
    };

    // state should be the correct lines and bounds
    std::stack<subvolume> split_stack;

    auto complete_volume = subvolume{bounds, nullptr,    dir::left,  all_lines,
                                     0,      processors, max_epsilon};
    split_stack.push(complete_volume);

    int total_overlap = 0;

    while (!split_stack.empty()) {
        const auto& sub = split_stack.top();

        int sub_depth = sub.depth;
        auto sub_bounds = sub.bounds;
        auto sub_p = sub.p;
        auto sub_epsilon = sub.epsilon;

        auto q = (int)log2(sub.p);

        // containers to hold the left and right lines
        std::vector<math::line<D, T>> left;
        std::vector<math::line<D, T>> right;

        // we now have the subvolume, and want to find the best split
	split_t best_split;
	int best_overlap;
        std::tie(best_split, best_overlap)  = find_split<D, T>(
            sub.lines, left, right, sub_bounds, ws, sub.epsilon / q);
        total_overlap += best_overlap;

        // add the split to the tree
        node* current_node = result.add(sub.parent, sub.direction, best_split);

        split_stack.pop();

        // now right and left contain lines for recursion
        // we move these into the 'lower splits'
        if (sub_depth + 1 < depth) {
            auto max_weight = ((T)1 + sub_epsilon) *
                              (bounds_weight<D, T>(sub_bounds, ws) / sub_p);

            auto bounds_left = sub_bounds;
            bounds_left[best_split.d][1] = best_split.a;
            auto eps_left =
                max_weight / bounds_weight<D, T>(bounds_left, ws) * sub_p / 2 -
                1;
            split_stack.push(subvolume{bounds_left, current_node, dir::left,
                                       std::move(left), sub_depth + 1,
                                       sub_p / 2, eps_left});

            auto bounds_right = sub_bounds;
            bounds_right[best_split.d][0] = best_split.a + 1;
            auto eps_right =
                max_weight / bounds_weight<D, T>(bounds_right, ws) * sub_p / 2 -
                1;
            split_stack.push(subvolume{bounds_right, current_node, dir::right,
                                       std::move(right), sub_depth + 1,
                                       sub_p / 2, eps_right});
        }
    }

    std::cout << "Partitioning succeeded with volume: " << total_overlap << "\n";

    return result;
}

} // namespace distributed
} // namespace tomo
