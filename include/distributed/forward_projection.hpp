#pragma once

#include <string>
#include <vector>

#include "common.hpp"
#include "image.hpp"
#include "partitioned_image.hpp"
#include "projector.hpp"
#include "projectors/linear.hpp"
#include "sinogram.hpp"
#include "volume.hpp"

#include "bulk/bulk.hpp"

namespace tomo {
namespace distributed {

template <dimension D, typename T, typename Geometry, typename World>
class partitioned_sinogram {
  public:
    struct exchange {
        int line;
        int target;
    };

    using value_type = T;

    partitioned_sinogram(World& world, bulk::partitioning<D, 1>& part,
                         Geometry& geometry)
        : world_(world), data_(world, geometry.lines()), part_(part),
          geometry_(geometry) {
        clear();
    }

    void clear() { std::fill(data_.begin(), data_.end(), 0); }

    int communication_volume() { return exchanges_.size(); }

    void harmonize() {
        auto sino_q = bulk::create_queue<int, T>(world_);
        // make exchange queue
        // exchange all overlaps
        for (auto msg : exchanges_) {
            sino_q(msg.target).send(msg.line, data_[msg.line]);
        }

        world_.sync();

        for (auto msg : sino_q) {
            data_[msg.tag] += msg.content;
        }
    }

    template <typename Projector>
    void compute_overlap(Projector& proj) {
        auto cyclic = bulk::cyclic_partitioning<1, 1>(
            {world_.active_processors()}, {geometry_.lines()});

        auto q = bulk::create_queue<int, int>(world_);
        int idx = 0;
        for (auto l : geometry_) {
            for (auto elem : proj(l)) {
                (void)elem;
                auto owner = cyclic.owner({idx})[0];
                q(owner).send(idx, world_.processor_id());
                break;
            }
            ++idx;
        }

        world_.sync();

        // targets[a] = list of procs for line a
        std::vector<std::vector<int>> targets(
            cyclic.local_extent({world_.processor_id()})[0]);

        for (auto& msg : q) {
            targets[msg.tag / world_.active_processors()].push_back(
                msg.content);
        }

        auto exchange_queue = bulk::create_queue<int, int>(world_);

        idx = 0;
        for (auto& line : targets) {
            int line_idx =
                (idx++) * world_.active_processors() + world_.processor_id();
            int procs = line.size();

            for (int s = 0; s < procs; ++s) {
                for (int t = 0; t < procs; ++t) {
                    if (s == t) {
                        continue;
                    }
                    exchange_queue(line[s]).send(line_idx, line[t]);
                }
            }
        }

        world_.sync();

        for (auto msg : exchange_queue) {
            exchanges_.push_back({msg.tag, msg.content});
        }
    }

    // FIXME do we want this here
    auto as_image() {
        auto v = tomo::volume<D>(geometry_.groups());
        tomo::image<D, T> result(v);

        int i = 0;
        for (auto& x : data_) {
            result[i++] = x;
        }

        return result;
    }

    auto& operator[](int idx) { return data_[idx]; }
    const auto& operator[](int idx) const { return data_[idx]; }

    auto get_volume() const { return geometry_.get_volume(); }
    const std::vector<exchange>& exchanges() const { return exchanges_; }

    auto& world() { return world_; }

  private:
    World& world_;
    bulk::coarray<T, World> data_;
    bulk::partitioning<D, 1>& part_;
    Geometry& geometry_;

    std::vector<exchange> exchanges_;
};

/**
 * Perform a forward-projection of a given image.
 *
 * TODO: alternatively, we can separate the 'geometry communication' info, and
 * the partitioned sinogram itself, so that we can still return it from this
 * function.
 * */
template <dimension D, typename T, class Geometry, class World, class Projector>
void forward_project(tomo::distributed::partitioned_image<D, T, World>& f,
                     const Geometry& g, Projector& proj,
                     partitioned_sinogram<D, T, Geometry, World>& sino) {
    int line_number = 0;
    for (auto line : g) {
        for (auto elem : proj(line)) {
            sino[line_number] += f[elem.index] * elem.value;
        }
        ++line_number;
    }

    sino.harmonize();
}

/**
 * Perform a back-projection of a given image.
 *
 * TODO: alternatively, we can separate the 'geometry communication' info, and
 * the partitioned sinogram itself, so that we can still return it from this
 * function.
 * */
template <dimension D, typename T, class Geometry, class World, class Projector>
void back_project(tomo::distributed::partitioned_image<D, T, World>& f,
                  const Geometry& g, Projector& proj,
                  partitioned_sinogram<D, T, Geometry, World>& sino) {
    int line_number = 0;
    for (auto line : g) {
        for (auto elem : proj(line)) {
            f[elem.index] += sino[line_number] * elem.value;
        }
        ++line_number;
    }
}

} // namespace distributed
} // namespace tomo
