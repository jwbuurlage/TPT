#pragma once

#include "bulk/bulk.hpp"

namespace tomo {
namespace distributed {

template <dimension D, typename T>
class partitioned_projection_stack {
  public:
    struct exchange {
        int line;
        int target;
    };

    using value_type = T;

    partitioned_projection_stack(bulk::world& world, bulk::partitioning<D>& part,
                         tomo::geometry::base<D, T>& geometry)
        : world_(world), data_(world, geometry.lines()), part_(part),
          geometry_(geometry) {
        clear();
    }

    void clear() { std::fill(data_.begin(), data_.end(), 0); }

    int communication_volume() { return exchanges_.size(); }

    void harmonize() {
        auto sino_q = bulk::queue<int, T>(world_);
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

    // TODO:
    // We may not want to redundantly send everything to everyone,
    // but rather sum on one of the owners, which then broadcasts
    // to a subset
    template <typename Projector>
    void compute_overlap(Projector& proj) {
        auto cyclic = bulk::cyclic_partitioning<1, 1>(
            world_, {geometry_.lines()}, {world_.active_processors()});

        auto q = bulk::queue<int, int>(world_);
        int idx = 0;
        for (auto l : geometry_) {
            for (auto elem : proj(l)) {
                (void)elem;
                auto owner = cyclic.owner({idx});
                q(owner).send(idx, world_.processor_id());
                break;
            }
            ++idx;
        }

        world_.sync();

        // targets[a] = list of procs for line a
        std::vector<std::vector<int>> targets(
            cyclic.local_count(world_.processor_id()));

        for (auto& msg : q) {
            targets[msg.tag / world_.active_processors()].push_back(
                msg.content);
        }

        auto exchange_queue = bulk::queue<int, int>(world_);

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
    bulk::world& world_;
    bulk::coarray<T> data_;
    bulk::partitioning<D>& part_;
    tomo::geometry::base<D, T>& geometry_;

    std::vector<exchange> exchanges_;
};

} // namespace distributed
} // namespace tomo
