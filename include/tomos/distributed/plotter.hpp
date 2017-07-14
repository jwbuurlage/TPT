#pragma once

#include <functional>

#include "../util/plotter.hpp"
#include "forward_projection.hpp"
#include "partitioned_projection_stack.hpp"

namespace tomo {
namespace distributed {

template <dimension D, typename T>
class ext_plotter {
  public:
    ext_plotter(std::string address, std::string name, bulk::world& world)
        : plotter_() {
        if (world.processor_id() != 0)
            return;

        plotter_.connect(address);
        plotter_.make_scene(name);
    }

    template <int G>
    void plot(partitioned_image<D, G, T>& f) {
        auto gathered = f.gather();

        if (f.world().processor_id() != 0)
            return;

        plotter_.plot(gathered);
    }

    template <typename Geometry>
    void plot(partitioned_projection_stack<D, T>& ps) {
        if (ps.world().processor_id() != 0)
            return;

        plotter_.plot(ps.as_image());
    }

  private:
    tomo::util::ext_plotter<D, T> plotter_;
};

} // namespace distributed
} // namespace tomo
