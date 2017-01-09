#pragma once

#include <functional>

#include "../util/plotter.hpp"
#include "forward_projection.hpp"

namespace tomo {
namespace distributed {

template <dimension D, typename T>
class ext_plotter {
  public:
    ext_plotter(std::string address, std::string name, bulk::world& world)
        : plotter_() {
        if (world.processor_id() != 0)
            return;

        plotter_.connect(address, name);
    }

    template <int G>
    void plot(partitioned_image<D, G, T>& f) {
        auto gathered = f.gather();

        if (f.world().processor_id() != 0)
            return;

        plotter_.plot(gathered);
    }

    template <typename Geometry>
    void plot(partitioned_sinogram<D, T, Geometry>& sino) {
        if (sino.world().processor_id() != 0)
            return;

        plotter_.plot(sino.as_image());
    }

  private:
    tomo::ext_plotter<D, T> plotter_;
};

} // namespace distributed
} // namespace tomo
