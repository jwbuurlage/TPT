#pragma once

namespace tomo {
namespace distributed {

class object {};

class geometry : public object {
  public:
    geometry(tomo::geometry& g, partitioned_volume part_vol)
        : part_vol_(part_vol) {
        // construct from partitioned volume
    }

  private:
    partitioned_volume part_vol_;
};

class image : public object {

  private:
    partitioned_volume part_vol_;
};

} // namespace distributed
} // namespace tomo
