// http://blog.notdot.net/2009/11/Damn-Cool-Algorithms-Spatial-indexing-with-Quadtrees-and-Hilbert-Curves

#include <cmath>

#include "../math.hpp"

namespace tpt {
namespace util {

class hilbert_curve {
  public:
    hilbert_curve(math::vec2<int> grid_size) : grid_size_(grid_size) {
        assert(grid_size[0] == grid_size[1]);
        assert(math::is_power_of_two(grid_size[0]));
        depth_ = std::log2<int>(grid_size[0]);
    }

    int to(math::vec2<int> idx) { return to_(idx, grid_size_, 0); }

    math::vec2<int> from(int idx) {
        int x = 0;
        int y = 0;
        from_(idx, grid_size_[0], 0, x, y, depth_);
        return {x, y};
    }

  private:
    int to_(math::vec2<int> idx, math::vec2<int> grid_size, int type) {
        (void)idx;
        (void)grid_size;
        (void)type;
        assert(0); // NOT IMPLEMENTED
        return 0;
    }

    void from_(int idx, int grid_size, int type, int& x, int& y, int depth) {
        if (depth == 0) {
            return;
        }

        // recurse
        auto lookup = table[type];

        // the local index at the current level
        int local = std::get<0>(lookup)[(idx >> 2 * (depth - 1)) & 3];

        // update x and y accordingly
        // FIXME (check me rather)
        y += (local / 2) * (grid_size / 2);
        x += (local % 2) * (grid_size / 2);

        // recurse
        from_(idx, grid_size / 2, std::get<1>(lookup)[local], x, y, depth - 1);
    }

    math::vec2<int> grid_size_;
    int depth_;

    /*
     * We index a 2x2 grid as follows:
     *  +-------+
     *  | 0 | 1 |
     *  +---|---+
     *  | 2 | 3 |
     *  +-------+
     *
     *  The following table defines the hilbert curve recursively. table[i]
     *  contains the information for when a point falls into quad i.
     *
     *  - The first entry of table[i] are 4 integers, representing the order in
     *  which the quads are followed
     *  - The second entry of table[i] are 4 integers, representing the
     * recursive quad type, for the index falls into the cell
     */
    using afi = std::array<int, 4>;
    using quad = std::tuple<afi, afi>;
    static constexpr auto table =
        std::array<quad, 4>{quad{afi{0, 2, 3, 1}, afi{3, 1, 0, 0}},
                            quad{afi{3, 2, 0, 1}, afi{1, 0, 1, 2}},
                            quad{afi{3, 1, 0, 2}, afi{2, 2, 3, 1}},
                            quad{afi{0, 1, 3, 2}, afi{0, 3, 2, 3}}};
};

} // namespace util
} // namespace tpt
