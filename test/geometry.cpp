#include "catch.hpp"
#include "tomo.hpp"

using T = float;

TEST_CASE("Trajectory based geometry", "[geometry]") {
    using namespace tomo;

    class static_trajectory_test : public geometry::trajectory<2_D, float> {
      public:
        using super = geometry::trajectory<2_D, float>;
        using super::super;

        math::vec<2_D, float> source_location(int) const override final {
            return {-10.0f, 5.0f};
        }

        math::vec<2_D, float> detector_location(int) const override final {
            return {10.0f, 5.0f};
        }

        std::array<math::vec<2_D, float>, 1>
        detector_tilt(int) const override final {
            return {math::vec<2_D, float>{0.0f, 1.0f}};
        }
    };

    SECTION("Basic trajectories") {
        int k = 10;
        auto v = tomo::volume<2_D>(k);
        auto g = static_trajectory_test(v, 10);
        CHECK(g.lines() == 10);
        bool static_lines_correct = true;
        for (auto line : g) {
            if (!math::approx_equal(line.origin.x, 0.0f) ||
                !math::approx_equal(line.origin.y, 5.0f) ||
                !math::approx_equal(line.delta.x, 1.0f) ||
                !math::approx_equal(line.delta.y, 0.0f)) {
                static_lines_correct = false;
                break;
            }
        }
        CHECK(static_lines_correct);
    }

    class static_trajectory_test_3d : public geometry::trajectory<3_D, float> {
      public:
        using super = geometry::trajectory<3_D, float>;
        using super::super;

        math::vec<3_D, float> source_location(int) const override final {
            return {-10.0f, 5.0f, 5.0f};
        }

        math::vec<3_D, float> detector_location(int) const override final {
            return {10.0f, 5.0f, 5.0f};
        }

        std::array<math::vec<3_D, float>, 2>
        detector_tilt(int) const override final {
            return {math::vec<3_D, float>{0.0f, 1.0f, 0.0f},
                    math::vec<3_D, float>{0.0f, 0.0f, 1.0f}};
        }
    };

    SECTION("Basic trajectories (3D)") {
        int k = 10;
        auto v = tomo::volume<3_D>(k);
        auto g = static_trajectory_test_3d(v, 10);
        CHECK(g.lines() == 10);

        bool static_lines_correct = true;
        for (auto line : g) {
            if (!math::approx_equal(line.origin.x, 0.0f) ||
                !math::approx_equal(line.origin.y, 5.0f) ||
                !math::approx_equal(line.origin.z, 5.0f) ||
                !math::approx_equal(line.delta.x, 1.0f) ||
                !math::approx_equal(line.delta.y, 0.0f) ||
                !math::approx_equal(line.delta.z, 0.0f)) {
                static_lines_correct = false;
                break;
            }
        }
        CHECK(static_lines_correct);
    }

    SECTION("Detector tilts computed correctly (3D)") {
        int k = 10;
        auto v = tomo::volume<3_D>(k);
        auto g = geometry::cone_beam<T>(v, 20);

        auto axes = g.detector_tilt(0);
        CHECK(math::approx_equal(math::norm<3_D, T>(axes[0]), 1.0f));
        CHECK(math::approx_equal(math::norm<3_D, T>(axes[1]), 1.0f));

        auto rotated_axes = g.detector_tilt(1);
        CHECK(math::approx_equal(math::norm<3_D, T>(rotated_axes[0]), 1.0f));
        CHECK(math::approx_equal(math::norm<3_D, T>(rotated_axes[1]), 1.0f));
    }

    SECTION("Dual-axis PB") {
        int k = 10;
        int steps_per_axis = 10;
        int detector_width = 2;
        auto v = tomo::volume<3_D>(k);
        auto g = geometry::dual_axis_parallel<T>(
            v, steps_per_axis, (T)1, math::vec<2_D, int>{detector_width});
        CHECK(g.lines() ==
              2 * steps_per_axis * detector_width * detector_width);

        bool lines_make_sense = true;
        for (auto line : g) {
            (void)line;
        }
        CHECK(lines_make_sense);

    }
}
