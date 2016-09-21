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
            if (!math::approx_equal(line.origin.x, -10.0f) ||
                !math::approx_equal(line.origin.y, 5.0f) ||
                !math::approx_equal(line.delta.x, 1.0f) ||
                !math::approx_equal(line.delta.y, 0.0f)) {
                static_lines_correct = false;
                break;
            }
        }
        CHECK(static_lines_correct);
    }

    SECTION("Detector arrays computed correctly") {
        int k = 10;
        auto v = tomo::volume<2_D>(k);
        auto g = static_trajectory_test(v, 1, (T)1, math::vec<1_D, int>{2}, 2);
        CHECK(g.lines() == 2);

        CHECK(g.get_line(0).origin == g.get_line(1).origin);

        auto detector_diff =
            (g.get_line(0).origin + 20.0f * g.get_line(0).delta) -
            (g.get_line(1).origin + 20.0f * g.get_line(1).delta);

        CHECK(math::approx_equal(math::norm<2_D, float>(detector_diff), 1.0f));
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
            if (!math::approx_equal(line.origin.x, -10.0f) ||
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

    SECTION("Detector arrays computed correctly (3D)") {
        int k = 10;
        auto v = tomo::volume<3_D>(k);
        auto g =
            static_trajectory_test_3d(v, 1, (T)1, math::vec<2_D, int>{2}, 4);
        CHECK(g.lines() == 4);

        CHECK(g.get_line(0).origin == g.get_line(1).origin);

        auto detector_diff =
            (g.get_line(0).origin + 20.0f * g.get_line(0).delta) -
            (g.get_line(1).origin + 20.0f * g.get_line(1).delta);

        CHECK(math::approx_equal(math::norm<3_D, float>(detector_diff), 1.0f));

        auto detector_diag_diff =
            (g.get_line(0).origin + 20.0f * g.get_line(0).delta) -
            (g.get_line(3).origin + 20.0f * g.get_line(3).delta);

        CHECK(math::approx_equal(math::norm<3_D, float>(detector_diag_diff),
                                 math::sqrt2<float>));
    }

    SECTION("Detector tilts computed correctly (3D)") {
        // ... TODO
    }

    SECTION("Dual-axis PB") {
        int k = 10;
        auto v = tomo::volume<3_D>(k);
        auto g = geometry::dual_axis_parallel<T>(v, 20, (T)1,
                                                 math::vec<2_D, int>{2}, 4);
        CHECK(g.lines() == 80);
    }
}
