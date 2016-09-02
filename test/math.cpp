#include "catch.hpp"
#include "tomo.hpp"

using T = float;

TEST_CASE("Basic operations on vectors", "[math]") {}

TEST_CASE("Intersection and box checking", "[math]") {
    SECTION("2D") {
        int k = 4;
        tomo::volume<2_D> v(k);
        std::vector<tomo::math::vec2<T>> xs = {
            {0, 0}, {k, k}, {k / 2, k / 2}, {k / 3, k / 2},
            {0, k}, {k, 0}, {1, 1}};
        std::vector<bool> is_inside = {false, false, true, true,
                                    false, false, true};
        assert(is_inside.size() == xs.size());
        for (std::size_t i = 0; i < is_inside.size(); ++i) {
            bool inside = tomo::math::inside<2_D, T>(xs[i], v);
            CHECK(inside == is_inside[i]);
        }
    }

    SECTION("3D") {
        int k = 4;
        tomo::volume<3_D> v(k);
        std::vector<tomo::math::vec3<T>> xs = {
            {0, 0, 0}, {k, k, 0}, {k / 2, k / 2, k / 2}, {k / 3, k / 2, k / 4},
            {0, k, 0}, {k, 0, 0}, {1, 1, 1}};
        std::vector<bool> is_inside = {false, false, true, true,
                                    false, false, true};
        assert(is_inside.size() == xs.size());
        for (std::size_t i = 0; i < is_inside.size(); ++i) {
            bool inside = tomo::math::inside<3_D, T>(xs[i], v);
            CHECK(inside == is_inside[i]);
        }
    }
}

TEST_CASE("Interpolation", "[math]") {
    SECTION("2D") {}

    SECTION("3D") {}
}
