#pragma once

#include "partitioning.hpp"
#include <boost/geometry.hpp>

namespace tpt::grcb::plot {

/**
   From a collection of:
   xs_det, ys_det: projected points
   hull_xs, hull_ys: convex hulls of xs, ys
   det: polygon of the detector,
   name: filename to save

   writes an SVG {name}.svg
 */
template <typename T>
void overlap(std::array<math::vec2<T>, 8> xs_det,
             std::array<math::vec2<T>, 8> ys_det,
             bg::model::polygon<math::vec2<T>> hull_xs,
             bg::model::polygon<math::vec2<T>> hull_ys,
             bg::model::polygon<math::vec2<T>> det, std::string name,
             bool svg = true, bool wkt = false) {
    if (svg) {
        std::ofstream svg(name + ".svg");

        auto mapper = bg::svg_mapper<math::vec2<T>, false>(svg, 500, 500);
        mapper.add(det);
        mapper.add(hull_xs);
        mapper.add(hull_ys);

        mapper.map(det,
                   "fill-opacity:1.0;fill:rgb(255,255,255);stroke:rgb(0,0,0);"
                   "stroke-width:2");
        for (auto x : xs_det) {
            mapper.map(x, "fill-opacity:0.5;fill:rgb(0,0,255);", 5);
        }

        for (auto y : ys_det) {
            mapper.map(y, "fill-opacity:0.5;fill:rgb(255,0,0);", 5);
        }

        mapper.map(hull_xs, "fill-opacity:0.2;fill:rgb(0,0,255);");
        mapper.map(hull_ys, "fill-opacity:0.2;fill:rgb(255,0,0);");
    }

    if (wkt) {
        std::ofstream wkt(name + ".wkt");
        for (auto x : xs_det) {
            wkt << bg::wkt(x) << "\n";
        }
        for (auto y : ys_det) {
            wkt << bg::wkt(y) << "\n";
        }
        wkt << bg::wkt(hull_xs) << "\n";
        wkt << bg::wkt(hull_ys) << "\n";
        wkt << bg::wkt(det) << "\n";
    }
}

namespace detail {
template <typename T>
void add_shadow(bg::svg_mapper<math::vec2<T>, false>& mapper,
                const std::unique_ptr<node<T>>& x, std::array<vec3<T>, 8> v,
                geometry::projection<3_D, T> pi, int& i, int p) {
    auto [a, b] = min_max_cube<T>(v);
    auto [vl, vr] = split_at(v, (T)0.5 * (a[x->a] + b[x->a]), x->a);
    if (!x->left) {
        auto [pl, shadowl] = shadow(pi, vl);
        mapper.map(shadowl, std::string("fill-opacity:0.5;fill:hsl(") +
                                std::to_string(i * (256 / p)) + ",100%,50%);");

        ++i;
        auto [pr, shadowr] = shadow(pi, vr);
        mapper.map(shadowr, std::string("fill-opacity:0.5;fill:hsl(") +
                                std::to_string(i * (256 / p)) + ",100%,50%);");
        ++i;

        (void)pl;
        (void)pr;
    } else {
        add_shadow(mapper, x->left, vl, pi, i, p);
        add_shadow(mapper, x->right, vr, pi, i, p);
    }
}
} // namespace detail

/**
   From a collection of:
   xs_det, ys_det: projected points
   hull_xs, hull_ys: convex hulls of xs, ys
   det: polygon of the detector,
   name: filename to save

   writes an SVG {name}.svg
 */
template <typename T>
void overlaps(geometry::projection<3_D, T> pi,
              const std::unique_ptr<node<T>>& part, std::array<vec3<T>, 8> v,
              std::string name, int p) {
    std::ofstream svg(name + ".svg");

    auto mapper = bg::svg_mapper<math::vec2<T>, false>(svg, 500, 500);
    auto det = detector(pi);
    mapper.add(det);
    mapper.map(det, "fill-opacity:1.0;fill:rgb(255,255,255);stroke:rgb(0,0,0);"
                    "stroke-width:2");

    auto i = 0;
    detail::add_shadow(mapper, part, v, pi, i, p);
}

template <typename T>
void shadow(geometry::projection<3_D, T> pi, std::array<vec3<T>, 8> v,
            std::string name) {
    std::ofstream svg(name + ".svg");

    auto mapper = bg::svg_mapper<math::vec2<T>, false>(svg, 500, 500);
    auto det = detector(pi);
    mapper.add(det);
    mapper.map(det, "fill-opacity:1.0;fill:rgb(255,255,255);stroke:rgb(0,0,0);"
                    "stroke-width:2");
    auto [xs, shadow] = grcb::shadow(pi, v);

    mapper.map(shadow, std::string("fill-opacity:0.5;fill:hsl(") +
                           std::to_string(0) + ",100%,50%);");

    for (auto x : xs) {
        mapper.map(x, "fill-opacity:1.0;fill:rgb(0,0,0);", 10);
    }
}

} // namespace tpt::grcb::plot
