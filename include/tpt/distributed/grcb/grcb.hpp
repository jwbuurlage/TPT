#pragma once

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/register/multi_point.hpp>
#include <boost/geometry/geometries/register/point.hpp>

#include "datatypes.hpp"

namespace bg = boost::geometry;

// register data structures with boost_geometry
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
BOOST_GEOMETRY_REGISTER_POINT_2D(tpt::math::vec2<float>, float,
                                 bg::cs::cartesian, x, y);
BOOST_GEOMETRY_REGISTER_POINT_3D(tpt::math::vec3<float>, float,
                                 bg::cs::cartesian, x, y, z);
BOOST_GEOMETRY_REGISTER_MULTI_POINT(std::vector<tpt::math::vec2<float>>);

template <typename T>
using mvec2 = std::array<tpt::math::vec2<T>, 8>;

BOOST_GEOMETRY_REGISTER_MULTI_POINT(mvec2<float>);
#pragma GCC diagnostic pop

#include "cube.hpp"
#include "geometry.hpp"
#include "integration.hpp"
#include "load_balance.hpp"
#include "partitioning.hpp"
#include "plot.hpp"
