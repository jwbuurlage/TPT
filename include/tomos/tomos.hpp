/*
 * MIT License
 *
 * Copyright (c) 2016 Jan-Willem Buurlage
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "common.hpp"
#include "geometry.hpp"
#include "image.hpp"
#include "logging.hpp"
#include "math.hpp"
#include "operations.hpp"
#include "phantoms.hpp"
#include "projector.hpp"
#include "projections.hpp"
#include "utilities.hpp"
#include "volume.hpp"

#include "util/bench.hpp"
#include "util/report.hpp"
#include "util/tomo_args.hpp"
#include "util/read_metadata.hpp"
#include "util/read_tiff.hpp"

#include "algorithms/art.hpp"
#include "algorithms/sart.hpp"
#include "algorithms/sirt.hpp"

#include "distributed/forward_projection.hpp"
#include "distributed/partitioned_image.hpp"
#include "distributed/partitioned_projection_stack.hpp"
#include "distributed/recursive_bisectioning.hpp"
#include "distributed/trivial_partitioning.hpp"
#include "distributed/partitioning_statistics.hpp"
#include "distributed/restricted_geometry.hpp"

#include "geometries/cone.hpp"
//#include "geometries/dual_axis_parallel.hpp"
#include "geometries/dynamic_cone_beam.hpp"
#include "geometries/fan.hpp"
#include "geometries/helical_cone_beam.hpp"
#include "geometries/laminography.hpp"
#include "geometries/list.hpp"
#include "geometries/parallel.hpp"
#include "geometries/tomosynthesis.hpp"
#include "geometries/trajectory.hpp"

#include "projectors/closest.hpp"
#include "projectors/joseph.hpp"
#include "projectors/linear.hpp"

/** The overarching namespace for the tomos library. */
namespace tomo {

/** Reconstruction algorithms and methods. */
namespace reconstruction {}

/** Definition of common geometries used in experimental tomography. */
namespace geometry {}

/** Implementation of different discrete integration methods (DIM). */
namespace dim {}

/** Support for distributed tomography. */
namespace distributed {}

} // namespace tomo

#ifndef TOMO_NO_LITERALS
// This exposes user-defined literals,
// TODO do we want to disallow this by default?
using namespace tomo::literals;
#endif
