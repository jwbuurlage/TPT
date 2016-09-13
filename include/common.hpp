#pragma once

namespace tomo {

/** A compile-time type used to hold the dimension of a problem */
using dimension = int;

/** The default scalar type to use. */
using default_scalar_type = double;

} // namespace tomo

// FIXME: maybe only expose this in namespace, i.e. when using namespace tomo
/** A user defined literal for dimensions. */
constexpr tomo::dimension operator"" _D(unsigned long long d) { return d; }
