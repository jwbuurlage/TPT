#pragma once

namespace tomo {

/** A compile-time type used to hold the dimension of a problem */
using dimension = int;

/** The default scalar type to use. */
using default_scalar_type = double;

/** User defined literals for the library. */
namespace literals {
    /** A user defined literal for dimensions. */
    constexpr tomo::dimension operator"" _D(unsigned long long d) { return d; }
}

} // namespace tomo

