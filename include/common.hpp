#pragma once

namespace tomo {

using dimension = int;
using default_scalar_type = double;

} // namespace tomo

// FIXME: maybe only expose this in namespace, i.e. when using namespace tomo
constexpr int operator"" _D(unsigned long long d) { return d; }
