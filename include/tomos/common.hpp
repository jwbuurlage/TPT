#pragma once

#include <memory>
#include <vector>

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

/* Meta programming hacks */
template <bool...>
struct bool_pack {};

template <bool... bs>
using all_true = std::is_same<bool_pack<bs..., true>, bool_pack<true, bs...>>;

template <class R, class... Ts>
using are_all_convertible = all_true<std::is_convertible<Ts, R>::value...>;

template <int count, class R, class... Ts>
using count_of_type =
    all_true<sizeof...(Ts) == count, are_all_convertible<R, Ts...>::value>;

template <int D, typename... Ts>
using check_dim =
    typename std::enable_if<count_of_type<D, int, Ts...>::value>::type;

} // namespace tomo
