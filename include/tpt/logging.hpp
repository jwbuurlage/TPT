#pragma once

#include <iostream>
#include <vector>

#include "fmt/format.h"
#include "fmt/ostream.h"

#include "geometry.hpp"
#include "math.hpp"

namespace tpt {

template <dimension Dimension, typename Type>
std::ostream& operator<<(std::ostream& os,
                         const math::line<Dimension, Type>& line) {
    os << "{{" << line.origin[0] << ',' << line.origin[1] << "}, {"
       << line.delta[0] << ',' << line.delta[1] << "}}";
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os,
                         const math::matrix_element<T>& elem) {
    os << "[" << elem.index << ": " << elem.value << "]";
    return os;
}

template <dimension D, typename T>
std::ostream& operator<<(std::ostream& os, const math::vec<D, T>& vec) {
    os << "{" << vec[0] << ", " << vec[1] << "}";
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    os << "[";
    auto sep = "";
    for (auto& elem : vec) {
        os << sep << elem;
        sep = ", ";
    }
    os << "]";
    return os;
}

#define TOMO_LOG_VAR(VAR) (fmt::print(#VAR " = {}\n", VAR))

} // namespace tpt
