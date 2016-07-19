#pragma once

#include <iostream>
#include <vector>

#include "math.hpp"
#include "geometry.hpp"

namespace tomo {

template <dimension Dimension, typename Type>
std::ostream& operator<<(std::ostream& os, const line<Dimension, Type>& line) {
    os << "{{" << line.origin[0] << ',' << line.origin[1] << "}, {"
       << line.delta[0] << ',' << line.delta[1] << "}}";
    return os;
}


template <typename T>
std::ostream& operator<<(std::ostream& os, const math::matrix_element<double>& elem) {
    os << "[" << elem.index << ": " << elem.value << "]";
    return os;
}


template <typename T>
std::ostream& operator<<(std::ostream& os, const math::vec2<double>& vec) {
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

#define TOMO_LOG_VAR(VAR) (std::cout << #VAR << " = " << VAR << std::endl)

} // namespace tomo
