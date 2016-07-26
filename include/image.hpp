#pragma once

#include <algorithm>
#include <stddef.h>
#include <type_traits>
#include <vector>

#include "common.hpp"
#include "volume.hpp"

namespace tomo {

template <dimension Dimension, typename Type = double>
class image {
  public:
    image(volume<Dimension> v) : v_(v) {
        int cells = 1;
        for (int i = 0; i < Dimension; ++i)
            cells *= v[i];
        data_.resize(cells);
    }

    size_t index(std::array<int, Dimension> xs) const {
       return v_.index(xs);
    }

    Type& operator[](size_t index) { return data_[index]; }
    const Type& operator[](size_t index) const { return data_[index]; }

    Type& operator()(std::array<int, Dimension> xs) {
        return (*this)[index(xs)];
    }

    const Type& operator()(std::array<int, Dimension> xs) const {
        return (*this)[index(xs)];
    }

    int size(int i) const { return v_[i]; }

    std::vector<Type>& mutable_data() { return data_; };
    const std::vector<Type>& data() const { return data_; };

    volume<Dimension> get_volume() const { return v_; }

    std::array<int, Dimension> dimensions() const { return v_.dimensions(); }

    void clear() { std::fill(data_.begin(), data_.end(), 0); }

  private:
    volume<Dimension> v_;
    std::vector<Type> data_;
};

} // namespace tomo
