#pragma once

#include <array>

#include <glm/glm.hpp>

#include "../common.hpp"

namespace tpt {

using namespace tpt::literals;

namespace math {

/** The type to use for D-dimensional vectors. */
template <dimension D, typename T>
struct vec_type {
    // FIXME: for higher dimensional support, add own vector class for non-glm
    // sized vectors
    /* using type = my_own_vector_class<D, T>; */
};


/** The type to use for D-dimensional square matrices. */
template <dimension D, typename T>
struct mat_type {
    // FIXME: for higher dimensional support, add own mat class for non-glm
    // sized mats
    /* using type = my_own_mat_class<D, T>; */
};

/** The type to use for 1-dimensional vectors. */
template <typename T>
struct vec_type<1_D, T> {
    using type = glm::tvec1<T>;
};

/** The type to use for 2-dimensional vectors. */
template <typename T>
struct vec_type<2_D, T> {
    using type = glm::tvec2<T>;
};

/** The type to use for 3-dimensional vectors. */
template <typename T>
struct vec_type<3_D, T> {
    using type = glm::tvec3<T>;
};

/** The type to use for 4-dimensional vectors. */
template <typename T>
struct vec_type<4_D, T> {
    using type = glm::tvec4<T>;
};

/** The type to use for 4x4 matrices. */
template <typename T>
struct mat_type<4_D, T> {
    using type = glm::tmat4x4<T>;
};

/**
 * Generic vector type for any dimension, specializing where necessary.
 * \tparam D the dimension of the vector
 * \tparam T the scalar type to use
 */
template <dimension D, typename T>
using vec = typename vec_type<D, T>::type;

/**
 * Generic square matrix type for any dimension, specializing where necessary.
 * \tparam D the dimension of the matrix
 * \tparam T the scalar type to use
 */
template <dimension D, typename T>
using mat = typename mat_type<D, T>::type;

/** Short-hand for 1-dimensional vectors. */
template <typename T>
using vec1 = vec<1_D, T>;

/** Short-hand for 2-dimensional vectors. */
template <typename T>
using vec2 = vec<2_D, T>;

/** Short-hand for 3-dimensional vectors. */
template <typename T>
using vec3 = vec<3_D, T>;

/** Short-hand for 4-dimensional vectors. */
template <typename T>
using vec4 = vec<4_D, T>;

} // namespace math
} // namespace tpt
