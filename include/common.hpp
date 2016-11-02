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

namespace core {

template <typename T>
struct binary_tree {
    binary_tree() = default;
    binary_tree(binary_tree&& other) : root(std::move(other.root)) {}

    enum class dir { left, right };

    struct node {
        node(T value_) : value(value_) {}

        std::unique_ptr<node> left = nullptr;
        std::unique_ptr<node> right = nullptr;
        T value;
    };

    node* add(node* parent, dir direction, T value) {
        if (parent != nullptr) {
            if (direction == dir::left) {
                parent->left = std::make_unique<node>(value);
                return parent->left.get();
            } else {
                parent->right = std::make_unique<node>(value);
                return parent->right.get();
            }
        }
        else {
            root = std::make_unique<node>(value);
            return root.get();
        }
    }

    std::unique_ptr<node> root;
};
}

} // namespace tomo
