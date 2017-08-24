/* This ile contains some utilities for working with partitioning trees. */

#include <fstream>
#include <memory>
#include <sstream>
#include <string>

#include <bulk/bulk.hpp>

namespace tomo {

template <typename T>
struct tree_node {
    // axis along which we split
    int d;

    // point along axis d at which the *second* subvolume begins
    // normalized to [0,1]
    T a;

    tree_node() : d(-1), a(0) {}
    tree_node(int d_, T a_) : d(d_), a(a_) {}

    tree_node(tree_node&& other)
        : left(std::move(other.left)), right(std::move(other.right)) {
        d = other.d;
        a = other.a;
    }

    std::unique_ptr<tree_node> left;
    std::unique_ptr<tree_node> right;
};

template <typename T>
void add_to_tree(tree_node<T>* neutral,
                 bulk::util::binary_tree<bulk::util::split>::node* original,
                 tomo::volume<3_D, T> v) {
    if (neutral && original) {
        neutral->d = original->value.d;
        neutral->a = (T)(original->value.a + 1) / v.voxels()[neutral->d]; // we save the mathematical split
        if (original->left.get()) {
            neutral->left = std::make_unique<tree_node<T>>();
            add_to_tree(neutral->left.get(), original->left.get(), v);
        }
        if (original->right.get()) {
            neutral->right = std::make_unique<tree_node<T>>();
            add_to_tree(neutral->right.get(), original->right.get(), v);
        }
    }
}

template <typename T>
tree_node<T> to_neutral_tree(bulk::util::binary_tree<bulk::util::split>& splits,
                             tomo::volume<3_D, T> v) {
    tree_node<T> neutral_root;
    add_to_tree(&neutral_root, splits.root.get(), v);
    return neutral_root;
}

template <typename T>
void add_to_voxel_tree(bulk::util::binary_tree<bulk::util::split>::node* voxel,
                       tree_node<T>* neutral, tomo::volume<3_D, T> v) {
    if (neutral && voxel) {
        voxel->value.d = neutral->d;
        voxel->value.a = (int)(neutral->a * v.voxels()[neutral->d] + 0.5) - 1;
        if (neutral->left.get()) {
            voxel->left = std::make_unique<
                bulk::util::binary_tree<bulk::util::split>::node>(
                bulk::util::split{0, 0});
            add_to_voxel_tree<T>(voxel->left.get(), neutral->left.get(), v);
        }
        if (neutral->right.get()) {
            voxel->right = std::make_unique<
                bulk::util::binary_tree<bulk::util::split>::node>(
                bulk::util::split{0, 0});
            add_to_voxel_tree<T>(voxel->right.get(), neutral->right.get(), v);
        }
    }
}

template <typename T>
bulk::util::binary_tree<bulk::util::split>
from_neutral_tree(tree_node<T>& splits, tomo::volume<3_D, T> v) {
    bulk::util::binary_tree<bulk::util::split> result;
    result.root =
        std::make_unique<bulk::util::binary_tree<bulk::util::split>::node>(
            bulk::util::split{0, 0});
    add_to_voxel_tree(result.root.get(), &splits, v);
    return result;
}

template <typename T>
void add_to_output(tree_node<T>* tree, std::stringstream& result) {
    if (tree == nullptr) {
        result << "[]";
        return;
    }

    result << "[[" << tree->d << ", " << tree->a << "], ";
    add_to_output(tree->left.get(), result);
    result << ", ";
    add_to_output(tree->right.get(), result);
    result << "]";
}

template <typename T>
void serialize_tree(tree_node<T>& tree, std::string filename) {
    std::stringstream result;
    add_to_output(&tree, result);

    std::ofstream of(filename, std::ios::out);
    of << result.str();

    std::cout << "Saved file: " << filename << "\n";
}

template <typename T>
std::unique_ptr<tree_node<T>> deserialize_node(std::string input, int depth) {
    if (depth <= 0) {
        return nullptr;
    }
    if (input.size() < 2 || input[0] != '[') {
        return nullptr;
    }
    if (input[1] == ']') {
        return nullptr;
    }

    std::vector<int> start_positions;
    std::vector<int> delim_positions;

    int imbalance = 0;
    for (int i = 1; i < (int)input.size(); ++i) {
        if (input[i] == '[') {
            imbalance++;
            if (imbalance == 1) {
                start_positions.push_back(i);
            }
        } else if (input[i] == ']') {
            imbalance--;
            if (imbalance == 0) {
                delim_positions.push_back(i);
            }
        }
    }

    assert(start_positions.size() == 3);
    assert(start_positions.size() == delim_positions.size());

    auto get_part = [&](int i) -> std::string {
        return input.substr(start_positions[i],
                            delim_positions[i] - start_positions[i] + 1);
    };

    auto split_value = [&](std::string value_string) -> std::pair<int, T> {
        std::pair<int, T> value;
        auto comma = value_string.find(", ");
        std::stringstream d_stream(value_string.substr(1, comma + 1));
        d_stream >> value.first;
        std::stringstream a_stream(value_string.substr(comma + 2));
        a_stream >> value.second;
        return value;
    };

    auto spl = split_value(get_part(0));
    auto node = std::make_unique<tree_node<T>>(spl.first, spl.second);
    node->left = std::move(deserialize_node<T>(get_part(1), depth - 1));
    node->right = std::move(deserialize_node<T>(get_part(2), depth - 1));

    return node;
}

template <typename T>
std::unique_ptr<tree_node<T>> deserialize_tree(std::string filename,
                                               int depth) {
    std::ifstream input(filename, std::ios::in);
    std::string tree_line;
    std::getline(input, tree_line);
    return deserialize_node<T>(tree_line, depth);
}

void print_bulk_node_(bulk::util::binary_tree<bulk::util::split>::node* node,
                      int depth = 0) {
    if (!node)
        return;

    for (int i = 0; i < depth; ++i) {
        std::cout << "> ";
    }
    std::cout << "[" << node->value.d << ", " << node->value.a << "]\n";
    if (auto left = node->left.get()) {
        print_bulk_node_(left, depth + 1);
    }
    if (auto right = node->right.get()) {
        print_bulk_node_(right, depth + 1);
    }
}

void print_tree(bulk::util::binary_tree<bulk::util::split>& tree) {
    print_bulk_node_(tree.root.get());
}

template <typename T>
void print_neutral_node_(tree_node<T>* node, int depth = 0) {
    for (int i = 0; i < depth; ++i) {
        std::cout << "> ";
    }
    std::cout << "[" << node->d << ", " << node->a << "]\n";
    if (auto left = node->left.get()) {
        print_neutral_node_(left, depth + 1);
    }
    if (auto right = node->right.get()) {
        print_neutral_node_(right, depth + 1);
    }
}

template <typename T>
void print_neutral_tree(tree_node<T>& node) {
    print_neutral_node_(&node);
}

template <typename T>
std::unique_ptr<bulk::tree_partitioning<3_D>>
load_partitioning(std::string filename, tomo::volume<3_D, T> v,
                  int depth = -1) {
    auto reloaded = tomo::deserialize_tree<T>(filename, depth);
    if (!reloaded) {
        std::cout << "Could not deserialize: " << filename << "\n";
        return nullptr;
    }

    auto large_voxel = tomo::from_neutral_tree<T>(*reloaded, v);

    return std::make_unique<bulk::tree_partitioning<3_D>>(
        math::vec_to_array<3_D, int>(v.voxels()), large_voxel.size() + 1,
        std::move(large_voxel));
}

} // namespace tomo
