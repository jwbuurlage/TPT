#include <cmath>
#include <iostream>
#include <memory>

#include <glm/gtx/string_cast.hpp>

#include "tomos/tomos.hpp"
#include "tomos/util/trees.hpp"

using T = float;
constexpr tomo::dimension D = 3_D;

int main() {
    namespace td = tomo::distributed;

    int size = 32;
    int procs = 8;
    T epsilon = 0.1;

    auto v = tomo::volume<D, T>(size);
    auto g = tomo::geometry::cone_beam<T>(v, size, {(T)2.0, (T)2.0},
                                          {size, size}, (T)4.0, (T)4.0);

    // [x] partition
    auto tree = td::partition_bisection(g, v, procs, epsilon);
    tomo::print_tree(tree);

    // [x] convert to neutral
    auto neutral = tomo::to_neutral_tree<T>(tree, v);
    tomo::print_neutral_tree(neutral);

    // [x] (write to file)
    tomo::serialize_tree(neutral, "data/partitionings/cone_test.toml");

    // [x] (load from file)
    auto reloaded =
        tomo::deserialize_tree<T>("data/partitionings/cone_test.toml");
    std::cout << "RELOADED: \n";
    tomo::print_neutral_tree(*reloaded);

    // [x] convert to non-neutral with larger size
    int large_size = 128;
    auto large_v = tomo::volume<D, T>(large_size);
    auto large_g = tomo::geometry::cone_beam<T>(large_v, large_size, {(T)2.0, (T)2.0},
                                          {large_size, large_size}, (T)4.0, (T)4.0);
    auto kernel = tomo::dim::joseph<3_D, T>(large_v);

    auto large_voxel = tomo::from_neutral_tree<T>(*reloaded, large_v);
    tomo::print_tree(large_voxel);

    // run in parallel with partitioning

    return 0;
}
