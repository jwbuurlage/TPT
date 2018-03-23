#include <fstream>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "tomos/tomos.hpp"
#include "tomos/util/simple_args.hpp"
#include "tomos/util/trees.hpp"

using T = float;

int main(int argc, char* argv[]) {
    auto opts = tomo::options{argc, argv};

    for (auto filename : opts.args("--in")) {
        auto v = tomo::volume<3_D, T>(512);
        auto root = tomo::deserialize_tree<T>(filename, log2(512));
        auto large_voxel = tomo::from_neutral_tree<T>(*root, v);
        auto neutral = tomo::to_neutral_tree<T>(large_voxel, v);
        tomo::serialize_tree<T>(neutral, filename + "_modified");
    }

    return 0;
}
