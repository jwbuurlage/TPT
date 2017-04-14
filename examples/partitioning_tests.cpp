#include "tomos/tomos.hpp"
#include "tomos/util/trees.hpp"
#include "tomos/util/simple_args.hpp"

using T = float;
constexpr tomo::dimension D = 3_D;

int main(int argc, char* argv[]) {
    auto opts = tomo::options{argc, argv};

    int processors = 16;

    if (opts.passed("-p")) {
        processors = opts.arg_as<int>("-p");
    }

    (void)processors;

    return 0;
}
