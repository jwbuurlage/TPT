#include "tomos/tomos.hpp"

using T = float;
constexpr tomo::dimension D = 2_D;

int main(int argc, char* argv[]) {
    if (argc > 0) {
        auto f = tomo::tiff_to_image<2_D, T>(argv[1]);
        tomo::ascii_plot(f);
    }

    return 0;
}
