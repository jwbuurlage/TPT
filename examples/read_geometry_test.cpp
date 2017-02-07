#include "tomo.hpp"

int main() {
    try {
        auto geom_and_volume =
            tomo::read_configuration<3_D, float>("data/geometries/parallel.toml");
        auto v = geom_and_volume.second;
        std::cout << "origin: " << glm::to_string(v.origin()) << "\n";
        std::cout << "lengths: " << glm::to_string(v.physical_lengths())
                  << "\n";
        std::cout << "voxels: " << glm::to_string(v.voxels()) << "\n";
    } catch (const std::exception& e) {
        std::cout << "Reading configuration failed: " << e.what() << "\n";
    }

    return 0;
}
