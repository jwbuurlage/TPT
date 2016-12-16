#include "tomo.hpp"
#include "util/plotter.hpp"

using T = float;

int main() {
    // request a plotter scene
    auto plotter =
        tomo::ext_plotter<2_D, T>("tcp://localhost:5555", "Test plot");

    int k = 256;
    auto v = tomo::volume<2_D>(k);
    auto f = tomo::modified_shepp_logan_phantom<T>(v);
    auto g = tomo::geometry::parallel<2_D, T>(k, k, v);

    auto proj = tomo::dim::linear<2_D, T>(v);
    auto dose = tomo::image<2_D, T>(v);

    for (int i = 0; i < 10; ++i) {
        for (auto line : g) {
            auto linear_attenuation = (T)0.0;
            auto intensity = (T)1.0;
            for (auto attenuation_point : proj(line)) {
                auto factor =
                    attenuation_point.value * f[attenuation_point.index];
                linear_attenuation += factor;
                auto new_intensity = tomo::math::exp(-linear_attenuation);
                dose[attenuation_point.index] += intensity - new_intensity;
                intensity = new_intensity;
            }
        }
        plotter.plot(dose);
    }

    return 0;
}
