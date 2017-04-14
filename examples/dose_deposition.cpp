#include "tomos/tomos.hpp"

using T = float;

int main() {
    int k = 256;
    auto v = tomo::volume<2_D, T>(k);
    auto f = tomo::modified_shepp_logan_phantom<T>(v);
    auto g = tomo::geometry::parallel<2_D, T>(v, k, k);

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
        tomo::ascii_plot(dose);
    }

    return 0;
}
