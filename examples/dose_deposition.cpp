#include "tomo.hpp"

using T = float;

int main() {
    int k = 32;
    auto v = tomo::volume<3_D>(k);
    auto f = tomo::modified_shepp_logan_phantom<T>(v);
    auto g = tomo::geometry::cone_beam<T>(v, k, 1.0, {k, k}, 2.0, 2.0);

    auto proj = tomo::dim::linear<3_D, T>(v);
    auto dose = tomo::image<3_D, T>(v);

    tomo::ascii_plot(f, 10, 1);

    for (auto line : g) {
        auto linear_attenuation = (T)0.0;
        auto intensity = (T)1.0;
        for (auto attenuation_point : proj(line)) {
            auto factor = attenuation_point.value * f[attenuation_point.index];
            linear_attenuation += factor;
            auto new_intensity = tomo::math::exp(-linear_attenuation);
            dose[attenuation_point.index] += intensity - new_intensity;
            intensity = new_intensity;
        }
    }

    tomo::ascii_plot(dose, 10, 1);
}
