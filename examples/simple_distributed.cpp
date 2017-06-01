#include "bulk/backends/cpp/cpp.hpp"
#include "bulk/bulk.hpp"

#include "tomos/tomos.hpp"
#include "tomos/util/plotter.hpp"

using T = float;

/* class restricted_geometry {
    struct window {
        int x_min;
        int x_max;
        int y_min;
        int y_max;
    };

  public:
    restricted_geometry(tomo::geometry::trajectory& geometry,
                        tomo::volume<3_D> limited_volume)
        : geometry_(geometry) {
        project();
    }

    void project() {
        std::vector<tomo::vec3<T>> corners;
        int pixels = geometry_.detector_pixel_count();
        windows_.resize(geometry_.projection_count());
        for (int i = 0; i < geometry_.projection_count(); ++i) {
            for (int c = 0; c < 8; ++c) {
                auto rvec =
                    intersect({geometry_.source_position(i), corners[c]},
                              geometry_.detector_plane(i));
                // update xmin xmax ymin ymax
            }
            // project all corners
            // gather all mins maxes
            windows_[i].xmin = 0;
        }
    }

    int global_number(int line_number) {
     //...
    }

  private:
    tomo::geometry::trajectory& geometry_;
    std::vector<window> windows_;
}; */

template <tomo::dimension D>
void run(tomo::util::args opt) {
    bulk::cpp::environment env;

    // xeon phi in comment
    int x = 2; // 4
    int y = 2; // 17
    int z = 2; // 4

    auto v = tomo::volume<D, T>(opt.k);
    auto g = tomo::geometry::cone_beam<T>(v, opt.k, {1.5, 1.5}, {opt.k, opt.k},
                                          10.0, 2.0);
    auto sino = tomo::projections<D, T>(g);

    env.spawn(x * y * z, [&](auto& world, int s, int p) {
        (void)world;
        (void)p;

        tomo::math::vec3<T> proc;
        proc[0] = s % x;
        proc[1] = (s / x) % y;
        proc[2] = (s / (x * y)) % z;

        tomo::util::ext_plotter<D, T> plotter("tcp://localhost:5555",
                                              "Sequential serve-and-plot");

        tomo::math::vec3<T> block = {(T)1 / x, (T)1 / y, (T)1 / z};
        auto lv = tomo::volume<D, T>({opt.k / x, opt.k / y, opt.k / z},
                                     proc * block, block);
        auto f = tomo::image<D, T>(lv);
        tomo::fill_ellipsoids_(f, tomo::mshl_ellipsoids_<T>(), lv, v);

        auto rg = restricted_geometry(g, lv);
        auto rsino = std::vector<T>(rg.lines());

        auto proj = tomo::dim::joseph<3_D, T>(lv);
        int line_number = 0;
        for (auto line : rg) {
            for (auto elem : proj(line)) {
                rsino[line_number] += f[elem.index] * elem.value;
            }
            ++line_number;
        }

        // now write to sino
        auto sinos = bulk::gather_all(rsino.data());

        plotter.plot(f);
        plotter.send_projection_data(g, sino, lv);
    });
}

int main(int argc, char* argv[]) {
    auto opt = tomo::util::args(argc, argv);
    run<3_D>(opt);

    return 0;
}
