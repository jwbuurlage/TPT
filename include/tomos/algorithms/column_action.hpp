#pragma once

#include <optional>
#include <random>
#include <string>

#include "../common.hpp"
#include "../geometry.hpp"
#include "../image.hpp"
#include "../math.hpp"
#include "../operations.hpp"
#include "../projections.hpp"
#include "../util/column_iterator.hpp"
#include "../util/image_processing.hpp"
#include "../util/read_tiff.hpp"
#include "../volume.hpp"

using namespace tomo::img;

namespace tomo {

struct index_space {
    virtual int operator()(int /* s */, int i) { return i; }

    virtual ~index_space() = default;
};

struct reverse_index_space : index_space {
    reverse_index_space(int n_) : n(n_) {}
    int n;

    int operator()(int /* s */, int i) override { return n - 1 - i; }
};

struct back_forth_index_space : index_space {
    back_forth_index_space(int n_) : n(n_) {}
    int n;

    int operator()(int s, int i) override {
        return s % 2 == 0 ? i : (n - 1 - i);
    }
};

struct random_index_space : index_space {
    random_index_space(int n_) : n(n_) {
        indices = std::vector<int>(n);
        std::iota(indices.begin(), indices.end(), 0);

        prev_s = -1;
    }

    void reset() {
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(indices.begin(), indices.end(), g);
    }

    int n;
    int prev_s;
    std::vector<int> indices;

    int operator()(int s, int i) override {
        if (prev_s != s) {
            reset();
        }
        prev_s = s;
        return indices[i];
    }
};

namespace reconstruction {

using namespace std::string_literals;

template <dimension D, typename T>
image<D, T> column_action_cyclic(
    volume<D, T>& v, const tomo::geometry::base<D, T>& g,
    tomo::dim::base<D, T>& kernel, const projections<D, T>& b,
    double beta = 0.5, int sweeps = 10,
    std::optional<image<D, T>> x0 = {}, std::optional<index_space*> idxs = {},
    std::function<void(const image<D, T>&, int, const projections<D, T>&)>
        callback = {}) {
    auto x = x0.value_or(tomo::image<D, T>(v, 0));

    tomo::write_png(x, "fan_beam_initial");
    auto ax = tomo::forward_projection(x, g, kernel);

    auto r = b - ax;
    auto column = tomo::column<D, T>(g, kernel);

    auto cs = tomo::image<D, T>(v);
    for (auto j = 0u; j < v.cells(); ++j) {
        auto voxel_index = v.unroll(j);
        for (auto[line_idx, value] : column(voxel_index)) {
            (void)line_idx;
            cs[j] += value * value;
        }
    }

    for (auto k = 0; k < sweeps; ++k) {
        for (auto q = 0u; q < v.cells(); ++q) {
            auto j = q;
            if (idxs) {
                j = (*idxs.value())(k, q);
            }
            auto voxel_index = v.unroll(j);
            if (cs[j] < math::epsilon<T>) {
                continue;
            }
            auto delta = (T)0;
            for (auto[line_idx, value] : column(voxel_index)) {
                delta += value * r[line_idx];
            }
            delta /= cs[j];
            delta *= beta;
            for (auto[line_idx, value] : column) {
                r[line_idx] -= value * delta;
            }
            x[j] += delta;
        }

        if (callback) {
            callback(x, k, r);
        }
    }

    return x;
}

template <dimension D, typename T>
image<D, T>
column_action_block(volume<D, T>& v, const tomo::geometry::base<D, T>& g,
                    tomo::dim::base<D, T>& kernel, const projections<D, T>& b,
                    std::function<std::vector<uint64_t>(uint64_t)> block,
                    uint64_t block_count, double beta = 0.5, int sweeps = 10) {
    auto x = tomo::image<D, T>(v);
    auto r = b;
    auto column = tomo::column<D, T>(g, kernel);

    auto cs = tomo::image<D, T>(v);
    for (auto j = 0u; j < v.cells(); ++j) {
        auto voxel_index = v.unroll(j);
        for (auto[line_idx, value] : column(voxel_index)) {
            (void)line_idx;
            cs[j] += value * value;
        }
    }

    auto delta = std::vector<T>();
    for (auto k = 0; k < sweeps; ++k) {
        for (auto i = 0u; i < block_count; ++i) {
            auto vals = block(i);
            delta.clear();
            delta.resize(vals.size());

            for (auto idx = 0u; idx < vals.size(); ++idx) {
                auto j = vals[idx];
                delta[idx] = (T)0;
                auto voxel_index = v.unroll(j);
                if (cs[j] < math::epsilon<T>) {
                    continue;
                }
                for (auto[line_idx, value] : column(voxel_index)) {
                    delta[idx] += value * r[line_idx];
                }
                delta[idx] /= cs[j];
                delta[idx] *= beta;
            }

            for (auto idx = 0u; idx < vals.size(); ++idx) {
                auto j = vals[idx];
                auto voxel_index = v.unroll(j);
                for (auto[line_idx, value] : column(voxel_index)) {
                    r[line_idx] -= value * delta[idx];
                }
                x[j] += delta[idx];
            }
        }
    }

    return x;
}

} // namespace reconstruction
} // namespace tomo
