#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/QR>
#include <eigen3/Eigen/SVD>

#include "../geometry.hpp"
#include "../image.hpp"
#include "../volume.hpp"

#include <random>

namespace tpt::util {

using namespace Eigen;

/**
 * Implementation of truncated SVD.
 *
 * Returns:
 *  - A tuple with three matrices
 *    - U (left singular vectors, singular data)
 *    - SIGMA (diagonal matrix with 10 largest singular values)
 *    - V (right singular vectors, singular images)
 *
 * Reference:
 * - <https://research.fb.com/fast-randomized-svd/>
 * - <https://arxiv.org/pdf/0909.4061.pdf>
 */
template <dimension D, typename T>
std::tuple<std::vector<projections<D, T>>, std::vector<T>,
           std::vector<image<D, T>>>
truncated_svd(geometry::base<D, T>& g, dim::base<D, T>& k, volume<D, T> v,
              int r, volume<D, T> lv, int reprojections = 2) {
    (void)reprojections;
    assert(reprojections == 2);

    // desired rank must be positive
    assert(r > 0);

    // 2. MASK THE RECONSTRUCTION ON THE LOCAL VOLUME
    auto mask = [&](image<D, T> img) -> image<D, T> {
        auto imgsize = img.get_volume().voxels()[0];
        for (int i = lv.origin()[0] * imgsize;
             i < (lv.origin()[0] + lv.physical_lengths()[0]) * imgsize; ++i) {
            for (int j = lv.origin()[1] * imgsize;
                 j < (lv.origin()[1] + lv.physical_lengths()[1]) * imgsize;
                 ++j) {
                img({i, j}) = (T)0;
            }
        }
        return img;
    };

    std::vector<T> sigma(r);

    std::random_device rd;
    std::mt19937 gen(rd());
    auto rng = std::uniform_real_distribution<T>((T)0, (T)1);

    // 1) generate `rank` random vectors
    auto Omega = std::vector<image<D, T>>(r, image<D, T>(v));
    auto AOmega = std::vector<projections<D, T>>();

    for (auto& img : Omega) {
        std::generate(img.begin(), img.end(), [&]() { return rng(gen); });
    }

    // 2) project everything in Omega to get a vector of projections
    // corresponding to Y
    for (auto& img : Omega) {
        auto y = forward_projection<D, T>(mask(img), g, k);
        auto bandf = [&](const auto& p) {
            auto x = mask(back_projection<D, T>(p, g, k, v));
            return forward_projection<D, T>(mask(x), g, k);
        };
        AOmega.push_back(bandf(bandf(y)));
    }

    // 3) Turn AOmega into a Eigen matrix Y
    auto Y = MatrixXd(g.lines(), r);
    for (auto i = 0u; i < g.lines(); ++i) {
        for (int j = 0; j < r; ++j) {
            Y(i, j) = AOmega[j][i];
        }
    }

    // 4) Compute QR factorization
    auto Q = HouseholderQR<decltype(Y)>(Y).householderQ() *
             MatrixXd::Identity(g.lines(), r);

    for (auto i = 0u; i < g.lines(); ++i) {
        for (int j = 0; j < r; ++j) {
            AOmega[j][i] = Q(i, j);
        }
    }

    // 5) Compute Omega <- A^T Q
    for (int i = 0; i < r; ++i) {
        Omega[i] = back_projection(AOmega[i], g, k, v);
    }

    // 6) Store Omega^T as Eigen matrix B
    auto B = MatrixXd(r, v.cells());
    for (int i = 0; i < r; ++i) {
        for (auto j = 0u; j < v.cells(); ++j) {
            B(i, j) = Omega[i][j];
        }
    }

    // 7) Compute SVD of B
    auto svd = BDCSVD(B, ComputeThinU | ComputeThinV);
    auto U = svd.matrixU();
    auto S = svd.singularValues();
    auto V = svd.matrixV();

    U = Q * U;

    for (int i = 0; i < r; ++i) {
        for (auto j = 0u; j < g.lines(); ++j) {
            AOmega[i][j] = U(j, i);
        }
        sigma[i] = S[i];
        for (auto j = 0u; j < v.cells(); ++j) {
            Omega[i][j] = V(j, i);
        }
    }

    return {std::move(AOmega), std::move(sigma), std::move(Omega)};
}

} // namespace tpt::util
