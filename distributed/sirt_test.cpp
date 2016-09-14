#include <random>
#include <set>
#include <vector>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "tomo.hpp"

#include "bulk/backends/mpi/mpi.hpp"
#include "bulk/bulk.hpp"

using T = double;

struct options {
    int k = 0;
    int size = 0;
    int x = 0;
    int y = 0;
    int iterations = 0;
    float beta = 0.0f;
};

int main(int argc, char* argv[]) {
    // Set up the program
    options opt;

    po::options_description desc("Allowed arguments");
    desc.add_options()("help,h", "show help message")(
        "lines,k", po::value<int>(&opt.k)->default_value(256),
        "number of lines per processor")(
        "size,s", po::value<int>(&opt.size)->default_value(100),
        "size of the phantom")(
        "iterations,i", po::value<int>(&opt.iterations)->default_value(10),
        "number of iterations")(
        "beta,b", po::value<float>(&opt.beta)->default_value(0.5f),
        "value for update relaxation");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return -1;
    }

    auto env = bulk::environment<bulk::mpi::provider>();

    // Spawn the distributed environment
    env.spawn(env.available_processors(), [opt](auto world, int s, int p) {
        // distribute the volume, construct \psi. Here, just compute
        // the total number of elements per processor to balance, and use index
        // to resolve.
        auto v = tomo::volume<2_D>(opt.size, opt.size);
        auto phantom = tomo::modified_shepp_logan_phantom<T>(v);
        auto proj = tomo::dim::closest<2_D, T>(v);
        auto image = tomo::image<2_D, T>(v);

        auto g = tomo::random_list_geometry<2_D, T>(opt.k, v);
        auto measurements = tomo::forward_projection(phantom, g, proj);

        int max_local_voxel_count = (v.x() * v.y() + p - 1) / p;
        int& m = max_local_voxel_count;

        auto obtain_owner = [m](int index) -> int { return index / m; };

        auto available_locally = [s, &obtain_owner](int index) -> bool {
            return obtain_owner(index) == s;
        };

        struct element {
            int owner = 0;
            int index = 0;
        };

        // the indices of elements outside of the range of this processor
        std::set<int> extent;
        for (auto l : g) {
            for (auto touch : proj(l)) {
                if (!available_locally(touch.index)) {
                    extent.insert(touch.index);
                }
            }
        }

        // obtain the owners, normally possibly distributed information
        std::vector<element> halo;
        for (auto& idx : extent) {
            halo.push_back({obtain_owner(idx), idx});
        }

        // now we compute R and bC
        std::vector<T> rs(g.lines());
        std::vector<T> cs(v.cells());

        int row = 0;
        for (auto l : g) {
            for (auto touch : proj(l)) {
                rs[row] += touch.value;
                cs[touch.index] += touch.value;
            }
            row++;
        }

        // communicate bc elements in halo
        auto c_queue = bulk::create_queue<int, T>(world);
        for (auto entry : halo) {
            c_queue(entry.owner).send(entry.index, cs[entry.index]);
        }
        world.sync();

        for (auto msg : c_queue) {
            cs[msg.tag] += msg.content;
        }

        // invert the sums
        for (auto& r : rs)
            r = 1.0 / r;
        for (auto& c : cs)
            c = opt.beta / c;

        auto register_queue = bulk::create_queue<int, int>(world);
        for (auto entry : halo) {
            register_queue(entry.owner).send(s, entry.index);
        }

        world.sync();

        std::vector<element> requests;
        for (auto msg : register_queue) {
            requests.push_back({msg.tag, msg.content});
        }

        // iterate until convergence
        tomo::sinogram<2_D, T, decltype(g), decltype(proj)> sinogram_buffer(g);
        tomo::image<2_D, T> image_buffer(v);

        auto min_local_index = s * m;
        auto max_local_index = std::min((s + 1) * m, v.cells());

        auto exchange_queue = bulk::create_queue<int, T>(world);

        for (int k = 0; k < opt.iterations; ++k) {
            sinogram_buffer.clear();

            // request halo of image from remote
            for (auto request : requests) {
                exchange_queue(request.owner)
                    .send(request.index, image[request.index]);
            }

            world.sync();

            for (auto msg : exchange_queue) {
                image[msg.tag] = msg.content;
            }

            // do a forward-projection
            int row_fp = 0;
            for (auto line : g) {
                for (auto elem : proj(line)) {
                    sinogram_buffer[row_fp] += elem.value * image[elem.index];
                }
                ++row_fp;
            }

            // compute R(p - Wx)
            for (int j = 0; j < g.lines(); ++j) {
                sinogram_buffer[j] =
                    (measurements[j] - sinogram_buffer[j]) * rs[j];
            }

            // zero the image buffer
            image_buffer.clear();

            // perform the back-projection
            int row_bp = 0;
            for (auto line : g) {
                for (auto elem : proj(line)) {
                    image_buffer[elem.index] +=
                        elem.value * sinogram_buffer[row_bp];
                }
                ++row_bp;
            }

            // update the image while scaling
            for (int j = min_local_index; j < max_local_index; ++j) {
                image[j] += cs[j] * image_buffer[j];
            }

            // exchange and sum image here
            for (auto entry : halo) {
                c_queue(entry.owner)
                    .send(entry.index,
                          cs[entry.index] * image_buffer[entry.index]);
            }
            world.sync();

            for (auto msg : c_queue) {
                image[msg.tag] += msg.content;
            }
        }

        // add communication queue for image exchanges
        auto image_queue = bulk::create_queue<int, T>(world);

        // image receive queue for processor 0 to show the entire image:
        for (int j = min_local_index; j < max_local_index; ++j) {
            image_queue(0).send(j, image[j]);
        }

        world.sync();
        if (s == 0) {
            for (auto msg : image_queue) {
                image[msg.tag] = msg.content;
            }
            tomo::ascii_plot(image);
        }
    });

    return 0;
}
