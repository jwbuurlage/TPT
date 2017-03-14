#pragma once

#include <functional>
#include <iostream>

#include <zmq.hpp>

#include "tomovis.hpp"

#include "../common.hpp"
#include "../image.hpp"
#include "../utilities.hpp"

namespace tomo {

template <dimension D, typename T>
class ext_plotter;

template <typename T>
std::vector<unsigned char> pack_image(image<2_D, T> f) {
    std::vector<unsigned char> grayscale_image(f.get_volume().cells());

    T max = (T)0;
    for (int k = 0; k < f.get_volume().cells(); ++k) {
        if (f[k] > max) {
            max = f[k];
        }
    }

    for (int k = 0; k < f.get_volume().cells(); ++k) {
        grayscale_image[k] =
            (unsigned char)((T)255 * std::max(f[k], (T)0) / max);
    }

    return grayscale_image;
}

template <typename T>
std::vector<unsigned char> pack_image(image<3_D, T>& f,
                                      std::vector<int> volume_size) {
    assert(volume_size.size() == 3);

    std::vector<unsigned char> grayscale_image(volume_size[0] * volume_size[1] *
                                               volume_size[2]);

    T max = (T)0;
    for (int k = 0; k < f.get_volume().cells(); ++k) {
        if (f[k] > max) {
            max = f[k];
        }
    }

    auto voxels = f.get_volume().voxels();
    auto stride = math::vec3<T>{voxels[0] / (T)volume_size[0],
                                voxels[1] / (T)volume_size[1],
                                voxels[2] / (T)volume_size[2]};

    std::cout << stride[0] << " " << stride[1] << " " << stride[2] << "\n";

    int idx = 0;
    for (int k = 0; k < volume_size[2]; ++k) {
        for (int j = 0; j < volume_size[1]; ++j) {
            for (int i = 0; i < volume_size[0]; ++i) {
                grayscale_image[idx++] =
                    (unsigned char)((T)255 * std::max(f({(int)(i * stride[0]),
                                                         (int)(j * stride[1]),
                                                         (int)(k)*stride[2]}),
                                                      (T)0) /
                                    max);
                std::cout << idx << " " << (int)grayscale_image[idx] << " "
                          << f({(int)(i * stride[0]), (int)(j * stride[1]),
                                (int)(k)*stride[2]})
                          << "\n";
            }
        }
    }

    std::cout << "10 = " << (int)grayscale_image[10] << " / max = " << max
              << "\n";
    return grayscale_image;
}

template <dimension D>
class ext_plotter_base {
  public:
    ext_plotter_base() : context_(1), socket_(context_, ZMQ_REQ) {}

    void connect(std::string address, std::string name) {
        using namespace std::chrono_literals;

        // FIXME merge this and 3d
        // set socket timeout to 200 ms
        socket_.setsockopt(ZMQ_LINGER, 200);
        socket_.connect(address);

        // fixme see if plotter is up
        auto packet = tomovis::MakeScenePacket(name, D);

        packet.send(socket_);

        //  Initialize poll set
        zmq::pollitem_t items[] = {{socket_, 0, ZMQ_POLLIN, 0}};
        auto poll_result = zmq::poll(items, 1, 1000ms);

        if (poll_result <= 0) {
            std::cout << "Could not connect to server\n";
            throw;
        } else {
            //  Get the reply.
            zmq::message_t reply;
            socket_.recv(&reply);
            scene_id_ = *(int*)reply.data();
            std::cout << "Connected with scene id: " << scene_id_ << "\n";
        }
    }

  protected:
    zmq::context_t context_;
    zmq::socket_t socket_;
    int scene_id_;
};

template <typename T>
class ext_plotter<2_D, T> : public ext_plotter_base<2_D> {
  public:
    ext_plotter() = default;

    ext_plotter(std::string address, std::string name = "Anonymous")
        : ext_plotter_base<2_D>() {
        connect(address, name);
    }

    void plot(tomo::image<2_D, T> f) {
        std::vector<int> image_size(2_D);
        for (int d = 0; d < 2_D; ++d) {
            image_size[d] = f.size(d);
        }

        auto upd_packet = tomovis::SliceDataPacket(scene_id_, 0, image_size,
                                                   std::move(pack_image(f)));

        upd_packet.send(socket_);

        zmq::message_t reply;
        socket_.recv(&reply);
    }

  private:
};

template <typename T>
class ext_plotter<3_D, T> : public ext_plotter_base<3_D> {
  public:
    ext_plotter() = default;

    ext_plotter(std::string address, std::string name = "Anonymous")
        : ext_plotter_base<3_D>() {
        connect(address, name);
    }

    void plot(tomo::image<3_D, T> f) {
        for (int axis = 0; axis < 3; ++axis) {
            auto slice =
                tomo::slice(f, f.get_volume().voxels()[axis] / 2, axis);

            std::vector<int> image_size(2_D);
            for (int d = 0; d < 2_D; ++d) {
                image_size[d] = slice.size(d);
            }

            auto upd_packet = tomovis::SliceDataPacket(
                scene_id_, axis, image_size, std::move(pack_image(slice)));

            upd_packet.send(socket_);

            zmq::message_t reply;
            socket_.recv(&reply);
        }

        // lets also send low resoution volume data
        int downsample_size = 32;
        std::vector<int> volume_size(3, downsample_size);

        auto vol_packet = tomovis::VolumeDataPacket(
            scene_id_, volume_size, std::move(pack_image(f, volume_size)));

        vol_packet.send(socket_);

        zmq::message_t reply;
        socket_.recv(&reply);
    }
};

} // namespace tomo
