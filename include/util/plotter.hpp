#pragma once

#include <functional>

#include <zmq.hpp>

#include "packets.hpp"

#include "../common.hpp"
#include "../image.hpp"

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
    ext_plotter(std::string address, std::string name = "Anonymous")
        : ext_plotter_base<3_D>() {
        connect(address, name);
    }

    void plot(tomo::image<3_D, T> f) {
        for (int axis = 0; axis < 3; ++axis) {
            auto slice = tomo::slice(f, f.get_volume()[axis] / 2, axis);

            std::vector<int> image_size(2_D);
            for (int d = 0; d < 2_D; ++d) {
                image_size[d] = slice.size(d);
            }

            std::cout << "size: " << image_size[0] << ", " << image_size[1]
                      << "\n";

            auto upd_packet = tomovis::SliceDataPacket(
                scene_id_, axis, image_size, std::move(pack_image(slice)));

            upd_packet.send(socket_);

            zmq::message_t reply;
            socket_.recv(&reply);
        }
    }
};

} // namespace tomo
