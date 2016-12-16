#pragma once

#include <functional>

#include <zmq.hpp>

#include "packets.hpp"

#include "../common.hpp"
#include "../image.hpp"

namespace tomo {

template <dimension D, typename T>
class ext_plotter {
  public:
    ext_plotter(std::string address, std::string name = "Anonymous")
        : context_(1), socket_(context_, ZMQ_REQ) {
        socket_.connect(address);
        // fixme see if plotter is up
        auto packet = tomovis::MakeScenePacket(name);

        // FIXME make part of packet
        packet.send(socket_);


        //  Get the reply.
        zmq::message_t reply;
        socket_.recv(&reply);
        scene_id_ = *(int*)reply.data();
    }

    void plot(tomo::image<D, T> f) {
        std::vector<int> image_size(D);
        for (int d = 0; d < D; ++d) {
            image_size[d] = f.size(d);
        }

        std::vector<unsigned char> grayscale_image(f.get_volume().cells());

        T max = (T)0;
        for (int k = 0; k < f.get_volume().cells(); ++k) {
            if (f[k] > max) {
                max = f[k];
            }
        }

        for (int k = 0; k < f.get_volume().cells(); ++k) {
            grayscale_image[k] = (unsigned char)((T)255 * std::max(f[k], (T)0) / max);
        }

        auto upd_packet = tomovis::UpdateImagePacket(
            scene_id_, image_size, std::move(grayscale_image));

        upd_packet.send(socket_);

        zmq::message_t reply;
        socket_.recv(&reply);
    }

  private:
    zmq::context_t context_;
    zmq::socket_t socket_;
    int scene_id_;
};

} // namespace tomo
