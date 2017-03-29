#pragma once

#include <functional>
#include <iostream>
#include <thread>
#include <vector>

#include <zmq.hpp>

#include "tomovis.hpp"

#include "../common.hpp"
#include "../image.hpp"
#include "../utilities.hpp"
#include "../volume.hpp"
#include "geometries/trajectory.hpp"
#include "reconstructor.hpp"

namespace tomo {
namespace util {

template <dimension D, typename T>
class ext_plotter;

template <dimension D, typename T>
std::vector<unsigned char> pack_image(image<D, T> f) {
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
std::vector<unsigned char> downsample_pack_image(image<3_D, T>& f,
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
            }
        }
    }

    return grayscale_image;
}

template <dimension D>
class ext_plotter_base {
  public:
    ext_plotter_base() : context_(1), socket_(context_, ZMQ_REQ) {}

    void connect(std::string address, std::string name) {
        using namespace std::chrono_literals;

        // set socket timeout to 200 ms
        socket_.setsockopt(ZMQ_LINGER, 200);
        socket_.connect(address);

        // FIXME see if plotter is up
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
    int scene_id_ = -1;
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
class ext_plotter<3_D, T> : public ext_plotter_base<3_D>,
                            public data_update_observer<T> {
  public:
    ext_plotter()
        : ext_plotter_base<3_D>(), subscribe_socket_(context_, ZMQ_SUB) {
        slices_.push_back(std::make_pair(0, math::slice<T>(0)));
        slices_.push_back(std::make_pair(1, math::slice<T>(1)));
        slices_.push_back(std::make_pair(2, math::slice<T>(2)));
    }

    ext_plotter(std::string address, std::string name = "Anonymous")
        : ext_plotter() {
        connect(address, name);
    }

    virtual ~ext_plotter() {
        if (serve_thread_.joinable()) {
            serve_thread_.join();
        }
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
            scene_id_, volume_size,
            std::move(downsample_pack_image(f, volume_size)));

        vol_packet.send(socket_);

        zmq::message_t reply;
        socket_.recv(&reply);
    }

    void set_reconstructor(std::shared_ptr<on_demand_reconstructor<T>> recon) {
        reconstructor_ = recon;
        recon->add_observer(this);
    }

    void new_data_available(on_demand_reconstructor<T>&) override {
        std::cout << "new data available\n";
        this->send_slices_();
        this->send_volume_();
    }

    void subscribe(std::string subscribe_host) {
        if (scene_id_ < 0) {
            std::cout << "subscribing before having a scene id... returning\n";
            return;
        }

        // set socket timeout to 200 ms
        socket_.setsockopt(ZMQ_LINGER, 200);

        //  Socket to talk to server
        subscribe_socket_.connect(subscribe_host);

        int filter[] = {(std::underlying_type<tomovis::packet_desc>::type)
                            tomovis::packet_desc::set_slice,
                        scene_id_};
        subscribe_socket_.setsockopt(ZMQ_SUBSCRIBE, filter,
                                     sizeof(decltype(filter)));

        int remove_filter[] = {
            (std::underlying_type<tomovis::packet_desc>::type)
                tomovis::packet_desc::remove_slice,
            scene_id_};

        subscribe_socket_.setsockopt(ZMQ_SUBSCRIBE, remove_filter,
                                     sizeof(decltype(remove_filter)));
    }

    void serve() {
        serve_thread_ = std::thread([&]() {
            while (true) {
                zmq::message_t update;
                subscribe_socket_.recv(&update);

                auto desc = ((tomovis::packet_desc*)update.data())[0];
                auto buffer =
                    tomovis::memory_buffer(update.size(), (char*)update.data());

                switch (desc) {
                case tomovis::packet_desc::set_slice: {
                    auto packet = std::make_unique<tomovis::SetSlicePacket>();
                    packet->deserialize(std::move(buffer));

                    int update_slice_index = -1;
                    int i = 0;
                    for (auto& id_and_slice : slices_) {
                        if (id_and_slice.first == packet->slice_id) {
                            update_slice_index = i;
                            break;
                        }
                        ++i;
                    }

                    // FIXME mutex for slices_
                    auto new_slice = math::slice<T>(packet->orientation);
                    if (update_slice_index >= 0) {
                        slices_[update_slice_index] =
                            std::make_pair(packet->slice_id, new_slice);
                    } else {
                        slices_.push_back(
                            std::make_pair(packet->slice_id, new_slice));
                    }

                    send_slices_();

                    break;
                }
                case tomovis::packet_desc::remove_slice: {
                    auto packet =
                        std::make_unique<tomovis::RemoveSlicePacket>();
                    packet->deserialize(std::move(buffer));

                    auto to_erase = std::find_if(
                        slices_.begin(), slices_.end(),
                        [&](auto x) { return x.first == packet->slice_id; });
                    slices_.erase(to_erase);
                    break;
                }
                default:
                    break;
                }
            }
        });
    }

    // FIXME
    // void send_partition_information() {

    // }

    void send_projection_data(
        const geometry::trajectory<3_D, T>& acquisition_geometry,
        const projections<3_D, T>& proj_stack, volume<3_D, T> volume) {
        // send geometry specification and wait for reply
        auto geo_spec_packet = tomovis::GeometrySpecificationPacket(
            scene_id_, false, acquisition_geometry.projection_count());
        geo_spec_packet.volume_min_point =
            math::vec_to_array<3_D, float>(volume.origin());
        geo_spec_packet.volume_max_point = math::vec_to_array<3_D, float>(
            volume.origin() + volume.physical_lengths());
        geo_spec_packet.send(socket_);

        zmq::message_t reply;
        socket_.recv(&reply);

        // send the projections
        for (int i = 0; i < acquisition_geometry.projection_count(); ++i) {
            int projection_id = i;
            std::array<float, 3> source_position =
                math::vec_to_array<3_D, float>(
                    acquisition_geometry.source_location(i));

            std::array<float, 9> detector_orientation;

            // detector_tilt
            auto detector_tilt = acquisition_geometry.detector_tilt(i);
            detector_orientation[0] = (T)2.0 * detector_tilt[0][0];
            detector_orientation[1] = (T)2.0 * detector_tilt[0][1];
            detector_orientation[2] = (T)2.0 * detector_tilt[0][2];
            detector_orientation[3] = (T)2.0 * detector_tilt[1][0];
            detector_orientation[4] = (T)2.0 * detector_tilt[1][1];
            detector_orientation[5] = (T)2.0 * detector_tilt[1][2];

            // detector_location
            std::array<float, 3> detector_position =
                math::vec_to_array<3_D, float>(
                    acquisition_geometry.detector_location(i));
            detector_orientation[6] =
                detector_position[0] -
                (T)0.5 * (detector_orientation[0] + detector_orientation[3]);
            detector_orientation[7] =
                detector_position[1] -
                (T)0.5 * (detector_orientation[1] + detector_orientation[4]);
            detector_orientation[8] =
                detector_position[2] -
                (T)0.5 * (detector_orientation[2] + detector_orientation[5]);

            std::array<int, 2> detector_shape = math::vec_to_array<2_D, int>(
                acquisition_geometry.detector_shape());
            std::vector<unsigned char> data =
                pack_image(proj_stack.get_projection(i));

            auto projection_packet = tomovis::ProjectionDataPacket(
                scene_id_, projection_id, source_position, detector_orientation,
                detector_shape, data);
            projection_packet.send(socket_);
            zmq::message_t proj_reply;
            socket_.recv(&proj_reply);
        }
    }

  private:
    void send_slices_() {
        if (!reconstructor_) {
            std::cout << "ERROR: Can not send slices without reconstructor.\n";
            return;
        }

        for (auto id_slice : slices_) {
            auto the_id = id_slice.first;
            auto the_slice = id_slice.second;

            auto image_data = reconstructor_->get_slice_data(the_slice);

            std::vector<int> image_size(2_D);
            for (int d = 0; d < 2_D; ++d) {
                image_size[d] = image_data.size(d);
            }

            auto upd_packet =
                tomovis::SliceDataPacket(scene_id_, the_id, image_size,
                                         std::move(pack_image(image_data)));

            upd_packet.send(socket_);

            zmq::message_t reply;
            socket_.recv(&reply);
        }
    }

    void send_volume_() {
        if (!reconstructor_) {
            std::cout << "ERROR: Can not send volume without reconstructor.\n";
            return;
        }

        int downsample_size = 32;
        auto vol_image = reconstructor_->get_volume_data(downsample_size);
        std::vector<int> volume_size(3, downsample_size);
        auto vol_packet = tomovis::VolumeDataPacket(
            scene_id_, volume_size, std::move(pack_image(vol_image)));

        vol_packet.send(socket_);

        zmq::message_t reply;
        socket_.recv(&reply);
    }

    std::shared_ptr<on_demand_reconstructor<T>> reconstructor_;
    std::vector<std::pair<int, math::slice<T>>> slices_;

    std::thread serve_thread_;
    zmq::socket_t subscribe_socket_;
};

} // namespace util
} // namespace tomo
