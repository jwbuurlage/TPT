#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <thread>
#include <vector>

#include "bulk/bulk.hpp"
namespace bulk {
using namespace experimental;
}
#include "tomop/tomop.hpp"
#include <zmq.hpp>

#include "../common.hpp"
#include "../distributed/restricted_geometry.hpp"
#include "../geometries/trajectory.hpp"
#include "../image.hpp"
#include "../utilities.hpp"
#include "../volume.hpp"
#include "reconstructor.hpp"

namespace tpt {
namespace util {

template <dimension D, typename T>
class ext_plotter;

template <dimension D, typename T>
std::vector<float> pack_image(image<D, T> f) {
    std::vector<float> grayscale_image(f.get_volume().cells());

    for (auto k = 0u; k < f.get_volume().cells(); ++k) {
        grayscale_image[k] = (float)f[k];
    }

    return grayscale_image;
}

template <typename T>
std::vector<float> downsample_pack_image(image<3_D, T>& f,
                                         std::array<int32_t, 3> volume_size) {
    assert(volume_size.size() == 3);

    std::vector<float> grayscale_image(volume_size[0] * volume_size[1] *
                                       volume_size[2]);

    auto voxels = f.get_volume().voxels();
    auto stride = math::vec3<T>{voxels[0] / (T)volume_size[0],
                                voxels[1] / (T)volume_size[1],
                                voxels[2] / (T)volume_size[2]};

    int idx = 0;
    for (int k = 0; k < volume_size[2]; ++k) {
        for (int j = 0; j < volume_size[1]; ++j) {
            for (int i = 0; i < volume_size[0]; ++i) {
                grayscale_image[idx++] =
                    (float)(f({(int)(i * stride[0]), (int)(j * stride[2]),
                               (int)(k * stride[2])}));
            }
        }
    }

    return grayscale_image;
}

template <dimension D>
class ext_plotter_base {
  public:
    ext_plotter_base() : context_(1), socket_(context_, ZMQ_REQ) {}

    void connect(std::string address) {
        // set socket timeout to 200 ms
        socket_.setsockopt(ZMQ_LINGER, 200);
        socket_.connect(address);
        // FIXME see if plotter is up
    }

    void make_scene(std::string name) {
        using namespace std::chrono_literals;

        auto packet = tomop::MakeScenePacket(name, D);

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

    int scene_id() const { return scene_id_; }
    void set_scene_id(int scene_id) { scene_id_ = scene_id; }

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
        connect(address);
        make_scene(name);
    }

    void plot(image<2_D, T> f) {
        std::array<int, 2> image_size;
        for (int d = 0; d < 2_D; ++d) {
            image_size[d] = f.size(d);
        }

        auto upd_packet = tomop::SliceDataPacket(
            scene_id_, 0, image_size, false, std::move(pack_image(f)));

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
        : ext_plotter_base<3_D>(), subscribe_socket_(context_, ZMQ_SUB) {}

    ext_plotter(std::string address) : ext_plotter() { connect(address); }

    ext_plotter(std::string address, std::string name) : ext_plotter(address) {
        make_scene(name);
    }

    virtual ~ext_plotter() {
        if (serve_thread_.joinable()) {
            serve_thread_.join();
        }
    }

    void plot(image<3_D, T> f) {
        for (int axis = 0; axis < 3; ++axis) {
            auto slice =
                slice(f, f.get_volume().voxels()[axis] / 2, axis);

            std::array<int, 2> image_size;
            for (int d = 0; d < 2; ++d) {
                image_size[d] = slice.size(d);
            }

            auto upd_packet =
                tomop::SliceDataPacket(scene_id_, axis, image_size, false,
                                       std::move(pack_image(slice)));

            upd_packet.send(socket_);

            zmq::message_t reply;
            socket_.recv(&reply);
        }

        // also send low resolution volume data
        int downsample_size = 32;
        std::array<int32_t, 3> volume_size{downsample_size, downsample_size,
                                           downsample_size};

        auto vol_packet = tomop::VolumeDataPacket(
            scene_id_, volume_size,
            std::move(downsample_pack_image(f, volume_size)));

        vol_packet.send(socket_);

        zmq::message_t reply;
        socket_.recv(&reply);
    }

    void partial_plot(image<3_D, T> f, volume<3_D, T> global,
                      volume<3_D, T> local) {
        // for (int axis = 0; axis < 3; ++axis) {
        //    auto slice =
        //        slice(f, f.get_volume().voxels()[axis] / 2, axis);

        //    std::array<int, 2> image_size;
        //    for (int d = 0; d < 2; ++d) {
        //        image_size[d] = slice.size(d);
        //    }

        //    auto upd_packet =
        //        tomop::PartialSliceDataPacket(scene_id_, axis, image_size,
        //        false,
        //                               std::move(pack_image(slice)));

        //    upd_packet.send(socket_);

        //    zmq::message_t reply;
        //    socket_.recv(&reply);
        //}

        // also send low resolution volume data
        int downsample_size = 64;
        std::array<int32_t, 3> volume_size{downsample_size, downsample_size,
                                           downsample_size};

        auto vol_packet = tomop::VolumeDataPacket(
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

        int filter[] = {(std::underlying_type<tomop::packet_desc>::type)
                            tomop::packet_desc::set_slice,
                        scene_id_};
        subscribe_socket_.setsockopt(ZMQ_SUBSCRIBE, filter,
                                     sizeof(decltype(filter)));

        int remove_filter[] = {(std::underlying_type<tomop::packet_desc>::type)
                                   tomop::packet_desc::remove_slice,
                               scene_id_};

        subscribe_socket_.setsockopt(ZMQ_SUBSCRIBE, remove_filter,
                                     sizeof(decltype(remove_filter)));
    }

    void serve() {
        serve_thread_ = std::thread([&]() {
            while (true) {
                zmq::message_t update;
                subscribe_socket_.recv(&update);

                auto desc = ((tomop::packet_desc*)update.data())[0];
                auto buffer =
                    tomop::memory_buffer(update.size(), (char*)update.data());

                switch (desc) {
                case tomop::packet_desc::set_slice: {
                    auto packet = std::make_unique<tomop::SetSlicePacket>();
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
                case tomop::packet_desc::remove_slice: {
                    auto packet = std::make_unique<tomop::RemoveSlicePacket>();
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

    template <int G>
    void
        send_partition_information(bulk::rectangular_partitioning<3_D, G>& part,
                                   int processors, volume<3_D, T> volume) {
        auto voxels = volume.voxels();
        auto lengths = volume.physical_lengths();
        auto volume_origin = volume.origin();
        for (int s = 0; s < processors; ++s) {
            auto origin = part.origin(s);
            auto local_size = part.local_size(s);

            auto min_pt = math::vec3<float>{(float)origin[0] / voxels[0],
                                            (float)origin[1] / voxels[1],
                                            (float)origin[2] / voxels[2]} *
                              lengths +
                          volume_origin;

            auto max_pt =
                math::vec3<float>{
                    (float)(origin[0] + local_size[0]) / voxels[0],
                    (float)(origin[1] + local_size[1]) / voxels[1],
                    (float)(origin[2] + local_size[2]) / voxels[2]} *
                    lengths +
                volume_origin;

            auto set_part_packet = tomop::SetPartPacket(
                scene_id_, s, math::vec_to_array<3_D, float>(min_pt),
                math::vec_to_array<3_D, float>(max_pt));

            set_part_packet.send(socket_);

            zmq::message_t reply;
            socket_.recv(&reply);
        }
    }

    void send_geometry(const geometry::base<3_D, T>& acquisition_geometry,
                       volume<3_D, T> volume) {
        // send geometry specification and wait for reply
        auto geo_spec_packet = tomop::GeometrySpecificationPacket(
            scene_id_, false, acquisition_geometry.projection_count());
        geo_spec_packet.volume_min_point =
            math::vec_to_array<3_D, float>(volume.origin());
        geo_spec_packet.volume_max_point = math::vec_to_array<3_D, float>(
            volume.origin() + volume.physical_lengths());
        geo_spec_packet.send(socket_);

        zmq::message_t reply;
        socket_.recv(&reply);
    }

    void group_request(int number) {
        auto grp = tomop::GroupRequestSlicesPacket(scene_id_, number);
        grp.send(socket_);
        zmq::message_t reply;
        socket_.recv(&reply);
    }

    void
    send_partial_projection_data(const distributed::restricted_geometry<T>& g,
                                 const projections<3_D, T>& p) {
        for (int i = 0; i < g.projection_count(); ++i) {
            int projection_id = i;
            std::array<float, 3> source_position =
                math::vec_to_array<3_D, float>(g.source_location(i));

            std::array<float, 9> detector_orientation;

            // detector_tilt
            auto detector_tilt = g.global_geometry().projection_delta(i);
            detector_tilt[0] *= (T)g.global_geometry().projection_shape(i)[0];
            detector_tilt[1] *= (T)g.global_geometry().projection_shape(i)[1];
            detector_orientation[0] = detector_tilt[0][0];
            detector_orientation[1] = detector_tilt[0][1];
            detector_orientation[2] = detector_tilt[0][2];
            detector_orientation[3] = detector_tilt[1][0];
            detector_orientation[4] = detector_tilt[1][1];
            detector_orientation[5] = detector_tilt[1][2];

            // detector_location
            std::array<float, 3> detector_position =
                math::vec_to_array<3_D, float>(
                    g.global_geometry().detector_corner(i));
            detector_orientation[6] = detector_position[0];
            detector_orientation[7] = detector_position[1];
            detector_orientation[8] = detector_position[2];
            std::array<int, 2> detector_shape = math::vec_to_array<2_D, int>(
                g.global_geometry().projection_shape(i));
            std::vector<float> data = pack_image(p.get_projection(i));
            auto sh = g.local_shadow(i);
            auto partial_shape = math::vec_to_array<2_D, int>(
                sh.max_pt - sh.min_pt + math::vec<2_D, int>{1, 1});
            auto partial_offset = math::vec_to_array<2_D, int>(sh.min_pt);

            auto projection_packet = tomop::PartialProjectionDataPacket(
                scene_id_, projection_id, source_position, detector_orientation,
                detector_shape, partial_offset, partial_shape, data);
            projection_packet.send(socket_);
            zmq::message_t proj_reply;
            socket_.recv(&proj_reply);
        }
    }

    void
    send_projection_data(const geometry::base<3_D, T>& acquisition_geometry,
                         const projections<3_D, T>& proj_stack,
                         volume<3_D, T> volume) {
        send_geometry(acquisition_geometry, volume);
        // send the projections
        for (int i = 0; i < acquisition_geometry.projection_count(); ++i) {
            int projection_id = i;
            std::array<float, 3> source_position =
                math::vec_to_array<3_D, float>(
                    acquisition_geometry.source_location(i));

            std::array<float, 9> detector_orientation;

            // detector_tilt
            auto detector_tilt = acquisition_geometry.projection_delta(i);
            detector_tilt[0] *= (T)acquisition_geometry.projection_shape(i)[0];
            detector_tilt[1] *= (T)acquisition_geometry.projection_shape(i)[1];
            detector_orientation[0] = detector_tilt[0][0];
            detector_orientation[1] = detector_tilt[0][1];
            detector_orientation[2] = detector_tilt[0][2];
            detector_orientation[3] = detector_tilt[1][0];
            detector_orientation[4] = detector_tilt[1][1];
            detector_orientation[5] = detector_tilt[1][2];

            // detector_location
            std::array<float, 3> detector_position =
                math::vec_to_array<3_D, float>(
                    acquisition_geometry.detector_corner(i));
            detector_orientation[6] = detector_position[0];
            detector_orientation[7] = detector_position[1];
            detector_orientation[8] = detector_position[2];
            std::array<int, 2> detector_shape = math::vec_to_array<2_D, int>(
                acquisition_geometry.projection_shape(i));
            std::vector<float> data = pack_image(proj_stack.get_projection(i));

            auto projection_packet = tomop::ProjectionDataPacket(
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

            std::array<int, 2> image_size;
            for (int d = 0; d < 2; ++d) {
                image_size[d] = image_data.size(d);
            }

            auto upd_packet =
                tomop::SliceDataPacket(scene_id_, the_id, image_size, false,
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
        std::array<int32_t, 3> volume_size{downsample_size, downsample_size,
                                           downsample_size};
        auto vol_packet = tomop::VolumeDataPacket(
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
} // namespace tpt
