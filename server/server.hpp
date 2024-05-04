#include <cinttypes>
#include "pic.hpp"

template<class T>
inline T extract(char* &msg, bool increment=true)
{
    T arg = *(reinterpret_cast<T*>(msg));
    if (increment)
        msg += sizeof(T);
    return arg;
}


class Server
{
public:
    static void operation_handler(char* msg)
    {
#ifndef NDEBUG
        CkPrintf("Operation handler called\n");
#endif
        char* cmd = msg + CmiMsgHeaderSizeBytes;
        uint8_t name = extract<uint8_t>(cmd);
        uint32_t epoch = extract<uint32_t>(cmd);
        uint32_t cmd_size = extract<uint32_t>(cmd);
#ifndef NDEBUG
        CkPrintf("%" PRIu8 ", %u, %" PRIu8 "\n", name, epoch);
#endif
        // std::string new_cmd = generate_code(cmd, it->second.odf,
        //    it->second.ndims, it->second.num_fields, it->second.dims,
        //    it->second.ghost_depth);
        // //CkPrintf("Epoch for handler %u\n", epoch);
        // st.receive_graph(epoch, new_cmd.size(), new_cmd.c_str());
    }

    static void create_handler(char* msg)
    {
#ifndef NDEBUG
        CkPrintf("Create handler called\n");
#endif
        char* cmd = msg + CmiMsgHeaderSizeBytes;
        uint8_t sim_box_ndims = extract<uint8_t>(cmd);
        uint32_t sim_box_dims[sim_box_ndims];

        for (int i = 0; i < sim_box_ndims; i++) {
            sim_box_dims[i] = extract<uint32_t>(cmd);
            CkPrintf("<create_handler> sim_box_dims[%d] = %u\n", i, sim_box_dims[i]);
        }

        uint8_t odf = extract<uint8_t>(cmd);
        float time_delta = extract<float>(cmd);
        uint32_t iterations = extract<uint32_t>(cmd);
        uint8_t boundary_conditions = extract<uint8_t>(cmd);

        // Species parameters
        uint32_t init_particle_count = extract<uint32_t>(cmd);
        uint8_t position_distribution = extract<uint8_t>(cmd);

        CkPrintf("<create_handler> %u %u %f %u\n", sim_box_ndims, odf, time_delta, init_particle_count);

        create_pic(sim_box_ndims, sim_box_dims, odf, init_particle_count, iterations);
        uint8_t ret = 1;
        CcsSendReply(1, (void*) &ret);
    }

    static void delete_handler(char* msg)
    {
        char* cmd = msg + CmiMsgHeaderSizeBytes;
        uint8_t name = extract<uint8_t>(cmd);
    }

    inline static uint8_t get_client_id() {
        return 0;
    }

    static void connection_handler(char* msg)
    {
        uint8_t client_id = get_client_id();
        CcsSendReply(1, (void*) &client_id);
    }

    static void disconnection_handler(char* msg)
    {
        char* cmd = msg + CmiMsgHeaderSizeBytes;
        uint8_t client_id = extract<uint8_t>(cmd);
#ifndef NDEBUG
        CkPrintf("Disconnected %" PRIu8 " from server\n", client_id);
#endif
    }

    inline static void exit_server(char* msg)
    {
        CkExit();
    }

    static void create_pic(uint32_t ndims, uint32_t* dims,
            uint32_t odf, uint32_t initial_particle_count, uint32_t max_iterations)
    {
        uint32_t num_chare_x, num_chare_y, num_chare_z;

        uint32_t sim_box_length = dims[0];
        num_chare_x = num_chare_y = num_chare_z = 1;

        uint32_t total_chares = odf * CkNumPes();

        if(ndims == 1)
            num_chare_x = total_chares; //std::min(odf, dims[0]);
        if(ndims == 2)
        {
            num_chare_x = sqrt(total_chares); //std::min(odf, dims[0]);
            num_chare_y = sqrt(total_chares); //std::min(odf, dims[1]);
        }
        if(ndims == 3)
        {
            num_chare_x = cbrt(total_chares); //std::min(odf, dims[0]);
            num_chare_y = cbrt(total_chares); //std::min(odf, dims[1]);
            num_chare_z = cbrt(total_chares); //std::min(odf, dims[1]);
        }

#ifndef NDEBUG
        CkPrintf("Creating PiC of size (%u, %u, %u), max_iter=%u, sim_box=%ux%u\n",
                num_chare_x, num_chare_y, num_chare_z, max_iterations, sim_box_length, sim_box_length);
#endif

        pic_proxy.run_pic(initial_particle_count, num_chare_x, max_iterations, sim_box_length);
    }
};
