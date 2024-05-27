#include "pic.hpp"
#include <cinttypes>

template <class T> inline T extract(char *&msg, bool increment = true) {
  T arg = *(reinterpret_cast<T *>(msg));
  if (increment)
    msg += sizeof(T);
  return arg;
}

class Server {
public:
  static void operation_handler(char *msg) {
#ifndef NDEBUG
    CkPrintf("Operation handler called\n");
#endif
    char *cmd = msg + CmiMsgHeaderSizeBytes;
    uint8_t name = extract<uint8_t>(cmd);
    uint32_t epoch = extract<uint32_t>(cmd);
    uint32_t cmd_size = extract<uint32_t>(cmd);
#ifndef NDEBUG
    CkPrintf("%" PRIu8 ", %u, %" PRIu8 "\n", name, epoch);
#endif
  }

  static void create_handler(char *msg) {
#ifndef NDEBUG
    CkPrintf("Create handler called\n");
#endif
    char *cmd = msg + CmiMsgHeaderSizeBytes;
    uint8_t sim_box_ndims = extract<uint8_t>(cmd);
    uint32_t sim_box_dims[sim_box_ndims];

    for (int i = 0; i < sim_box_ndims; i++) {
      sim_box_dims[i] = extract<uint32_t>(cmd);
      // CkPrintf("<create_handler> sim_box_dims[%d] = %u\n", i, sim_box_dims[i]);
    }

    // Simulation parameters
    uint8_t odf = extract<uint8_t>(cmd);
    float time_delta = extract<float>(cmd);
    uint32_t iterations = extract<uint32_t>(cmd);
    uint8_t boundary_conditions = extract<uint8_t>(cmd);
    uint8_t geometry = extract<uint8_t>(cmd);
    uint32_t lb_freq = extract<uint32_t>(cmd);
    uint32_t migrate_freq = extract<uint32_t>(cmd);

    // Species parameters
    uint32_t init_particle_count = extract<uint32_t>(cmd);
    uint8_t position_distribution = extract<uint8_t>(cmd);
    float alpha = extract<float>(cmd);
    float beta = extract<float>(cmd);
    float mass = extract<float>(cmd);
    float charge = extract<float>(cmd);

    // CkPrintf("<create_handler> %u %u %f %u mass=%E\n", sim_box_ndims, odf, time_delta,
             // init_particle_count, mass);

    PicParams pp = PicParams(position_distribution, geometry, odf,
                             init_particle_count, iterations, time_delta, lb_freq, migrate_freq);

    pp.mass = mass; pp.charge = charge;
    pp.alpha = alpha; pp.beta = beta;

    create_pic(sim_box_ndims, sim_box_dims, pp);
    uint8_t ret = 1;
    CcsSendReply(1, (void *)&ret);
  }

  static void delete_handler(char *msg) {
    char *cmd = msg + CmiMsgHeaderSizeBytes;
    uint8_t name = extract<uint8_t>(cmd);
  }

  inline static uint8_t get_client_id() { return 0; }

  static void connection_handler(char *msg) {
    uint8_t client_id = get_client_id();
    CcsSendReply(1, (void *)&client_id);
  }

  static void disconnection_handler(char *msg) {
    char *cmd = msg + CmiMsgHeaderSizeBytes;
    uint8_t client_id = extract<uint8_t>(cmd);
#ifndef NDEBUG
    CkPrintf("Disconnected %" PRIu8 " from server\n", client_id);
#endif
  }

  inline static void exit_server(char *msg) { CkExit(); }

  static void create_pic(uint32_t ndims, uint32_t *dims, PicParams pp) {
    uint32_t num_chare_x, num_chare_y, num_chare_z;

    pp.sim_box_length = dims[0];
    num_chare_x = num_chare_y = num_chare_z = 1;

    uint32_t total_chares = pp.odf * CkNumPes();

    if (ndims == 1)
      num_chare_x = total_chares; // std::min(odf, dims[0]);
    if (ndims == 2) {
      num_chare_x = sqrt(total_chares); // std::min(odf, dims[0]);
      num_chare_y = sqrt(total_chares); // std::min(odf, dims[1]);
    }
    if (ndims == 3) {
      num_chare_x = cbrt(total_chares); // std::min(odf, dims[0]);
      num_chare_y = cbrt(total_chares); // std::min(odf, dims[1]);
      num_chare_z = cbrt(total_chares); // std::min(odf, dims[1]);
    }

    pp.box_count = num_chare_x;
    pic_proxy.run_pic(pp);
  }
};
