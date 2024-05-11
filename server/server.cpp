#include "server.hpp"
#include "conv-ccs.h"
#include "converse.h"
// #include "server.decl.h"

class Main : public CBase_Main {
private:
  int reductions, prints;
  int max_particles, total_particles;
  uint32_t max_iter;

public:
  Main(CkArgMsg *msg) : prints(0), reductions(0) {
    register_handlers();
#ifndef NDEBUG
    CkPrintf("Initialization done\n");
#endif
    pic_proxy = thisProxy; // for reductions
  }

  void pup(PUP::er &p) {
    p | reductions;
    p | prints;
  }

  void run_pic(PicParams params) {
    max_iter = params.max_iterations;
    CkPrintf("Running CharmPIC on %d processors with %d chares\n", CkNumPes(),
             params.box_count * params.box_count);

    // Create array of particles.
    CkArrayOptions opts(params.box_count, params.box_count, 1);

    CProxy_Particles particles_array =
        CProxy_Particles::ckNew(params, opts);

    // Create array of grids.
    opts.bindTo(particles_array);
    CProxy_Cell cell_array =
        CProxy_Cell::ckNew(params, opts);

    particles_array.set_cell_proxy(cell_array);
    cell_array.set_particles_proxy(particles_array);
    cell_array.run();
    particles_array.run();
  }

  void register_handlers() {
    CcsRegisterHandler("pic_connect", (CmiHandler)Server::connection_handler);
    CcsRegisterHandler("pic_disconnect",
                       (CmiHandler)Server::disconnection_handler);
    // CcsRegisterHandler("pic_operation", (CmiHandler)
    // Server::operation_handler); CcsRegisterHandler("pic_sync", (CmiHandler)
    // Server::sync_handler);
    CcsRegisterHandler("pic_delete", (CmiHandler)Server::delete_handler);
    CcsRegisterHandler("pic_exit", (CmiHandler)Server::exit_server);
    CcsRegisterHandler("pic_create", (CmiHandler)Server::create_handler);
  }

  void updateMax(int new_max) {
    max_particles = new_max;
    if (++reductions == 2)
      print_summary();
  }
  void updateTotal(int new_total) {
    total_particles = new_total;
    if (++reductions == 2)
      print_summary();
  }

  void print_summary() {
    reductions = 0;
    prints++;
    CkPrintf("[%d] Max Particles: %d, Total Particles: %d\n", prints,
             max_particles, total_particles);
    if (prints == max_iter / SUMMARY_FREQ)
      CkExit();
  }
};

#include "server.def.h"
