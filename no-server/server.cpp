#include "server.hpp"
#include "conv-ccs.h"
#include "converse.h"
// #include "server.decl.h"

class Main : public CBase_Main {
private:
  int reductions, prints;
  int max_particles, total_particles;
  uint32_t max_iter;
  double walltime_start;

public:
  Main(CkArgMsg *msg) : prints(0), reductions(0) {
    // register_handlers();
#ifndef NDEBUG
    CkPrintf("Initialization done\n");
#endif
    pic_proxy = thisProxy; // for reductions
    PicParams pp = PicParams(LINEAR, CART3D, 64, 1000, 1000, DT, 14, 1);
    pp.alpha = 1; pp.beta = 2; pp.mass = 1; pp.charge = 1;
    pp.sim_box_length = 10;
    pp.box_count = 5;
    run_pic(pp);
  }

  void run_pic(PicParams params) {
    max_iter = params.max_iterations;
    int box_count_z = (params.geometry == CART2D) ? 1 : params.box_count;
    CkPrintf("Running CharmPIC on %d processors with %d chares\n", CkNumPes(),
             params.box_count * params.box_count * box_count_z);

    CkPrintf("Migrate Freq = %d\n", params.migrate_freq);
    // Create array of particles.
    CkArrayOptions opts(params.box_count, params.box_count, box_count_z);

    CProxy_Particles particles_array =
        CProxy_Particles::ckNew(params, opts);

    // Create array of grids.
    opts.bindTo(particles_array);
    CProxy_Cell cell_array =
        CProxy_Cell::ckNew(params, opts);

    particles_array.set_cell_proxy(cell_array);
    cell_array.set_particles_proxy(particles_array);

    walltime_start = CkWallTimer();
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

  // void updateMax(int new_max) {
  //   max_particles = new_max;
  //   if (++reductions == 2)
  //     print_summary();
  // }
  // void updateTotal(int new_total) {
  //   total_particles = new_total;
  //   if (++reductions == 2)
  //     print_summary();
  // }

  void print_summary() {
    // reductions = 0;
    // prints++;
    // CkPrintf("[%d] Max Particles: %d, Total Particles: %d\n", prints,
    //          max_particles, total_particles);
    // if (prints == max_iter / SUMMARY_FREQ) {
      CkPrintf("Execution Time: %lf\n", CkWallTimer() - walltime_start);
      CkExit();
    // }
  }
};

#include "server.def.h"
