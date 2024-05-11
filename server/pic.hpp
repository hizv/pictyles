#include "defs.h"
#include "server.decl.h"
#include "initialise.hpp"
#include "deposit.hpp"
#include "interpolate.hpp"
#include "solver.hpp"

#include "pup_stl.h"
#include <cmath>
#include <random>

/* readonly */ CProxy_Main pic_proxy;

/*********************************************************************************************
 Message to carry Particle details
 *********************************************************************************************/
struct ParticleDataMsg : public CMessage_ParticleDataMsg {
  Particle *particles; // list of atoms
  int size;            // length of list

  ParticleDataMsg(int size_) : size(size_) {}

  void pup(PUP::er &p) {
    CMessage_ParticleDataMsg::pup(p);
    p | size;
    PUParray(p, particles, size);
  }
};

/*********************************************************************************************
 Cell of the grid which generates the ElectricField forcing function
 *********************************************************************************************/
class Cell : public CBase_Cell {
Cell_SDAG_CODE
private:
  CProxy_Particles particles_array;
  int iter, imsg, box_count, corners, inbrs;
  PicParams params;
  uint32_t max_iter, sim_box_length;
  Field fields[CORNERS_2D];
  vec3 *B_inext, *B_jnext, *E_iprev, *E_jprev;

public:
  double *charge_densities;
  vec3 *current_densities;
  Cell(CkMigrateMessage *m) {}

  Cell(PicParams pp)
      : iter(0), max_iter(pp.max_iterations), box_count(pp.box_count),
        sim_box_length(pp.sim_box_length), params(pp) {
    usesAtSync = true;

    corners = (pp.geometry == CART2D) ? CORNERS_2D : CORNERS_3D;
    inbrs = (pp.geometry == CART2D) ? NBRS_2D : NBRS_3D;

    charge_densities = new double[corners];
    current_densities = new vec3[corners];

    B_inext = new vec3[corners];
    B_jnext = new vec3[corners];
    E_iprev = new vec3[corners];
    E_jprev = new vec3[corners];
  }

  void set_particles_proxy(CProxy_Particles particles_array_) {
    particles_array = particles_array_;
  }

  void pup(PUP::er &p) {
    p | iter;
    p | imsg;
    p | particles_array;
    p | box_count;
    p | max_iter;
    p | sim_box_length;
    p | params;
    p | corners;
    p | inbrs;
    PUP::PUParray<Field>(p, fields, corners);

    if (p.isPacking()) {
      delete charge_densities;
      delete current_densities;
      delete B_inext;
      delete B_jnext;
      delete E_iprev;
      delete E_jprev;
    }
    if (p.isUnpacking()) {
      charge_densities = new double[corners];
      current_densities = new vec3[corners];
      B_inext = new vec3[corners];
      B_jnext = new vec3[corners];
      E_iprev = new vec3[corners];
      E_jprev = new vec3[corners];
    }
  }

  void process_particles();

  void broadcast_corner() {
    // pack corner
    double *data = new double[4];
    data[0] = charge_densities[1];
    data[1] = current_densities[1].x;
    data[2] = current_densities[1].y;
    data[3] = current_densities[1].z;

    thisProxy(MODZ(thisIndex.x - 1, box_count), thisIndex.y, thisIndex.z)
        .recv_properties(iter, RIGHT, NUM_GHOSTS, data);
    thisProxy(thisIndex.x, MODZ(thisIndex.y + 1, box_count), thisIndex.z)
        .recv_properties(iter, DOWN, NUM_GHOSTS, data);
    thisProxy(MODZ(thisIndex.x - 1, box_count),
              MODZ(thisIndex.y + 1, box_count), thisIndex.z)
        .recv_properties(iter, DOWN|RIGHT, NUM_GHOSTS, data);
  }

  void update_properties(int direction, int size, double *data) {
    switch (size) {
    case 4: {
      vec3 curr_density = vec3(data[1], data[2], data[3]);

      // CkPrintf("update_properties[%d, (%d, %d, %d)]: J = %E, %E, %E\n", iter,
      // thisIndex.x, thisIndex.y, thisIndex.z, data[1], data[2], data[3]);
      switch (direction) {
      case DOWN:
        current_densities[0] = curr_density;
        break;
      case RIGHT:
        current_densities[3] = curr_density;
        break;
      case DOWN|RIGHT:
        current_densities[2] = curr_density;
        break;
      default:
        CkPrintf("impossible!\n");
      }
    }
    case 1:
      switch (direction) {
      case DOWN:
        charge_densities[0] = data[0];
        break;
      case RIGHT:
        charge_densities[3] = data[0];
        break;
      case DOWN|RIGHT:
        charge_densities[2] = data[0];
        break;
      default:
        CkPrintf("impossible!\n");
      }
      break;
    }
  }

  void share_fields() {
    double B_inext[6] = vec3_pair_to_list(fields, magnetic_field, 2, 3);
    double B_jnext[6] = vec3_pair_to_list(fields, magnetic_field, 1, 3);
    double E_iprev[6] = vec3_pair_to_list(fields, electric_field, 0, 1);
    double E_jprev[6] = vec3_pair_to_list(fields, electric_field, 0, 2);

    int i = thisIndex.x, j = thisIndex.y, k = thisIndex.z;
    thisProxy(MODZ(i - 1, box_count), j, k).recv_vec3(iter, RIGHT, 2, B_inext);
    thisProxy(i, MODZ(j - 1, box_count), k).recv_vec3(iter, UP, 2, B_jnext);
    thisProxy(MODZ(i + 1, box_count), j, k).recv_vec3(iter, LEFT, 2, E_iprev);
    thisProxy(i, MODZ(j + 1, box_count), k).recv_vec3(iter, DOWN, 2, E_jprev);
  }

  void update_field_ghosts(int direction, int count, double *data) {
    if (count == 2) {
      switch (direction) {
      case LEFT:
        E_iprev[0] = pack_vec3(data, 0);
        E_iprev[1] = pack_vec3(data, 1);
        break;
      case UP:
        B_jnext[1] = pack_vec3(data, 0);
        B_jnext[3] = pack_vec3(data, 1);
        break;
      case RIGHT:
        B_inext[2] = pack_vec3(data, 0);
        B_inext[3] = pack_vec3(data, 1);
        break;
      case DOWN:
        E_jprev[0] = pack_vec3(data, 0);
        E_jprev[2] = pack_vec3(data, 1);
        break;
      default:
        CkPrintf("<Cell> Error: Unexpected field ghost direction");
      }
    }
  }

  void solve_fields() {
    Solver *solver = SolverFactory::create(params);
    // Solve poisson on first time step
    if (iter == 1) {
      for (int i = 0; i < corners; i++) {
        fields[i] = solver->poisson(charge_densities[i]);
      }
    } else {
      // add fields inside cell to buffer for completion
      B_inext[0] = fields[2].magnetic_field;
      B_inext[1] = fields[3].magnetic_field;

      B_jnext[0] = fields[1].magnetic_field;
      B_jnext[2] = fields[3].magnetic_field;

      E_iprev[2] = fields[0].electric_field;
      E_iprev[3] = fields[1].electric_field;

      E_jprev[1] = fields[0].electric_field;
      E_jprev[3] = fields[2].electric_field;

      for (int i = 0; i < corners; i++) {
        // CkPrintf("[%d] J = %E, %E, %E\n", iter, current_densities[i].x,
        // current_densities[i].y, current_densities[i].z);
        const vec3 context[5] = {current_densities[i], B_inext[i], B_jnext[i],
                                 E_iprev[i], E_jprev[i]};
        solver->field_solver(&fields[i], context);
      }
    }
    delete solver;
  }
};

/*********************************************************************************************
 Particles associated with a Cell in a Grid
 *********************************************************************************************/
class Particles : public CBase_Particles {
Particles_SDAG_CODE
private:
  CProxy_Cell cell_array;
  int iter, imsg, box_count, corners, inbrs;
  uint8_t geometry;
  uint32_t max_iter;
  float box_width;
  double ghost_data[NBRS * NUM_GHOSTS];

public:
  std::vector<Particle> particles;

  void pup(PUP::er &p) {
    p | cell_array;
    p | iter;
    p | imsg;
    p | box_count;
    p | max_iter;
    p | box_width;
    p | geometry;
    p | corners;
    p | inbrs;
    PUP::PUParray<double>(p, ghost_data, inbrs * NUM_GHOSTS);
    p | particles;
  }

  Particles(CkMigrateMessage *m) {}
  Particles() {}

  Particles(PicParams pp)
      : iter(0), max_iter(pp.max_iterations), box_count(pp.box_count),
        geometry(pp.geometry) {
    usesAtSync = true;
    std::random_device rd;
    std::mt19937 e2(rd());

    box_width = pp.sim_box_length / (float)box_count;

    corners = (geometry == CART2D) ? CORNERS_2D : CORNERS_3D;
    inbrs = (geometry == CART2D) ? NBRS_2D : NBRS_3D;

    // initialise random generators for coordinates
    double low_x = thisIndex.x * box_width;
    double high_x = low_x + box_width;

    double low_y = thisIndex.y * box_width;
    double high_y = low_y + box_width;

    double low_z = 0.0, high_z = 0.0;
    if (geometry == CART3D) {
      low_z = thisIndex.z * box_width;
      high_z = low_z + box_width;
    }
    std::uniform_real_distribution<double> dist_x(low_x, high_x);
    std::uniform_real_distribution<double> dist_y(low_y, high_y);
    std::uniform_real_distribution<double> dist_z(low_z, high_z);

    particles = std::vector<Particle>(0);

    int particle_count;
    switch (pp.pos_distribution) {
    case LINEAR:
      particle_count = linear_init2D(pp.initial_particle_count, box_count,
                                     thisIndex.x, pp.alpha, pp.beta);
      break;
    case SINE:
      particle_count =
          sinusoidal_init2D(pp.initial_particle_count, box_count, thisIndex.x);
      break;
    case GEOMETRIC:
      particle_count = geometric_init2D(pp.initial_particle_count, box_count,
                                        thisIndex.x, pp.alpha);
      break;
    }

    for (int i = 0; i < particle_count; i++) {
      double pos_x = dist_x(e2);
      double pos_y = dist_y(e2);
      double pos_z = dist_z(e2);

      double rel_x = std::fmod(pos_x, 1.0);
      double rel_y = std::fmod(pos_y, 1.0);
      double rel_z = std::fmod(pos_z, 1.0);

      double r1_sq = rel_z * rel_z + rel_y * rel_y + rel_x * rel_x;
      double r2_sq = rel_z * rel_z + rel_y * rel_y + (1.0 - rel_x) * (1.0 - rel_x);
      double cos_theta = rel_x / std::sqrt(r1_sq);
      double cos_phi = (1.0 - rel_x) / std::sqrt(r2_sq);

      double base_charge = pp.charge * CHARGE_ELECTRON /
                           (pp.time_delta * pp.time_delta *
                            (cos_theta / r1_sq + cos_phi / r2_sq));
      // dipole assumption
      double particle_charge =
          (thisIndex.x % 2 == 0) ? base_charge : -1.0 * base_charge;
      particles.push_back(
          Particle(pos_x, pos_y, pos_z, particle_charge, pp.mass * MASS_ELECTRON, geometry));
    }

    // thisProxy(thisIndex.x, thisIndex.y).run();
  }

  void set_cell_proxy(const CProxy_Cell &cell_array_) {
    cell_array = cell_array_;
  }

  void deposit_properties() {
    Cell *my_cell = cell_array[thisIndex].ckLocal();
    Depositor depositor = DepositorFactory::create(geometry);

    for (int i = 0; i < corners; i++) {
      my_cell->charge_densities[i] = 0;
      my_cell->current_densities[i] = vec3(0.0);
    }

    for (auto &p : particles) {
      // deposit charge
      depositor.deposit_property<double>(box_width, p.position, p.charge,
                        my_cell->charge_densities);

      // Compute current density and deposit
      vec3 current_density = p.velocity * p.charge;
      depositor.deposit_property<vec3>(box_width, p.position, current_density,
                           my_cell->current_densities);
    }
    // CkPrintf("%lf\n", my_cell->charge_densities[1]);
  }

  void broadcast_ghosts() {
    int i = thisIndex.x, j = thisIndex.y, k = thisIndex.z;
    Cell *my_cell = cell_array[thisIndex].ckLocal();

    double *data = new double[NUM_GHOSTS];

    switch (geometry) {
    case CART2D: {
      int corner_src[3] = {3, 0, 2};
      int dest_offset[3][3] = {{1, 0, 0}, {0, -1, 0}, {1, -1, 0}};
      int dir_src[3] = {LEFT, UP, UP | LEFT};
      for (int ci = 0; ci < corners - 1; ci++) {
        PACK_GHOST(data, corner_src[ci]);
        thisProxy(MODZ(i + dest_offset[ci][0], box_count),
                  MODZ(j + dest_offset[ci][1], box_count),
                  MODZ(k + dest_offset[ci][2], box_count))
            .receive_ghost(iter, dir_src[ci], NUM_GHOSTS, data);
      }
    } break;

    case CART3D: {
      int corner_src[7] = {0, 2, 3, 4, 5, 6, 7};
      int dest_offset[7][3] = {{0, 0, -1}, {0, 1, -1}, {0, 1, 0}, {1, 0, -1},
                               {1, 0, 0},  {1, 1, -1}, {1, 1, 0}};
      int dir_src[7] = {UP,   FRONT | UP,        FRONT,       LEFT | UP,
                        LEFT, LEFT | FRONT | UP, LEFT | FRONT};
      for (int ci = 0; ci < corners - 1; ci++) {
        PACK_GHOST(data, corner_src[ci]);
        thisProxy(MODZ(i + dest_offset[ci][0], box_count),
                  MODZ(j + dest_offset[ci][1], box_count),
                  MODZ(k + dest_offset[ci][2], box_count))
            .receive_ghost(iter, dir_src[ci], NUM_GHOSTS, data);
      }
    } break;
    }
  }

  void process_ghost(int direction, int size, double *data) {
    for (int i = 0; i < size; i++)
      ghost_data[GHOST_IDX(direction, i)] = data[i];
  }

  void reduce_ghosts() {
    Cell *my_cell = cell_array[thisIndex].ckLocal();
    switch (geometry) {
      case CART2D: {
        int dir_src[3] = {LEFT, UP, UP | LEFT};
        for (int &di : dir_src) {
          my_cell->charge_densities[1] += ghost_data[GHOST_IDX(di, 0)];
          my_cell->current_densities[1] += GHOSTARR_TO_VEC(ghost_data, di, 1);
        }
      }
        break;
      case CART3D: {
      int dir_src[7] = {UP,   FRONT | UP,        FRONT,       LEFT | UP,
                        LEFT, LEFT | FRONT | UP, LEFT | FRONT};
        for (int &di : dir_src) {
          my_cell->charge_densities[1] += ghost_data[GHOST_IDX(di, 0)];
          my_cell->current_densities[1] += GHOSTARR_TO_VEC(ghost_data, di, 1);
        }
      }
        break;
    }
  }

  // add out of chare particles to the message buffer for the
  // destination chare
  void send_displaced() {
    // std::unordered_map<idx3, std::vector<Particle>> bufs;
    std::vector<std::vector<Particle>> bufs(27);
    int i = 0;
    while (i < particles.size()) {
      int dest_chare_x = particles[i].position.x / (double)box_width;
      int dest_chare_y = particles[i].position.y / (double)box_width;
      int dest_chare_z = particles[i].position.z / (double)box_width;

      // If just touching the edge, then floating point precision rounds up to
      // the right boundary
      if (dest_chare_x == box_count &&
          std::fmod(particles[i].position.x, box_width) < 1e-5) {
        dest_chare_x -= 1;
      }
      if (dest_chare_y == box_count &&
          std::fmod(particles[i].position.y, box_width) < 1e-5) {
        dest_chare_y -= 1;
      }

      if (dest_chare_x == thisIndex.x && dest_chare_y == thisIndex.y) {
        i++; // Stay in same chare
      } else {
        int disp_x = dest_chare_x - thisIndex.x,
            disp_y = dest_chare_y - thisIndex.y,
            disp_z = dest_chare_z - thisIndex.z;

        // if periodic, borders are 1 unit apart
        apply_periodicity(disp_x, box_count);
        apply_periodicity(disp_y, box_count);
        apply_periodicity(disp_z, box_count);

        // only migrate to immediate neighbours
        if (std::abs(disp_x) > 1 || std::abs(disp_y) > 1 ||
            std::abs(disp_z) > 1) {
          CkPrintf("[iter=%d] Out of bounds error: %lf, %lf, %lf: box_width=%f\n", iter,
                   particles[i].position.x, particles[i].position.y,
                   particles[i].position.z, box_width);
          CkPrintf("[iter=%d] Out of bounds error: %d, %d, %d: [%d, %d, %d -> %d, %d, %d]\n",
                   iter, disp_x, disp_y, disp_z, thisIndex.x, thisIndex.y,
                   thisIndex.z, dest_chare_x, dest_chare_y, dest_chare_z);
          CkExit();
        }
        bufs[IDX3_HASH(disp_x, disp_y, disp_z)].push_back(particles[i]);

        particles[i] = particles.back();
        particles.pop_back();
      }
    }

    // send the messages
    for (int ni = 0; ni < bufs.size(); ni++) {
      int off_x = ni / 9 - 1, off_y = (ni % 9) / 3 - 1, off_z = (ni % 3) - 1;
      switch(geometry) {
        case CART2D:
          if (off_z != 0)
            continue;
          if (off_x == 0 && off_y == 0)
            continue;
        break;
        case CART3D:
          if (off_x == 0 && off_y == 0 && off_z == 0)
            continue;
      }
      int len = bufs[ni].size();
      ParticleDataMsg *msg = new (len) ParticleDataMsg(len);

      for (int pi = 0; pi < len; ++pi)
        msg->particles[pi] = bufs[ni][pi];

      thisProxy(MODZ(thisIndex.x + off_x, box_count),
                MODZ(thisIndex.y + off_y, box_count), thisIndex.z + off_z)
          .receive_particles_from_neighbour(msg);
    }
  }

  void append_particles(ParticleDataMsg *msg) {
    for (int i = 0; i < msg->size; i++) {
      particles.push_back(msg->particles[i]);
    }
  }

  void check_and_contribute() {
    if (iter % SUMMARY_FREQ == 0) {
      int particle_count = particles.size();
      CkCallback cb1 =
          CkCallback(CkReductionTarget(Main, updateMax), pic_proxy);
      CkCallback cb2 =
          CkCallback(CkReductionTarget(Main, updateTotal), pic_proxy);

      contribute(sizeof(int), &particle_count, CkReduction::max_int, cb1);
      contribute(sizeof(int), &particle_count, CkReduction::sum_int, cb2);
    }
  }
};

void Cell::process_particles() {
  // use bound array to get pointer
  Particles *my_particles = particles_array[thisIndex].ckLocal();
  std::vector<Particle> *particles = &my_particles->particles;
  float box_width = sim_box_length / (float)box_count;
  int particle_count = particles->size();

  Interpolator *interpolator = InterpolatorFactory::create(params.geometry);
  // if (particle_count != 0)
  //   CkPrintf("processing %d particles\n", particle_count);
  for (int i = 0; i < particle_count; i++) {
    Field field = interpolator->interpolated_field(
        box_width, particles->at(i).position, fields);
    particles->at(i).apply_force(field_to_force(field, particles->at(i)),
                                 sim_box_length);
  }

  my_particles->cell_done();
}
// #include "server.def.h"
