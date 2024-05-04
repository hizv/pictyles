#include "defs.h"
#include "interpolate.hpp"
#include "server.decl.h"

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
  int iter, imsg, box_count;
  uint32_t initial_particle_count, max_iter, sim_box_length;
  Field fields[CORNERS];

public:
  double charge_densities[CORNERS];
  vec3 current_densities[CORNERS];
  Cell(CkMigrateMessage *m) {}

  Cell(int box_count_, uint32_t initial_particle_count_, uint32_t max_iter_,
       uint32_t sim_box_length_)
      : iter(0), initial_particle_count(initial_particle_count_),
        max_iter(max_iter_), box_count(box_count_),
        sim_box_length(sim_box_length_) {
    usesAtSync = true;
  }

  void set_particles_proxy(CProxy_Particles particles_array_) {
    particles_array = particles_array_;
  }

  void pup(PUP::er &p) {
    p | iter;
    p | imsg;
    p | particles_array;
    p | box_count;
    p | initial_particle_count;
    p | max_iter;
    p | sim_box_length;
    PUP::PUParray<Field>(p, fields, CORNERS);
    PUP::PUParray<double>(p, charge_densities, CORNERS);
    PUP::PUParray<vec3>(p, current_densities, CORNERS);
  }

  void process_particles();

  void broadcast_corner() {
    double *data = new double[4];
    data[0] = charge_densities[1];
    data[1] = current_densities[1].x;
    data[2] = current_densities[1].y;
    data[3] = current_densities[1].z;
    thisProxy(MODZ(thisIndex.x - 1, box_count), thisIndex.y)
        .recv_properties(iter, RIGHT, NUM_GHOSTS, data);
    thisProxy(thisIndex.x, MODZ(thisIndex.y + 1, box_count))
        .recv_properties(iter, DOWN, NUM_GHOSTS, data);
    thisProxy(MODZ(thisIndex.x - 1, box_count),
              MODZ(thisIndex.y + 1, box_count))
        .recv_properties(iter, DOWN & RIGHT, NUM_GHOSTS, data);
  }

  void update_properties(int direction, int size, double *data) {
    switch (size) {
    case 4: {
      vec3 curr_density = vec3(data[1], data[2], data[3]);
      switch (direction) {
      case DOWN:
        current_densities[0] = curr_density;
        break;
      case RIGHT:
        current_densities[3] = curr_density;
        break;
      case DOWN &RIGHT:
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
      case DOWN &RIGHT:
        charge_densities[2] = data[0];
        break;
      default:
        CkPrintf("impossible!\n");
      }
      break;
    }
  }

  void solve_fields() {
    for (int i = 0; i < CORNERS; i++) {
      vec3 e_field(charge_densities[i] / EPSILON_FREE, 0.0, 0.0);
      vec3 b_field(0.0, 0.0, MU_FREE * current_densities[i].x);
      fields[i] = Field(e_field, b_field);
    }
  }
};

/*********************************************************************************************
 Particles associated with a Cell in a Grid
 *********************************************************************************************/
class Particles : public CBase_Particles {
Particles_SDAG_CODE
private:
  CProxy_Cell cell_array;
  int iter, imsg, box_count;
  uint32_t initial_particle_count, max_iter;
  float box_width;
  CkCallback cb1, cb2;
  ipair neighbours[NBRS];
  double ghost_data[NBRS * NUM_GHOSTS];

public:
  std::vector<Particle> particles;

  void pup(PUP::er &p) {
    p | particles;
    p | iter;
    p | imsg;
    p | cb1;
    p | cb2;
    p | box_count;
    p | box_width;
    p | max_iter;
    p | initial_particle_count;
    PUP::PUParray<ipair>(p, neighbours, NBRS);
    PUP::PUParray<double>(p, ghost_data, NBRS * NUM_GHOSTS);
  }


  Particles(CkMigrateMessage *m) {}
  Particles() {}
  Particles(uint32_t box_count_, uint32_t initial_particle_count_,
            uint32_t max_iter_, uint32_t sim_box_length_)
      : iter(0), initial_particle_count(initial_particle_count_),
        max_iter(max_iter_), box_count(box_count_) {
    usesAtSync = true;
    std::random_device rd;
    std::mt19937 e2(rd());

    box_width = sim_box_length_ / (float)box_count;

    // initialise random generators for coordinates
    double low_x = thisIndex.x * box_width;
    double high_x = low_x + box_width;

    double low_y = thisIndex.y * box_width;
    double high_y = low_y + box_width;

    std::uniform_real_distribution<double> dist_x(low_x, high_x);
    std::uniform_real_distribution<double> dist_y(low_y, high_y);

    particles = std::vector<Particle>(0);

    // Linearly distributed particle count for now
    double alpha = 1.0;
    double beta = 2.0;
    double step = 1.0 / (box_count - 1);
    double total_weight =
        beta * box_count - alpha * 0.5 * step * box_count * (box_count - 1);
    double current_weight = (beta - alpha * step * ((double)thisIndex.x));

    int particle_count =
        initial_particle_count * (current_weight / total_weight);

    for (int i = 0; i < particle_count; i++) {
      double pos_x = dist_x(e2);
      double pos_y = dist_y(e2);
      double rel_x = std::fmod(pos_x, 1.0);
      double rel_y = std::fmod(pos_y, 1.0);

      double r1_sq = rel_y * rel_y + rel_x * rel_x;
      double r2_sq = rel_y * rel_y + (1.0 - rel_x) * (1.0 - rel_x);
      double cos_theta = rel_x / std::sqrt(r1_sq);
      double cos_phi = (1.0 - rel_x) / std::sqrt(r2_sq);

      double base_charge =
          1e4 * UNIT_CHARGE / (DT * DT * (cos_theta / r1_sq + cos_phi / r2_sq));
      double particle_charge =
          (thisIndex.x % 2 == 0) ? base_charge : -1.0 * base_charge;
      particles.push_back(Particle(pos_x, pos_y, particle_charge, UNIT_MASS));
    }

    cb1 = CkCallback(CkReductionTarget(Main, updateMax), pic_proxy);
    cb2 = CkCallback(CkReductionTarget(Main, updateTotal), pic_proxy);

    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        if (i == 0 && j == 0)
          continue; // (0, 0) is chare itself
        int idx = 3 * (j + 1) + (i + 1);
        if (idx == 8)
          idx = 4; // map (1, 1) to (0, 0) to fit in array

        neighbours[idx] = {MODZ(thisIndex.x + i, box_count),
                           MODZ(thisIndex.y + j, box_count)};

        // CkPrintf("%d, %d -> %d, %d\n", thisIndex.x, thisIndex.y,
        // neighbours[idx].first, neighbours[idx].second);
      }
    }

    // thisProxy(thisIndex.x, thisIndex.y).run();
  }

  void set_cell_proxy(const CProxy_Cell &cell_array_) {
    cell_array = cell_array_;
  }

  void deposit_properties() {
    Cell *my_cell = cell_array[thisIndex].ckLocal();
    for (int i = 0; i < CORNERS; i++) {

      my_cell->charge_densities[i] = 0;
      my_cell->current_densities[i] = vec3(0.0);
    }
    for (auto &p : particles) {
      deposit_density2D(box_width, p.position, p.charge,
                        my_cell->charge_densities);

      // Compute current density and deposit
      vec3 current_density = p.position * dot(p.velocity, p.velocity) *
                             p.charge / dot(p.position, p.position);
      deposit_density2Dvec(box_width, p.position, current_density,
                           my_cell->current_densities);
    }
    // CkPrintf("%lf\n", my_cell->charge_densities[1]);
  }

  void broadcast_ghosts() {
    int i = thisIndex.x, j = thisIndex.y;
    Cell *my_cell = cell_array[thisIndex].ckLocal();

    double *data = new double[NUM_GHOSTS];
    // send my corners
    PACK_GHOST(data, 3)
    thisProxy(MODZ(i + 1, box_count), j)
        .receive_ghost(iter, LEFT, NUM_GHOSTS, data);

    PACK_GHOST(data, 0)
    thisProxy(i, MODZ(j - 1, box_count))
        .receive_ghost(iter, UP, NUM_GHOSTS, data);

    PACK_GHOST(data, 2)
    thisProxy(MODZ(i + 1, box_count), MODZ(j - 1, box_count))
        .receive_ghost(iter, UP & LEFT, NUM_GHOSTS, data);
  }

  void process_ghost(int direction, int size, double *data) {
    for (int i = 0; i < size; i++)
      ghost_data[GHOST_IDX(direction, i)] = data[i];
  }

  void reduce_ghosts() {
    Cell *my_cell = cell_array[thisIndex].ckLocal();
    my_cell->charge_densities[1] += ghost_data[GHOST_IDX(LEFT, 0)] +
                                    ghost_data[GHOST_IDX(UP & LEFT, 0)] +
                                    ghost_data[GHOST_IDX(UP, 0)];

    my_cell->current_densities[1] += ARR_TO_VEC(ghost_data, LEFT) +
                                     ARR_TO_VEC(ghost_data, UP & LEFT) +
                                     ARR_TO_VEC(ghost_data, UP);
  }

  // add out of chare particles to the message buffer for the
  // destination chare
  void send_displaced() {
    std::map<ipair, std::vector<Particle>> bufs;

    int i = 0;
    while (i < particles.size()) {
      int dest_chare_x = particles[i].position.x / box_width;
      int dest_chare_y = particles[i].position.y / box_width;

      if (dest_chare_x == thisIndex.x && dest_chare_y == thisIndex.y) {
        i++; // Stay in same chare
      } else {
        // buffer the particle in the destination neighbour's vectors
        ipair dest = {dest_chare_x, dest_chare_y};
        bool flag = true;
        for (ipair &nb : neighbours) {
          if (dest == nb)
            flag = false;
        }
        if (flag)
          CkPrintf("[%f] Dest: %d, %d (%lf, %lf), From: %d, %d\n", box_width, dest.first, dest.second,
                   particles[i].position.x, particles[i].position.y, thisIndex.x, thisIndex.y);
        bufs[dest].push_back(particles[i]);

        particles[i] = particles.back();
        particles.pop_back();
      }
    }

    // send the messages
    for (ipair &nb : neighbours) {
      int len = bufs[nb].size();
      ParticleDataMsg *msg = new (len) ParticleDataMsg(len);

      for (int i = 0; i < len; ++i)
        msg->particles[i] = bufs[nb][i];

      thisProxy(nb.first, nb.second).receive_particles_from_neighbour(msg);
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

  for (int i = 0; i < particle_count; i++) {
    Field field =
        interpolated_field(box_width, particles->at(i).position, fields);
    particles->at(i).apply_force(field_to_force(field, particles->at(i)),
                                 sim_box_length);
  }

  my_particles->cell_done();
}

class PicParams {
public:
  uint8_t ndims;
  uint32_t dims[3];
  uint8_t odf;
  uint32_t initial_particle_count;

  PicParams(uint8_t ndims_, uint32_t *dims_, uint8_t odf_,
            uint32_t initial_particle_count_)
      : ndims(ndims_), odf(odf_),
        initial_particle_count(initial_particle_count_) {}
};

// #include "server.def.h"
