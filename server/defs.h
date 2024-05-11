#ifndef DEFS_H_
#define DEFS_H_

#define DT 0.014

#define CHARGE_ELECTRON 1.6e-19
#define MASS_ELECTRON 9.109e-31

#define MAX_VELOCITY 0.1

#define NBRS 8
#define NBRS_2D 8
#define NBRS_3D 26

#define CORNERS 4
#define CORNERS_2D 4
#define CORNERS_3D 8

#define NUM_GHOSTS 4 // charge_density (1) + current_density (3)

#include "charm++.h"
#include "pup.h"
#include <cmath>
#include <cinttypes>

//#define SIM_BOX_LENGTH 50.0
//#define MAX_ITER 1000

#define SUMMARY_FREQ 1000
#define LB_FREQ 17

#define MODZ(x, a) (x + a) % a // modulo with negative wrap around
#define SAFEMOD(x, a) a ? x % a : x
#define ROUND(x) (int)std::floor(x)

#define EPSILON_FREE 8.85418782e-12
#define MU_FREE 1.25663706e-6

#define apply_periodicity(disp_, width_)                                       \
  disp_ = (std::abs(disp_) == width_ - 1) ? ((disp_ > 0) ? -1 : 1) : disp_
#define IDX3_HASH(x_, y_, z_) 3 * 3 * (x_ + 1) + 3 * (y_ + 1) + (z_ + 1)

#define GHOST_IDX(_NBR, _GHOST) _NBR*NUM_GHOSTS+_GHOST

#define PACK_GHOST(_DATA, _INDEX) \
  _DATA[0] = my_cell->charge_densities[_INDEX];    \
  _DATA[1] = my_cell->current_densities[_INDEX].x; \
  _DATA[2] = my_cell->current_densities[_INDEX].y; \
  _DATA[3] = my_cell->current_densities[_INDEX].z;

#define GHOSTARR_TO_VEC(_DATA, _DIR, _OFFSET)                                                \
  vec3(_DATA[GHOST_IDX(_DIR, 3*_OFFSET)], _DATA[GHOST_IDX(_DIR, 3*_OFFSET+1)],                   \
       _DATA[GHOST_IDX(_DIR, 3*_OFFSET+2)])

#define pack_vec3(_DATA, _OFFSET)                                              \
  vec3(_DATA[_OFFSET], _DATA[_OFFSET + 1], _DATA[_OFFSET + 2])

#define unpack_vec3(_v) _v.x, _v.y, _v.z
#define vec3_pair(_DATA, _TYPE, _i, _j)                                       \
  unpack_vec3(_DATA[_i]._TYPE), unpack_vec3(_DATA[_j]._TYPE)

typedef std::pair<int, int> ipair;

// i: x-axis, j: y-axis, k: z-axis
enum directions2D {
  LEFT = 1, // i - 1
  RIGHT = 2, // i + 1
  UP = 4, // j + 1
  DOWN = 8 // j - 1
};

enum directions3D {
  // LEFT: i - 1
  // RIGHT: i + 1
  // UP: k + 1
  // DOWN: k - 1
  BACK = 16, // j + 1
  FRONT = 32 // j - 1
};

enum properties {
  CHARGE_DENSITY = 1,
  CURRENT_DENSITY = 2
};

enum distribution {
  LINEAR,
  SINE,
  GEOMETRIC
};

enum geometry {
  CART2D,
  CART3D
};

/*********************************************************************************************
 3D vector
 *********************************************************************************************/
struct vec3 {
  double x, y, z;

  vec3(double d = 0.0) : x(d), y(d), z(d) {}
  vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

  inline vec3 &operator+=(const vec3 &rhs) {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
  }

  inline vec3 operator+(const vec3 &rhs) {
    return vec3(x + rhs.x, y + rhs.y, z + rhs.z);
  }

  inline vec3 &operator-=(const vec3 &rhs) { return *this += (rhs * -1.0); }

  inline vec3 operator*(const double d) const {
    return vec3(d * x, d * y, d * z);
  }

  inline vec3 operator/(const double d) const {
    return vec3(x/d, y/d, z/d);
  }

  inline vec3 operator-(const vec3 &rhs) const {
    return vec3(x - rhs.x, y - rhs.y, z - rhs.z);
  }
};
inline double dot(const vec3 &a, const vec3 &b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
PUPbytes(vec3)


#endif // DEFS_H_


/*********************************************************************************************
 Charged Particle (with Unit Mass)
 *********************************************************************************************/
struct Particle {
  vec3 position, velocity, acceleration;
  float charge, mass;
  uint8_t geometry;

  Particle(CkMigrateMessage *m) {}

  Particle() {}
  virtual ~Particle(){};

  Particle(double x, double y, double z, double v_x, double v_y, double v_z, float charge, float mass)
      : charge(charge), mass(mass) {
    position.x = x;
    position.y = y;
    position.z = z;
    velocity.x = v_x;
    velocity.y = v_y;
    velocity.z = v_z;
  }

  Particle(double x, double y, double z, float charge, float mass, uint8_t geometry)
    : charge(charge), mass(mass), geometry(geometry) {
    position.x = x;
    position.y = y;
    position.z = z;
    position.z = 0.0;
    velocity.x = 0.0;
    velocity.y = 0.0;
    velocity.z = 0.0;
  }

  void apply_force(const vec3 &force, const uint32_t sim_box_length) {
    acceleration = force / mass;

    velocity += acceleration * DT;
    limit_velocity();

    // assuming periodic boundaries
    position.x = wrap_around(position.x + velocity.x*DT, sim_box_length);
    position.y = wrap_around(position.y + velocity.y*DT, sim_box_length);
    if (geometry == CART3D)
      position.z = wrap_around(position.z + velocity.z*DT, sim_box_length);
  }

  double wrap_around(double t, float w) {
    if (t > w)
      t -= w;
    if (t < 0)
      t += w;

    return t;
  }

  double check_velocity(double in_vel) {
    if (std::fabs(in_vel) > MAX_VELOCITY) {
      if (in_vel < 0.0)
        return -MAX_VELOCITY;
      else
        return MAX_VELOCITY;
    } else {
      return in_vel;
    }
  }

  void limit_velocity() {
    velocity.x = check_velocity(velocity.x);
    velocity.y = check_velocity(velocity.y);
    velocity.z = check_velocity(velocity.z);
  }
};
PUPbytes(Particle)

/*********************************************************************************************
 Forcing Field associated with Cell
 *********************************************************************************************/
struct Field {
  vec3 electric_field;
  vec3 magnetic_field;

  Field(vec3 e_field, vec3 b_field)
      : electric_field(e_field), magnetic_field(b_field) {}

  Field(CkMigrateMessage *m) {}

  Field() {}
  virtual ~Field(){};
};
PUPbytes(Field)

struct PicParams {
  uint8_t pos_distribution, geometry;
  uint32_t odf, initial_particle_count, max_iterations;
  float mass, charge, time_delta, alpha, beta;
  uint32_t sim_box_length, box_count, lb_freq;
  PicParams() {}
  PicParams(uint8_t pos_distribution_, uint8_t geometry_, uint32_t odf_,
            uint32_t initial_particle_count_, uint32_t max_iterations_,
            float time_delta_, uint32_t lb_freq_)
      : pos_distribution(pos_distribution_), geometry(geometry_), odf(odf_),
        initial_particle_count(initial_particle_count_),
        max_iterations(max_iterations_), time_delta(time_delta_), lb_freq(lb_freq_) {
    alpha = 0.0;
    beta = 0.0;
  }
};
PUPbytes(PicParams)

// calculate force using electric and magnetic fields
inline vec3 field_to_force(const Field& field, const Particle& p) {
  vec3 force;
  force = field.electric_field * p.charge; // F_E = qE

  // F_B = q (v cross B)
  vec3 v = p.velocity;

  vec3 cross_product = {
      field.magnetic_field.y * v.z - field.magnetic_field.z * v.y,
      field.magnetic_field.z * v.x - field.magnetic_field.x * v.z,
      field.magnetic_field.x * v.y - field.magnetic_field.y * v.x};

  force += cross_product * p.charge;

  return force;
}
