#ifndef DEFS_H_
#define DEFS_H_

#define DT 0.014

#define UNIT_CHARGE 1.6e-19
#define UNIT_MASS 1.0

#define MAX_VELOCITY 0.01

#define NBRS 8
#define CORNERS 4
#define NUM_GHOSTS 4

#include "pup.h"
#include <cmath>
#include <cinttypes>

//#define SIM_BOX_LENGTH 50.0
//#define MAX_ITER 1000

#define SUMMARY_FREQ 100
#define LB_FREQ 10

#define MODZ(x, a) (x + a) % a // modulo with negative wrap around
#define SAFEMOD(x, a) a ? x % a : x
#define ROUND(x) (int)std::floor(x)

#define EPSILON_FREE 8.85418782e-12
#define MU_FREE 1.25663706e-6

#define GHOST_IDX(_NBR, _GHOST) _NBR*NUM_GHOSTS+_GHOST

#define PACK_GHOST(_DATA, _INDEX) \
  _DATA[0] = my_cell->charge_densities[_INDEX];    \
  _DATA[1] = my_cell->current_densities[_INDEX].x; \
  _DATA[2] = my_cell->current_densities[_INDEX].y; \
  _DATA[3] = my_cell->current_densities[_INDEX].z;

#define ARR_TO_VEC(_DATA, _DIR)                                                \
  vec3(_DATA[GHOST_IDX(_DIR, 1)], _DATA[GHOST_IDX(_DIR, 2)],                   \
       _DATA[GHOST_IDX(_DIR, 3)])

typedef std::pair<int, int> ipair;

enum directions {
  LEFT = 1,
  RIGHT = 2,
  UP = 4,
  DOWN = 8
};

enum properties {
  CHARGE_DENSITY = 1,
  CURRENT_DENSITY = 2
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

  Particle(CkMigrateMessage *m) {}

  virtual void pup(PUP::er &p) {
    p | position;
    p | velocity;
    p | acceleration;
    p | charge;
    p | mass;
  }
  Particle() {}
  virtual ~Particle(){};

  Particle(double x, double y, double v_x, double v_y, float charge, float mass)
      : charge(charge), mass(mass) {
    position.x = x;
    position.y = y;
    position.z = 0.0;
    velocity.x = v_x;
    velocity.y = v_y;
    velocity.z = 0.0;
  }

  Particle(double x, double y, float charge, float mass)
    : charge(charge), mass(mass) {
    position.x = x;
    position.y = y;
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
    check_velocity(velocity.x);
    check_velocity(velocity.y);
    check_velocity(velocity.z);
  }
};

/*********************************************************************************************
 Forcing Field associated with Cell
 *********************************************************************************************/
struct Field {
  vec3 electric_field;
  vec3 magnetic_field;

  Field(vec3 e_field, vec3 b_field)
      : electric_field(e_field), magnetic_field(b_field) {}

  Field(CkMigrateMessage *m) {}

  virtual void pup(PUP::er &p) {
    p | electric_field;
    p | magnetic_field;
  }
  Field() {}
  virtual ~Field(){};
};

