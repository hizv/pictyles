#ifndef SOLVER_H_
#define SOLVER_H_

struct Solver {
  float dt_ov_dx, dt_ov_dy, dt_ov_dz, dt;

  Solver(CkMigrateMessage *m) {}

  virtual void pup(PUP::er &p) {
    p | dt_ov_dx;
    p | dt_ov_dy;
    p | dt_ov_dz;
    p | dt;
  }

  Solver() {}
  virtual ~Solver(){};


  Solver(float box_width, float time_delta) {
    dt = time_delta;
    dt_ov_dx = time_delta / box_width;
    dt_ov_dy = time_delta / box_width;
    dt_ov_dz = time_delta / box_width;
  }

  inline Field poisson2D(float charge_density) {
    vec3 e_field(charge_density / EPSILON_FREE, 0.0, 0.0);
    vec3 b_field(0.0, 0.0, 0.0);
    return Field(e_field, b_field);
  }

  inline void maxwell_solver2D(vec3 J, Field *field, vec3 B_inext,
                              vec3 B_jnext, vec3 E_iprev, vec3 E_jprev) {
      vec3 B_cur = field->magnetic_field;
      vec3 E_cur = field->electric_field;

      field->electric_field.x +=
          -dt * J.x + dt_ov_dy * (B_jnext.z - B_cur.z);
      field->electric_field.y +=
          -dt * J.y - dt_ov_dx * (B_inext.z - B_cur.z);
      field->electric_field.z += -dt * J.z +
                                     dt_ov_dx * (B_inext.y - B_cur.y) -
                                     dt_ov_dy * (B_jnext.x - B_cur.x);

      field->magnetic_field.x -= dt_ov_dy * (E_cur.z - E_jprev.z);
      field->magnetic_field.y += dt_ov_dx * (E_cur.z - E_iprev.z);
      field->magnetic_field.z +=
          dt_ov_dy * (E_cur.x - E_jprev.x) - dt_ov_dx * (E_cur.y - E_iprev.y);
  }
};

#endif // SOLVER_H_
