#ifndef SOLVER_H_
#define SOLVER_H_

struct Solver {
  virtual void field_solver(Field*, const vec3*) = 0;
  virtual Field poisson(float) = 0;
};

struct Solver_Maxwell_Cart2D : public Solver {
  float dt_ov_dx, dt_ov_dy, dt_ov_dz, dt;

  Solver_Maxwell_Cart2D() {}

  virtual ~Solver_Maxwell_Cart2D()  {};

  Solver_Maxwell_Cart2D(float box_width, float time_delta) {
    dt = time_delta;
    dt_ov_dx = time_delta / box_width;
    dt_ov_dy = time_delta / box_width;
    dt_ov_dz = time_delta / box_width;
  }

  Field poisson(float charge_density) {
    vec3 e_field(charge_density / EPSILON_FREE, 0.0, 0.0);
    vec3 b_field(0.0, 0.0, 0.0);
    return Field(e_field, b_field);
  }

  void field_solver(Field *field, const vec3* context) {
      vec3 B_cur = field->magnetic_field;
      vec3 E_cur = field->electric_field;

      vec3 J = context[0], B_inext = context[1], B_jnext = context[2],
           E_iprev = context[3], E_jprev = context[4];
      field->electric_field.x += -dt * J.x + dt_ov_dy * (B_jnext.z - B_cur.z);
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

struct Solver_Maxwell_Cart3D : public Solver {
  float dt_ov_dx, dt_ov_dy, dt_ov_dz, dt;

  Solver_Maxwell_Cart3D() {}

  virtual ~Solver_Maxwell_Cart3D()  {};

  Solver_Maxwell_Cart3D(float box_width, float time_delta) {
    dt = time_delta;
    dt_ov_dx = time_delta / box_width;
    dt_ov_dy = time_delta / box_width;
    dt_ov_dz = time_delta / box_width;
  }

  Field poisson(float charge_density) {
    vec3 e_field(charge_density / EPSILON_FREE, 0.0, 0.0);
    vec3 b_field(0.0, 0.0, 0.0);
    return Field(e_field, b_field);
  }

  void field_solver(Field *field, const vec3* context) {
      vec3 B_cur = field->magnetic_field;
      vec3 E_cur = field->electric_field;

      vec3 J = context[0], B_inext = context[1], B_jnext = context[2],
           E_iprev = context[3], E_jprev = context[4];
      field->electric_field.x += -dt * J.x + dt_ov_dy * (B_jnext.z - B_cur.z);
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


struct SolverFactory {
  static Solver *create(PicParams& pp) {
    float box_width = pp.sim_box_length / (float)pp.box_count;
    switch(pp.geometry) {
      case CART2D:
        return new Solver_Maxwell_Cart2D(box_width, pp.time_delta);
      case CART3D:
        return new Solver_Maxwell_Cart3D(box_width, pp.time_delta);
    }
  }
};


#endif // SOLVER_H_
