#ifndef INTERPOLATE_H_
#define INTERPOLATE_H_

// Function to deposit scalar densities (charge etc) on 2D grid corners
inline void deposit_density2D(double grid_size, const vec3& position,
                        double particle_charge,
                        double* charge_densities) {
  // Calculate the grid indices for the given position
  int ix = static_cast<int>(floor(position.x / grid_size));
  int iy = static_cast<int>(floor(position.y / grid_size));

  // Calculate the fractional distances within the grid cell
  double dx = position.x / grid_size - ix;
  double dy = position.y / grid_size - iy;

  // Bilinear interpolation coefficients
  double w00 = (1.0 - dx) * (1.0 - dy);
  double w01 = (1.0 - dx) * dy;
  double w10 = dx * (1.0 - dy);
  double w11 = dx * dy;

  // Deposit charge densities on grid corners
  charge_densities[0] += w00 * particle_charge;
  charge_densities[1] += w01 * particle_charge;
  charge_densities[2] += w10 * particle_charge;
  charge_densities[3] += w11 * particle_charge;

}


// Function to deposit charge and current densities on 3D grid corners
inline void deposit_density3D(double grid_size, const vec3& position,
                        double particle_charge,
                        double* charge_densities) {
  // Calculate the grid indices for the given position
  int ix = static_cast<int>(floor(position.x / grid_size));
  int iy = static_cast<int>(floor(position.y / grid_size));
  int iz = static_cast<int>(floor(position.z / grid_size));

  // Calculate the fractional distances within the grid cell
  double dx = position.x / grid_size - ix;
  double dy = position.y / grid_size - iy;
  double dz = position.z / grid_size - iz;

  // Bilinear interpolation coefficients
  double w00 = (1.0 - dx) * (1.0 - dy);
  double w01 = (1.0 - dx) * dy;
  double w10 = dx * (1.0 - dy);
  double w11 = dx * dy;

  // Deposit charge densities on grid corners
  charge_densities[0] += w00 * particle_charge;
  charge_densities[1] += w01 * particle_charge;
  charge_densities[2] += w10 * particle_charge;
  charge_densities[3] += w11 * particle_charge;

}


// Function to deposit vector densities (current etc) on 2D grid corners
inline void deposit_density2Dvec(double grid_size, const vec3& position,
                        vec3 current_density,
                        vec3* current_densities) {
  // Calculate the grid indices for the given position
  int ix = static_cast<int>(floor(position.x / grid_size));
  int iy = static_cast<int>(floor(position.y / grid_size));

  // Calculate the fractional distances within the grid cell
  double dx = position.x / grid_size - ix;
  double dy = position.y / grid_size - iy;

  // Bilinear interpolation coefficients
  double w00 = (1.0 - dx) * (1.0 - dy);
  double w01 = (1.0 - dx) * dy;
  double w10 = dx * (1.0 - dy);
  double w11 = dx * dy;

  // Deposit charge densities on grid corners
  current_densities[0] += current_density * w00;
  current_densities[1] += current_density * w01;
  current_densities[2] += current_density * w10;
  current_densities[3] += current_density * w11;

}
// Function to calculate field on particle, interpolating fields from all corners
inline Field interpolated_field2D(double grid_size, const vec3& position,
                               const Field* field) {
  // Calculate the grid indices for the given position
  int ix = static_cast<int>(floor(position.x / grid_size));
  int iy = static_cast<int>(floor(position.y / grid_size));

  // Calculate the fractional distances within the grid cell
  double dx = position.x / grid_size - ix;
  double dy = position.y / grid_size - iy;

  // Bilinear interpolation coefficients
  double w00 = (1.0 - dx) * (1.0 - dy);
  double w01 = (1.0 - dx) * dy;
  double w10 = dx * (1.0 - dy);
  double w11 = dx * dy;

  Field field_at_point;

  // Interpolate electric field
  vec3 efield_at_point = field[0].electric_field * w00
                        + field[1].electric_field * w01
                        + field[2].electric_field * w10
                        + field[3].electric_field * w11;
  field_at_point.electric_field = efield_at_point;

  if (std::isnan(field_at_point.electric_field.x)) {
      CkPrintf("[interpolated_field2D] E.x is NaN!\n");
  }


  // Interpolate magnetic field
  vec3 bfield_at_point = field[0].magnetic_field * w00
                        + field[1].magnetic_field * w01
                        + field[2].magnetic_field * w10
                        + field[3].magnetic_field * w11;
  field_at_point.magnetic_field = bfield_at_point;

  if (std::isnan(field_at_point.magnetic_field.x)) {
      CkPrintf("[interpolated_field2D] B.x is NaN!\n");
  }


  return field_at_point;
}

// calculate force using electric and magnetic fields
inline vec3 field_to_force(const Field& field, const Particle& p) {
  vec3 force;
  force = field.electric_field * p.charge; // F_E = qE


  // F_B = q (v cross B)
  vec3 v = p.velocity;
  vec3 cross_product = { field.magnetic_field.y * v.z - field.magnetic_field.z * v.y,
                            field.magnetic_field.z * v.x - field.magnetic_field.x * v.z,
                            field.magnetic_field.x * v.y - field.magnetic_field.y * v.x };
  force += cross_product * p.charge;

  if (std::isnan(force.x)) {
      CkPrintf("[field_to_force] Force.x is NaN!, charge: %lf\n", p.charge * 1e12);
  }
  return force;
}
#endif // INTERPOLATE_H_
