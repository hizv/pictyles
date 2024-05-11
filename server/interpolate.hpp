#ifndef INTERPOLATE_H_
#define INTERPOLATE_H_

struct Interpolator {
  // Function to calculate field on particle, interpolating fields from all
  // corners
   virtual Field interpolated_field(double, const vec3 &, const Field *)  = 0;
};

struct Interpolator2D : public Interpolator {
   Field interpolated_field(double grid_size, const vec3 &position,
                           const Field *field) {
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
    vec3 efield_at_point =
        field[0].electric_field * w00 + field[1].electric_field * w01 +
        field[2].electric_field * w10 + field[3].electric_field * w11;
    field_at_point.electric_field = efield_at_point;

    // Interpolate magnetic field
    vec3 bfield_at_point =
        field[0].magnetic_field * w00 + field[1].magnetic_field * w01 +
        field[2].magnetic_field * w10 + field[3].magnetic_field * w11;
    field_at_point.magnetic_field = bfield_at_point;

    return field_at_point;
  }
};

struct Interpolator3D : public Interpolator {
   Field interpolated_field(double grid_size, const vec3 &position,
                           const Field *field) {
    // Calculate the grid indices for the given position
    int ix = static_cast<int>(floor(position.x / grid_size));
    int iy = static_cast<int>(floor(position.y / grid_size));
    int iz = static_cast<int>(floor(position.z / grid_size));

    // Calculate the fractional distances within the grid cell
    double dx = position.x / grid_size - ix;
    double dy = position.y / grid_size - iy;
    double dz = position.z / grid_size - dz;

    // Trilinear interpolation coefficients
    double w[8];
    w[0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
    w[1] = (1.0 - dx) * (1.0 - dy) * dz;
    w[2] = (1.0 - dx) * dy * (1.0 - dz);
    w[3] = (1.0 - dx) * dy * dz;
    w[4] = dx * (1.0 - dy) * (1.0 - dz);
    w[5] = dx * (1.0 - dy) * dz;
    w[6] = dx * dy * (1.0 - dz);
    w[7] = dx * dy * dz;

    Field field_at_point;

    vec3 efield_at_point = {0, 0, 0}, bfield_at_point = {0, 0, 0};

    for (int i = 0; i < 8; i++) {
      // Interpolate electric field
      efield_at_point += field[i].electric_field * w[i];

      // Interpolate magnetic field
      bfield_at_point += field[i].magnetic_field * w[i];
    }

    field_at_point.electric_field = efield_at_point;
    field_at_point.magnetic_field = bfield_at_point;

    return field_at_point;
  }
};

struct InterpolatorFactory {
  static Interpolator* create(uint8_t geometry) {
    switch (geometry) {
    case CART2D:
      return new Interpolator2D();
    case CART3D:
      return new Interpolator3D();
    }
  }
};

#endif // INTERPOLATE_H_
