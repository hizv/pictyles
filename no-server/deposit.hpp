#ifndef DEPOSIT_H_
#define DEPOSIT_H_

struct Depositor {
  // Function to deposit properties (charge etc) on 2D grid corners
  template <class T>
  void deposit_property(double, const vec3 &, T, T*) {}
};

struct Depositor2D : public Depositor {
  template <class T>
  void deposit_property(double grid_size, const vec3 &position,
                              T property, T *properties) {
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
    properties[0] += property * w00;
    properties[1] += property * w01;
    properties[2] += property * w10;
    properties[3] += property * w11;
  }
};

struct Depositor3D : public Depositor {
  template <class T>
  void deposit_property(double grid_size, const vec3 &position,
                              T property, T *properties) {
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

    // Deposit charge densities on grid corners
    for (int i = 0; i < 8; i++)
      properties[i] += property * w[i];
  }
};

struct DepositorFactory {
  static Depositor create(uint8_t geometry) {
    switch(geometry) {
      case CART2D:
        return Depositor2D();
      case CART3D:
        return Depositor3D();
    }
  }
};


#endif // DEPOSIT_H_
