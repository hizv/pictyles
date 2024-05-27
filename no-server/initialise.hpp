#ifndef INITIALISE_H_
#define INITIALISE_H_

#include <cmath>

// Linearly distributed particle count
int linear_init2D(int initial_particle_count, int box_count, int x,
                  double alpha, double beta) {
  double step = 1.0 / (box_count - 1);

  // The linear function is f(x) = -alpha * x + beta , x in [0,1]
  double total_weight =
      beta * box_count - alpha * 0.5 * step * box_count * (box_count - 1);
  double current_weight = beta - (alpha * step * x);

  return initial_particle_count * (current_weight / total_weight);
}

// Geometrically distributed particle count
int geometric_init2D(int initial_particle_count, int box_count, int x,
                     double rho) {

  // Each cell in the i-th column of cells contains p(i) = A * rho^i particles
  double A = initial_particle_count *
             ((1.0 - rho) / (1.0 - std::pow(rho, box_count))) /
             (double)box_count;
  return A * pow(rho, x);
}

// Geometrically distributed particle count
int sinusoidal_init2D(int initial_particle_count, int box_count, int x) {
  double step = M_PI / box_count;
  return 2.0 * std::cos(x * step) * std::cos(x * step) *
         initial_particle_count / (box_count * box_count);
}

#endif // INITIALISE_H_
