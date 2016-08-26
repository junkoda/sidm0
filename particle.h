#ifndef PARTICLE_H
#define PARTICLE_H 1

struct Particle {
  double v[3]; // Dimensionless such that one-dimensional velocity rms = 1
};

static inline double norm(const double v[])
{
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static inline double norm2(const double v[])
{
  return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}


#endif
