#ifndef SCATTER_H
#define SCATTER_H 1

#include <gsl/gsl_rng.h>
#include "particle.h"

struct ScatterEvent {
  double dtau;
  size_t i, j;
  double v1_before[3], v2_before[3];
  double v1_after[3], v2_after[3];
  double rv;
};

void scatter_init(gsl_rng* const random_generator, const double n_power_law);

void scatter(Particle* const p, const size_t n, const double dtau,
	     ScatterEvent* const ev);

size_t scatter_n_warning();




#endif
