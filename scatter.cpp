#include <iostream>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "scatter.h"

using namespace std;

static gsl_rng* rng= 0;
static double npow;
static size_t nwarn= 0;

static void random_direction(double e[]);
static int scat_pair(double* const v1, double* const v2, const double dtau,
		     ScatterEvent* const ev);


void scatter_init(gsl_rng* const random_generator, const double n_power_law)
{
  rng= random_generator;
  npow= n_power_law;

  printf("# npow = %lf\n", npow);
}

size_t scatter_n_warning()
{
  return nwarn;
}

void scatter(Particle* const p, const size_t n, const double dtau, ScatterEvent* const ev)
{
  // Scatter one pair of particles
  // dtau = dt/t_r; where t_r= 1/(rho sigma(v0) v0)
  // Returns:
  //    ev->dtau; time passed
  assert(rng);
  assert(n > 1);
  
  int ntrial= 0;
  int ret= 0;
  
  while(1) {
    const size_t i= (size_t) (n*gsl_rng_uniform(rng)); 
    const size_t j= (size_t) (n*gsl_rng_uniform(rng));

    if(i == n || j == n || i == j) {
      continue;
    }
    
    ret= scat_pair(p[i].v, p[j].v, dtau, ev);

    ntrial++;

    if(ret == 1) {
      if(ev) {
	ev->i= i;
	ev->j= j;
	ev->dtau= ntrial*dtau/n;
      }
      return;
    }
  };
}

int scat_pair(double* const v1, double* const v2, const double dtau, ScatterEvent* const ev)
{
  // Returns:
  //    1 if v1 and v2 scatter, 0 otherwise

  const double rv3[3]= {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]};
  const double rv= norm(rv3);

  // Probabily that the pair i-j scatter during time interval dtau
  double prob= pow(rv, npow + 1.0)*dtau;
  if(prob > 1.0) {
    fprintf(stderr, "Warning: dtau too large; prob= %f\n", prob);
    nwarn++;
  }
  
  if(gsl_rng_uniform(rng) > prob)
    return 0;
  
  const double vcm[]= {0.5*(v1[0] + v2[0]),
		       0.5*(v1[1] + v2[1]),
		       0.5*(v1[2] + v2[2])};
  
  double e[3];
  random_direction(e);

  if(ev) {
    for(int k=0; k<3; ++k) {
      ev->v1_before[k] = v1[k];
      ev->v2_before[k] = v2[k];
      ev->rv= rv;
    }
  }
  v1[0]= vcm[0] + 0.5*rv*e[0];
  v1[1]= vcm[1] + 0.5*rv*e[1];
  v1[2]= vcm[2] + 0.5*rv*e[2];

  v2[0]= vcm[0] - 0.5*rv*e[0];
  v2[1]= vcm[1] - 0.5*rv*e[0];
  v2[2]= vcm[2] - 0.5*rv*e[0];

  if(ev) {
    for(int k=0; k<3; ++k) {
      ev->v1_after[k] = v1[k];
      ev->v2_after[k] = v2[k];
    }
  }

  return 1;
}


  
void random_direction(double e[])
{
  double y1, y2, r2;

  do{
    y1= 1.0f-2.0*gsl_rng_uniform(rng);
    y2= 1.0f-2.0*gsl_rng_uniform(rng);
    r2= y1*y1 + y2*y2;
  } while(r2 > 1.0);

  double sq1r= sqrt(1.0 - r2);
  e[0]= 2*y1*sq1r;
  e[1]= 2*y2*sq1r;
  e[2]= 1.0 - 2.0*r2;
}
