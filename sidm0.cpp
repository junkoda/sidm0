//
// Zero dimensional SIDM toy model
//

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <boost/program_options.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "particle.h"
#include "scatter.h"

using namespace std;
using namespace boost::program_options;

static double energy_lost= 0.0;

static gsl_rng* rng;

static void rng_init();
static Particle* set_initial_gaussian(const size_t n);
static size_t remove_high_speed_particles_all(Particle* const p, const size_t n,
					      const double vmax2);

static size_t remove_high_speed_particles_pair(Particle* const p, const size_t n,
					  ScatterEvent const * const ev,
					  const double vmax2);
void test_nscat(Particle* const p, const size_t n, const double dtau,
		gsl_rng* const rng,
		const int nscat_max);
static void compute_rms(Particle const * const p, const size_t n);
static double compute_energy(Particle const * const p, const size_t n);

static inline double energy(double v[])
{
  return 0.5*(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

//
// main
//
int main(int argc, char* argv[])
{
  //
  // Read command-line options
  //
  variables_map vm;
  {
    options_description opt("sidm0 <n>");
    opt.add_options()
      ("help,h", "display this help")
      ("n", value<double>(), "VHS power-law index")
      ("n-particles", value<size_t>()->default_value(10000), "initial number of particles")
      ("vmax", value<double>()->default_value(2.0), "cutoff speed")
      ("test", "test number of scattering")
      ("final-ratio", value<double>()->default_value(0.8), "ratio of remaing particle")
      
      ;
    
    positional_options_description p;
    p.add("n", -1);
    
    store(command_line_parser(argc, argv).options(opt).positional(p).run(), vm);
    notify(vm);
    
    if(vm.count("help") || ! vm.count("n")) {
      cout << opt;
      return 0;
    }
  }
  
  size_t n= vm["n-particles"].as<size_t>();
  const double npow= vm["n"].as<double>();
  const double vmax= vm["vmax"].as<double>();
  const double vmax2= vmax*vmax;

  //
  // Initialise
  //
  unsigned int seed= (unsigned int) time(NULL);
  rng= gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rng, seed);

  scatter_init(rng, npow);

  Particle* const p= set_initial_gaussian(n);
  cerr << "initial energy " << compute_energy(p, n) << endl;

  double dtau= 0.01;
  if(npow >= 0) {
    dtau= 1.0/pow(2.0*vmax, npow + 1.0);
    cerr << "dtau = " << dtau << endl;
  }
  if(npow == -1.0) {
    dtau= 1.0;
  }
  if(npow == -3.0) {
    dtau= 0.001;
  }
      
  
  if(vm.count("test")) {
    compute_rms(p, n);

    test_nscat(p, n, dtau, rng, n);

    size_t nwarn= scatter_n_warning();
    fprintf(stderr, "warning %lu %e\n", nwarn, (double)nwarn/n);
    compute_rms(p, n);
    return 0;
  }
  
  n= remove_high_speed_particles_all(p, n, vmax2);
  cerr << n << " initial n-particles\n";
  size_t nmin= n*vm["final-ratio"].as<double>();
  const double energy_initial= compute_energy(p, n);
  
  //
  // main loop
  //
  double tau= 0.0;
  const double n0= n;
  printf("%e %e 1.0\n", tau, n/n0);
  
  ScatterEvent ev;
  double energy_remaining= energy_initial;
  cerr << n << " " << nmin << endl;
  
  while(n >= nmin && energy_remaining > 0.2*energy_initial) {
    scatter(p, n, dtau, &ev);
    tau += ev.dtau;
    size_t n_before= n;
    n= remove_high_speed_particles_pair(p, n, &ev, vmax2);
    energy_remaining= compute_energy(p, n);
    
    if(n < n_before)
      printf("%e %e %e\n",
	     tau,
	     n/n0,
	     (energy_remaining/n)/(energy_initial/n0)
	     );
    // Column 1: tau = t/tr
    // Column 2: n/n_initial
    // Column 3: temperature/temperature_initial
  }

  size_t nwarn= scatter_n_warning();
  fprintf(stderr, "warning %lu %e\n", nwarn, (double)nwarn/n);
  fprintf(stderr, "energy conservation %e/%e = %.4f\n",
	  energy_initial, energy_remaining + energy_lost,
	  (energy_remaining + energy_lost)/energy_initial);

  
  return 0;
}

Particle* set_initial_gaussian(const size_t n)
{
  // Allocate n particles with Maxwellian (Gaussian) velocity distribution
  // 1-dimensional velocity rms is 1.0
  Particle* const p= (Particle*) malloc(sizeof(Particle)*n);

  for(int i=0; i<n; ++i) {
    p[i].v[0]= gsl_ran_gaussian(rng, 1.0);
    p[i].v[1]= gsl_ran_gaussian(rng, 1.0);
    p[i].v[2]= gsl_ran_gaussian(rng, 1.0);
  }

  return p;
}

size_t remove_high_speed_particles_all(Particle* const p, const size_t n,
				   const double vmax2)
{
  // Remove all particles with speed^2 larger than vmax2
  
  size_t j=0;
  for(size_t i=0; i<n; ++i) {
    if(norm2(p[i].v) < vmax2)
      p[j++]= p[i];
  }

  return j;
}

static size_t remove_particle(Particle* const p, const size_t n, const size_t i)
{
  // Remove particle with index i
  // Shift all particles after i by one
  // Returns:
  //    new number of particles (n-1)
  for(size_t j=i+1; j<n; ++j)
    p[j-1]= p[j];

  return n-1;
}

size_t remove_high_speed_particles_pair(Particle* const p, size_t n,
					ScatterEvent const * const ev,
					const double vmax2)
{
  // Check a pair of particles ev->i or ev->j and
  // remove if they have speed^2 larger than vmax2
  
  if(norm2(p[ev->i].v) >= vmax2) {
    energy_lost += energy(p[ev->i].v);
    n= remove_particle(p, n, ev->i);
  }

  if(norm2(p[ev->j].v) >= vmax2) {
    energy_lost += energy(p[ev->j].v);
    n= remove_particle(p, n, ev->j);
  }

  return n;
}


void test_nscat(Particle* const p, const size_t n, const double dtau,
		gsl_rng* const rng,
		const int nscat_max)
{
  int nscat= 0;
  ScatterEvent ev;
  double tau= 0.0;
  
  for(nscat=0; nscat < nscat_max; ++nscat) {
    scatter(p, n, dtau, &ev);
    tau += ev.dtau;
    printf("%e %e %e\n", tau, (double)nscat/n, ev.rv);
  }
}

void compute_rms(Particle const * const p, const size_t n)
{
  double v2_sum[]= {0.0, 0.0, 0.0};
  for(int i=0; i<n; ++i) {
    v2_sum[0] += p[i].v[0]*p[i].v[0];
    v2_sum[1] += p[i].v[1]*p[i].v[1];
    v2_sum[2] += p[i].v[2]*p[i].v[2];
  }

  fprintf(stderr, "Velocity rms %.3f %.3f %.3f %.3e\n",
	  v2_sum[0]/n, v2_sum[1]/n, v2_sum[2]/n,
	  (v2_sum[0] + v2_sum[1] + v2_sum[2])/3.0/n);
}

double compute_energy(Particle const * const p, const size_t n)
{
  double v2_sum[]= {0.0, 0.0, 0.0};
  for(int i=0; i<n; ++i) {
    v2_sum[0] += p[i].v[0]*p[i].v[0];
    v2_sum[1] += p[i].v[1]*p[i].v[1];
    v2_sum[2] += p[i].v[2]*p[i].v[2];
  }

  return 0.5*(v2_sum[0] + v2_sum[1] + v2_sum[2]);
}
