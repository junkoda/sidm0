#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <boost/program_options.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


using namespace std;
using namespace boost::program_options;

struct Particle {
  double v[3];
};

struct Options {
  size_t n;
  double vmax;
};

static gsl_rng* rng;

static Options* parse_options(int argc, char* argv[]);
static void rng_init();
static Particle* set_initial_gaussian(const size_t n);
static size_t remove_high_speed_particles(Particle* const p, const size_t n,
					  const double vmax);

static inline double norm(const double v[])
{
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static inline double norm2(const double v[])
{
  return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

//
// main
//
int main(int argc, char* argv[])
{
  rng_init();
  Options* opt = parse_options(argc, argv);
  if(opt == 0)
    return 0;
  
  size_t n= opt->n;
  Particle* const p= set_initial_gaussian(n);

  n= remove_high_speed_particles(p, n, opt->vmax);
  cerr << n << " particles remaining\n";


  
  return 0;
}

Options* parse_options(int argc, char* argv[])
{
  Options* o= new Options();
  
  options_description opt("sidm0 <n>");
  opt.add_options()
    ("help,h", "display this help")
    ("n", value<size_t>()->default_value(10000), "initial number of particles")
    ("vmax", value<double>()->default_value(2.0), "cutoff speed")
    ;

  positional_options_description p;
  p.add("n", -1);

  variables_map vm;
  store(command_line_parser(argc, argv).options(opt).positional(p).run(), vm);
  notify(vm);

  if(vm.count("help") || ! vm.count("n")) {
    cout << opt;
    return 0;
  }
  
  o->n= vm["n"].as<size_t>();
  o->vmax= vm["vmax"].as<double>();

  return o;
}

void rng_init()
{
  unsigned int seed= (unsigned int) time(NULL);
  rng= gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rng, seed);
}  

Particle* set_initial_gaussian(const size_t n)
{
  Particle* const p= (Particle*) malloc(sizeof(Particle)*n);

  for(int i=0; i<n; ++i) {
    p[i].v[0]= gsl_ran_gaussian(rng, 1.0);
    p[i].v[1]= gsl_ran_gaussian(rng, 1.0);
    p[i].v[2]= gsl_ran_gaussian(rng, 1.0);
  }

  return p;
}

size_t remove_high_speed_particles(Particle* const p, const size_t n,
				   const double vmax)
{
  const double vmax2= vmax*vmax;
  
  size_t j=0;
  for(size_t i=0; i<n; ++i) {
    if(norm2(p[i].v) < vmax2)
      p[j++]= p[i];
  }

  return j;
}
