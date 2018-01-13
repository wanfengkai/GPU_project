#ifndef _PART_HEADER
#define _PART_HEADER

#include "common.hpp"
#include "util.hpp"
//#include "math_functions.h"
#include "curand_kernel.h"
// #include "helper_math.h"
//#include "cuda.h"


#define P_NUM 8192
#define TPB 256
// universal gravitational constant in km
#define G  6.67408E-20
#define Epsilon 47.0975
#define D 376.78
#define Mearth 5.97219E24
#define Msi Mearth*0.6/0.86/P_NUM //7.4161E21//19   
#define Mfe Mearth*0.4/0.14/P_NUM //1.9549E22//20
#define Ksi 2.9114E14//01//11 -11 +01 seems good as well
#define Kfe 5.8228E14//01//11 -11 +01 seems good as well
#define KRPsi 0.01
#define KRPfe 0.02
#define SDPsi 0.001
#define SDPfe 0.002
#define Time_step 5.8117
#define R 3185.5 
#define R1 3185.5 
#define R2 6371.0 
#define PI 3.14
#define Init_V -3.2416
#define CENTER_MASS_X 3185.5f
#define CENTER_MASS_Z 500.0f
//#define CENTER_MASS_X 2392.5
//#define CENTER_MASS_Z 9042.7
#define OMEGA 1 
#define softeningSquared 0.01f


struct vecfloat3
{
	float x;
	float y;
	float z;
};

struct Particle
{
  struct vecfloat3 position;
  struct vecfloat3 velocity;
  bool p_type;   // true :silicate  false :iron
};


void particle_init(unsigned seed, struct Particle *cpuP);
void particle_update(struct Particle *cpuP);
void cuda_release();


#endif

