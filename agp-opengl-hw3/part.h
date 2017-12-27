#ifndef _PART_HEADER
#define _PART_HEADER

#include "common.hpp"
#include "util.hpp"
//#include "math_functions.h"
#include "curand_kernel.h"
// #include "helper_math.h"
//#include "cuda.h"

#define P_NUM 100
#define TPB 256

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

//struct Particle
//{
//  float3 position;
//  float3 velocity;
//  bool p_type;   // true :silicate  false :iron
//};



//void launch_kernel_initial_position_velocity(unsigned seed, struct Particle *particles);
//void launch_kernel_calculate_forces(void *devP, void *devA);
//void launch_kernel_update_pos(void *devP, void *devA);
//void launch_kernel_update_vel(void *devP, void *devA, void *devA_next);

void particle_init(unsigned seed, struct Particle *cpuP);
void particle_update(struct Particle *cpuP);
void cuda_release();


#endif

