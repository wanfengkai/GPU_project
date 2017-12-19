#ifndef _PART_HEADER
#define _PART_HEADER

#include "common.hpp"
#include "util.hpp"

# define P_NUM 100

struct Particle
{
	glm::vec3  position;
	glm::vec3  velocity;
	bool p_type;   // True :silicate  False :iron
};



void launchkernel(Particle *cpuP,Particle *cpuA);

#endif

