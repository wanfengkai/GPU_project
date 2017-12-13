#include <stdio.h>
#include <iostream>
#include <sys/time.h>
#include <iomanip>
#include "part.h"

using namespace std; 

const int N=100;//NUM_PARTICLES
const int T=10000;//NUM_ITERATION
const int BLOCK_SIZE = 256;
const int dt=1;



__global__ void cal(Particle *Pd)
{
    int tid = blockIdx.x* blockDim.x +threadIdx.x;

    if (tid<N)
    {	
	Pd[i].v.x =  0.0001f* rand() / float(RAND_MAX) ;
        Pd[i].v.y =  0.0001f* rand() / float(RAND_MAX) ;
        Pd[i].v.z =  0.0f ;
        Pd[i].p.x += Pd[i].v.x * dt;
        Pd[i].p.y += Pd[i].v.y * dt;
        Pd[i].p.z += Pd[i].v.z * dt; 
    } 
}

void launchkernel(Particle *P)
{

    Particle *P = (Particle *)malloc(sizeof(Particle) * N);

    Particle *Pd;

    int gridsize = (N + BLOCK_SIZE-1) / BLOCK_SIZE ;
    dim3 dimBlock(BLOCK_SIZE, 1 ,1);
    dim3 dimGrid(gridsize, 1 ,1);

    
    cudaMalloc((void**)&Pd, sizeof(Particle) * N); 
    cudaMemcpy(Pd, P,sizeof(Particle) * N, cudaMemcpyHostToDevice); 

    cal<<<dimGrid, dimBlock>>>(Pd);
 

    cudaMemcpy( P, Pd, sizeof(Particle) * N, cudaMemcpyDeviceToHost );


    return 0;
 }
