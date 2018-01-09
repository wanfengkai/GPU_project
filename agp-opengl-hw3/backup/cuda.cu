#include <stdio.h>
#include "part.h"

const int N=100;//NUM_PARTICLES
const int BLOCK_SIZE = 256;
const int dt=1;
int gridsize = (N + BLOCK_SIZE-1) / BLOCK_SIZE ;
dim3 dimBlock(BLOCK_SIZE, 1 ,1);
dim3 dimGrid(gridsize, 1 ,1);

__global__ void update(Particle *Pd)
{
    int tid = blockIdx.x* blockDim.x +threadIdx.x;

    if (tid<N)
    {	
	Pd[tid].v.x =  0.0001f ;
        Pd[tid].v.y =  0.0001f ;
        Pd[tid].v.z =  0.0f ;
        Pd[tid].p.x += Pd[tid].v.x * dt;
        Pd[tid].p.y += Pd[tid].v.y * dt;
        Pd[tid].p.z += Pd[tid].v.z * dt; 
    } 
}

void launchkernel(Particle *P)
{

    Particle *Pd;

    cudaMalloc((void**)&Pd, sizeof(Particle) * N); 
    cudaMemcpy(Pd, P,sizeof(Particle) * N, cudaMemcpyHostToDevice); 

    update<<<dimGrid, dimBlock>>>(Pd);
 
    cudaMemcpy( P, Pd, sizeof(Particle) * N, cudaMemcpyDeviceToHost );
    cudaFree(Pd);

}
