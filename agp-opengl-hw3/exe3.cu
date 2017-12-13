#include <stdio.h>
#include <iostream>
#include <sys/time.h>
#include <iomanip>
using namespace std; 

const int N=10000;//NUM_PARTICLES
const int T=10000;//NUM_ITERATION
const int BLOCK_SIZE = 256;
const int dt=1;

struct Particle{
float px;
float py;
float pz;
float vx;
float vy;
float vz;};



__global__ void cal(Particle *Pd, float random, int N)
{
    int tid = blockIdx.x* blockDim.x +threadIdx.x;

    if (tid<N)
    {	
        Pd[tid].vx += random;
        Pd[tid].vy -= random;
        Pd[tid].vz = Pd[tid].vx - Pd[tid].vy;

        Pd[tid].px += Pd[tid].vx*dt;
        Pd[tid].py += Pd[tid].vy*dt;
        Pd[tid].pz += Pd[tid].vz*dt; 
    } 
}

void launchkernel(int argc, char *argv[])
{

    //int N=atoi(argv[1]);//NUM_PARTICLES
    //int T=atoi(argv[2]);//NUM_ITERATION
    //int BLOCK_SIZE=atoi(argv[3]);
    Particle *P = (Particle *)malloc(sizeof(Particle) * N);
    Particle *P1 = (Particle *)malloc(sizeof(Particle) * N);

    Particle *Pd;
    struct timeval t1,t2;
    double time;
    bool success = true;
    int gridsize = (N + BLOCK_SIZE-1) / BLOCK_SIZE ;
    dim3 dimBlock(BLOCK_SIZE, 1 ,1);
    dim3 dimGrid(gridsize, 1 ,1);

    
    printf("Computing Particles on the GPU...");

    cudaMalloc((void**)&Pd, sizeof(Particle) * N); 
    cudaMemcpy(Pd, P,sizeof(Particle) * N, cudaMemcpyHostToDevice); 

    // Make sure we use the same seed
    srand(5);
    
   
        float random = rand() / float(RAND_MAX);
        
        cal<<<dimGrid, dimBlock>>>(Pd, random, N);
    

    cudaMemcpy( P1, Pd, sizeof(Particle) * N, cudaMemcpyDeviceToHost );




    return 0;
 }
