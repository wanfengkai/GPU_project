// #include "math_functions.h"
#include "curand_kernel.h"
// #include "helper_math.h"
//#include "cuda.h"
#include <stdio.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>
#include "part.h"






// particles' position & velocity
__device__ void generate_uniform_random_number(curandState state, double* rho1, double* rho2, double* rho3 ){
    // 0 -1 range
    *rho1= curand_uniform(&state);
    *rho2= curand_uniform(&state);
    *rho3= curand_uniform(&state);
}


__global__ void initial_position_velocity(unsigned seed, struct Particle *particles)
{

    int i= blockIdx.x*blockDim.x + threadIdx.x;

	if(i<P_NUM){

    double rho1, rho2, rho3;
    double miu;

    curandState state;
    curand_init(seed, i, 0, &state);

    bool planet;  // true: Earth ; false: Moon

    if (curand_uniform(&state)>0.5)
        particles[i].p_type = true  ;  // silicate particle
    else
        particles[i].p_type = false ;  // iron particle

    if (curand_uniform(&state)>0.5)
        planet = true   ;  // Earth particle
    else
        planet = false ;  // Moon particle
    
    
    if (planet == true) {   
        if (particles[i].p_type == true) {
            // position initialization for outer shell, silicate particles
            generate_uniform_random_number(state, &rho1, &rho2, &rho3 );
            miu= 1- 2 * rho2;   
            particles[i].position.x = cbrt(pow(R1,3.0)+(pow(R2,3.0)-pow(R1,3.0))*rho1) * sqrt(1-pow(miu,2.0))*cos(2*PI*rho3)+ CENTER_MASS_X;
            particles[i].position.y = cbrt(pow(R1,3.0)+(pow(R2,3.0)-pow(R1,3.0))*rho1) * sqrt(1-pow(miu,2.0))*sin(2*PI*rho3);
            particles[i].position.z = cbrt(pow(R1,3.0)+(pow(R2,3.0)-pow(R1,3.0))*rho1) * miu + CENTER_MASS_Z;
        }
        else {
            // position initialization for inner core, iron particles   
            generate_uniform_random_number(state, &rho1, &rho2, &rho3 );
            miu= 1- 2 * rho2;
            particles[i].position.x = R * cbrt(rho1) * sqrt(1-pow(miu,2.0)) * cos(2*PI*rho3) + CENTER_MASS_X;
            particles[i].position.y = R * cbrt(rho1) * sqrt(1-pow(miu,2.0)) * sin(2*PI*rho3) ;
            particles[i].position.z = R * cbrt(rho1) * miu + CENTER_MASS_Z;
        }

        // velocity initialization
        particles[i].velocity.x = Init_V;
        particles[i].velocity.y = 0;
        particles[i].velocity.z = 0;
    }
    else {
        if (particles[i].p_type == true) {
            // position initialization for outer shell, silicate particles
            generate_uniform_random_number(state, &rho1, &rho2, &rho3 );
            miu= 1- 2 * rho2;   
            particles[i].position.x = cbrt(pow(R1,3.0)+(pow(R2,3.0)-pow(R1,3.0))*rho1) * sqrt(1-pow(miu,2.0))*cos(2*PI*rho3)- CENTER_MASS_X;
            particles[i].position.y = cbrt(pow(R1,3.0)+(pow(R2,3.0)-pow(R1,3.0))*rho1) * sqrt(1-pow(miu,2.0))*sin(2*PI*rho3);
            particles[i].position.z = cbrt(pow(R1,3.0)+(pow(R2,3.0)-pow(R1,3.0))*rho1) * miu;
        }
        else {
            // position initialization for inner core, iron particles   
            generate_uniform_random_number(state, &rho1, &rho2, &rho3 );
            miu= 1- 2 * rho2;
            particles[i].position.x = R * cbrt(rho1) * sqrt(1-pow(miu,2.0)) * cos(2*PI*rho3) - CENTER_MASS_X;
            particles[i].position.y = R * cbrt(rho1) * sqrt(1-pow(miu,2.0)) * sin(2*PI*rho3) ;
            particles[i].position.z = R * cbrt(rho1) * miu - CENTER_MASS_Z;
        }

        // velocity initialization
        particles[i].velocity.x = -Init_V;
        particles[i].velocity.y = 0;
        particles[i].velocity.z = 0;
    }
}

    
}

// interaction forces

__device__ double interactionForce(double radius, double sqr_r, int force_type, bool if_sep_dec){
// radius is the parapmeter 'r' in Table 2
    double M1, M2, K1, K2, KRP1, KRP2;
    double accel;
    if (force_type == 0) {  // accel for silicate particle, silicate-iron 
        M1 = Mfe;
		M2 = Msi;
        K1 = Ksi;
        K2 = Kfe;
        KRP1 = KRPsi;
        KRP2 = KRPfe;
    }
    else if (force_type == 1) { // accel for iron particle, iron-silicate
        M1 = Msi;
		M2 = Mfe;
        K1 = Ksi;
        K2 = Kfe;
        KRP1 = KRPsi;
        KRP2 = KRPfe;
    }
    else if (force_type == 2) { // accel for silicate particle, silicate-silicate
        M1 = Msi;
		M2 = Msi;
        K1 = Ksi;
        K2 = Ksi;
        KRP1 = KRPsi;
        KRP2 = KRPsi;
    }
    else if (force_type == 3) { // accel for iron particle, iron-iron
        M1 = Mfe;
		M2 = Mfe;
        K1 = Kfe;
        K2 = Kfe;
        KRP1 = KRPfe;
        KRP2 = KRPfe;
    }
    else {
      //  fprintf(stderr,"invalid force type");
    }


    if (radius >= D){
        accel = G * M1 / sqr_r;
    }  
    else if (radius >= D-D*SDPsi){
        accel = G * M1 / sqr_r - 0.5*(K1+K2)*(D*D-sqr_r)/M2;
    }

    else if (radius >= D-D*SDPfe && if_sep_dec==true){
        accel = G * M1 / sqr_r - 0.5*(K1+K2)*(D*D-sqr_r)/M2 ;
    }
    else if (radius >= D-D*SDPfe && if_sep_dec==false){
        accel = G * M1 / sqr_r - 0.5*(K1*KRP1+K2)*(D*D-sqr_r)/M2 ;
    }
    else if (radius >= Epsilon  && if_sep_dec==true){
        accel = G * M1 / sqr_r - 0.5*(K1+K2)*(D*D-sqr_r)/M2 ;
    }
    else if (radius >= Epsilon  && if_sep_dec==false){
        accel = G * M1 / sqr_r - 0.5*(K1*KRP1+K2*KRP2)*(D*D-sqr_r)/M2 ;
    }
    else if(radius < Epsilon && radius > 0.11){
//        radius=Epsilon;
		sqr_r = Epsilon*Epsilon;
        accel = G * M1 / sqr_r - 0.5*(K1+K2)*(D*D-sqr_r)/M2 ;
    }
	else {
		accel = 0;
	}

    return accel;
}


__device__ struct vecfloat3 bodyBodyInteraction(struct Particle pi, struct Particle pj, struct vecfloat3 acc)
{
  struct vecfloat3 r;
  struct vecfloat3 v;
  struct vecfloat3 r_next;
  double accel;
  // r_ij  [3 FLOPS]
  r.x = pj.position.x - pi.position.x;
  r.y = pj.position.y - pi.position.y;
  r.z = pj.position.z - pi.position.z;

  v.x = pj.velocity.x - pi.velocity.x;
  v.y = pj.velocity.y - pi.velocity.y;
  v.z = pj.velocity.z - pi.velocity.z;

  r_next.x = r.x + v.x * Time_step;
  r_next.y = r.y + v.y * Time_step;
  r_next.z = r.z + v.z * Time_step;

  int force_type;
  bool if_sep_dec;
  float sqr_r = r.x * r.x + r.y * r.y + r.z * r.z;
  float sqr_r_next = r_next.x * r_next.x + r_next.y * r_next.y + r_next.z * r_next.z;
  if (sqr_r_next < sqr_r) {
    if_sep_dec = 1;
  }
  sqr_r += softeningSquared;
  float radius = sqrtf(sqr_r);

  if ( pi.p_type && !pj.p_type ) {
    force_type = 0;
  } 
  else if ( !pi.p_type && pj.p_type ) {
    force_type = 1;
  }
  else if ( pi.p_type && pj.p_type ) {
    force_type = 2;
  }
  else if ( !pi.p_type && !pj.p_type ) {
    force_type = 3;
  }
  else {
    // fprintf(stderr,"invalid force type!");
  }
  accel = interactionForce(radius, sqr_r, force_type, if_sep_dec);
  
  acc.x += r.x / radius * accel;
  acc.y += r.y / radius * accel;
  acc.z += r.z / radius * accel;

  return acc;
}

__device__ struct vecfloat3 tile_calculation(struct Particle myParticle, struct vecfloat3 acc)
{
  int i;
  extern __shared__ struct Particle shParticles[];
  for (i = 0; i < blockDim.x; i++) {
	if (isnan(myParticle.position.x) || isnan(shParticles[i].position.x)) 
		printf("i is %d\n", i);
	else
    	acc = bodyBodyInteraction(myParticle, shParticles[i], acc);
  }
  return acc;
}

__global__ void calculate_forces(struct Particle *devP, struct vecfloat3 *devA)
{
  extern __shared__ struct Particle shParticles[]; 
  struct Particle myParticle;
  int i, tile;
  struct vecfloat3 acc; //= {0.0f, 0.0f, 0.0f}
  acc.x = 0.0f;
  acc.y = 0.0f;
  acc.z = 0.0f;
  int gtid = blockIdx.x * blockDim.x + threadIdx.x;
  if(gtid < P_NUM){
    myParticle = devP[gtid];
	if (isnan(myParticle.position.x)) 
		printf("gtid is %d\n", gtid);
  	for (i = 0, tile = 0; i < P_NUM; i += TPB, tile++) {
    	int idx = tile * blockDim.x + threadIdx.x;
    	shParticles[threadIdx.x] = devP[idx];
    	__syncthreads();
    	acc = tile_calculation(myParticle, acc);
    	__syncthreads();
  	}
	devA[gtid] = acc;

  } 

}

__global__ void update_pos(struct Particle *devP, struct vecfloat3 *devA)
{

  int i= blockIdx.x*blockDim.x + threadIdx.x;

  if(i<P_NUM){
  	devP[i].position.x = devP[i].position.x + devP[i].velocity.x * Time_step + devA[i].x / 2.0f * pow(Time_step, 2.0); 
  	devP[i].position.y = devP[i].position.y + devP[i].velocity.y * Time_step + devA[i].y / 2.0f * pow(Time_step, 2.0); 
  	devP[i].position.z = devP[i].position.z + devP[i].velocity.z * Time_step + devA[i].z / 2.0f * pow(Time_step, 2.0); 
  }
	
}

__global__ void update_vel(struct Particle *devP, struct vecfloat3 *devA, struct vecfloat3 *devA_next)
{

  int i= blockIdx.x*blockDim.x + threadIdx.x;

  if(i<P_NUM){
  	devP[i].velocity.x = devP[i].velocity.x + (devA_next[i].x + devA[i].x) / 2.0f * Time_step; 
  	devP[i].velocity.y = devP[i].velocity.y + (devA_next[i].y + devA[i].y) / 2.0f * Time_step; 
  	devP[i].velocity.z = devP[i].velocity.z + (devA_next[i].z + devA[i].z) / 2.0f * Time_step; 
  }


}

__global__ void print_devP(struct Particle *devP, struct vecfloat3 *devA)
{

  int i= blockIdx.x*blockDim.x + threadIdx.x;
  if(i<P_NUM)
	printf("GPU particle position is: %f, %f, %f\n", devP[i].position.x, devP[i].position.y, devP[i].position.z);
}


void particle_init(unsigned seed, struct Particle *cpuP)
{
	struct Particle *devP=0;

	cudaMalloc(&devP, P_NUM*sizeof(struct Particle));
	
    initial_position_velocity<<<(P_NUM+TPB-1)/TPB, TPB>>>(seed, devP);

    cudaDeviceSynchronize();

    cudaMemcpy(cpuP, devP, P_NUM*sizeof(struct Particle),cudaMemcpyDeviceToHost);
	cudaFree(devP);
}


void particle_update(struct Particle *cpuP)
{
	struct Particle *devP=0;
	struct vecfloat3 *devA=0;
	struct vecfloat3 *devA_next=0;
	checkCudaErrors(cudaMalloc(&devP, P_NUM*sizeof(struct Particle)));
	checkCudaErrors(cudaMalloc(&devA, P_NUM*sizeof(struct vecfloat3)));
	checkCudaErrors(cudaMalloc(&devA_next, P_NUM*sizeof(struct vecfloat3)));
 
	checkCudaErrors(cudaMemcpy(devP, cpuP, P_NUM*sizeof(struct Particle),cudaMemcpyHostToDevice));
	cudaDeviceSynchronize();

    // update position
    calculate_forces<<<(P_NUM+TPB-1)/TPB, TPB, TPB*sizeof(struct Particle)>>>(devP, devA);
    getLastCudaError("Kernel execution failed");  	// check if kernel execution generated and error
    cudaDeviceSynchronize();

    update_pos<<<(P_NUM+TPB-1)/TPB, TPB>>>(devP, devA);
	getLastCudaError("Kernel execution failed");  	// check if kernel execution generated and error
    cudaDeviceSynchronize();
	
    // update velocity
    calculate_forces<<<(P_NUM+TPB-1)/TPB, TPB, TPB*sizeof(struct Particle)>>>(devP, devA_next);
	getLastCudaError("Kernel execution failed");  	// check if kernel execution generated and error
    cudaDeviceSynchronize();

    update_vel<<<(P_NUM+TPB-1)/TPB, TPB>>>(devP, devA, devA_next);
	getLastCudaError("Kernel execution failed");  	// check if kernel execution generated and error
    cudaDeviceSynchronize();
	
	checkCudaErrors(cudaMemcpy(cpuP, devP, P_NUM*sizeof(struct Particle),cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaFree(devP));
    checkCudaErrors(cudaFree(devA));
    checkCudaErrors(cudaFree(devA_next));

}

//void cuda_release()
//{
//    cudaFree(devP);
//    cudaFree(devA);
//    cudaFree(devA_next);
//}




