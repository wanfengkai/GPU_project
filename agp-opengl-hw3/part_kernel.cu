// #include "math_functions.h"
#include "curand_kernel.h"
// #include "helper_math.h"
//#include "cuda.h"
#include <stdio.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>

// universal gravitational constant in km
#define G  6.67408E-20
#define Epsilon 47.0975
#define D 376.78
#define Msi 7.4161E19   
#define Mfe 1.9549E20
#define Ksi 2.9114E14
#define Kfe 5.8228E14
#define KRPsi 0.01
#define KRPfe 0.02
#define SDPsi 0.001
#define SDPfe 0.002
#define Time_step 5.8117
#define R 3185.5 
#define R1 3185.5 
#define R2 6371.0 
#define PI 3.14
#define Init_V 3.2416
#define CENTER_MASS_X 2392.5
#define CENTER_MASS_Z 9042.7
#define OMEGA 1 
#define P_NUM 10000
#define TPB 256
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


// particles' position & velocity
__device__ void generate_uniform_random_number(unsigned seed, int i, double* rho1, double* rho2, double* rho3 ){
    curandState state;
    curand_init(seed, i, 0, &state);
    // 0 -1 range
    *rho1= curand_uniform(&state);
    *rho2= curand_uniform(&state);
    *rho3= curand_uniform(&state);
}


__global__ void initial_position_velocity(unsigned seed, struct Particle *particles)
{

//  extern __shared__ struct *particles;

    int i= blockIdx.x*blockDim.x + threadIdx.x;
//	printf("init thread idx is: blockIdx=%d, blockDim=%d, threadIdx=%d\n", blockIdx.x, blockDim.x, threadIdx.x); 

	if(i<P_NUM){
//	printf("thread id is: %d\n", i);  

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
        planet = true   ;  // silicate particle
    else
        planet = false ;  // iron particle
    
    
    if (planet == true) {   
        if (particles[i].p_type == true) {
            // position initialization for outer shell, silicate particles
            generate_uniform_random_number(seed, i, &rho1, &rho2, &rho3 );
            miu= 1- 2 * rho2;   
            particles[i].position.x = cbrt(pow(R1,3.0)+(pow(R2,3.0)-pow(R1,3.0))*rho1) * sqrt(1-pow(miu,2.0))*cos(2*PI*rho3)+ CENTER_MASS_X;
            particles[i].position.y = cbrt(pow(R1,3.0)+(pow(R2,3.0)-pow(R1,3.0))*rho1) * sqrt(1-pow(miu,2.0))*sin(2*PI*rho3)+ CENTER_MASS_Z;
            particles[i].position.z = cbrt(pow(R1,3.0)+(pow(R2,3.0)-pow(R1,3.0))*rho1) * miu;
        }
        else {
            // position initialization for inner core, iron particles   
            generate_uniform_random_number(seed, i, &rho1, &rho2, &rho3 );
            miu= 1- 2 * rho2;
            particles[i].position.x = R * cbrt(rho1) * sqrt(1-pow(miu,2.0)) * cos(2*PI*rho3) + CENTER_MASS_X;
            particles[i].position.y = R * cbrt(rho1) * sqrt(1-pow(miu,2.0)) * sin(2*PI*rho3) + CENTER_MASS_Z;
            particles[i].position.z = R * cbrt(rho1) * miu;
        }

        // velocity initialization
        particles[i].velocity.x = Init_V;
        particles[i].velocity.y = 0;
        particles[i].velocity.z = 0;
        // calculate the distance r_xz on the plane xz from the center of mass for INER
        float r_xz = sqrt(pow((particles[i].position.x + CENTER_MASS_X),2.0) + pow((particles[i].position.z + CENTER_MASS_Z),2.0));
        float theta = atan((particles[i].position.z + CENTER_MASS_Z) / (particles[i].position.x + CENTER_MASS_X));
        particles[i].velocity.x +=  OMEGA * r_xz * sin(theta);
        particles[i].velocity.z +=  -OMEGA * r_xz * cos(theta);
        particles[i].velocity.y +=  0;
    }
    else {
        if (particles[i].p_type == true) {
            // position initialization for outer shell, silicate particles
            generate_uniform_random_number(seed, i, &rho1, &rho2, &rho3 );
            miu= 1- 2 * rho2;   
            particles[i].position.x = cbrt(pow(R1,3.0)+(pow(R2,3.0)-pow(R1,3.0))*rho1) * sqrt(1-pow(miu,2.0))*cos(2*PI*rho3)- CENTER_MASS_X;
            particles[i].position.y = cbrt(pow(R1,3.0)+(pow(R2,3.0)-pow(R1,3.0))*rho1) * sqrt(1-pow(miu,2.0))*sin(2*PI*rho3)- CENTER_MASS_Z;
            particles[i].position.z = cbrt(pow(R1,3.0)+(pow(R2,3.0)-pow(R1,3.0))*rho1) * miu;
        }
        else {
            // position initialization for inner core, iron particles   
            generate_uniform_random_number(seed, i, &rho1, &rho2, &rho3 );
            miu= 1- 2 * rho2;
            particles[i].position.x = R * cbrt(rho1) * sqrt(1-pow(miu,2.0)) * cos(2*PI*rho3) - CENTER_MASS_X;
            particles[i].position.y = R * cbrt(rho1) * sqrt(1-pow(miu,2.0)) * sin(2*PI*rho3) - CENTER_MASS_Z;
            particles[i].position.z = R * cbrt(rho1) * miu;
        }

        // velocity initialization
        particles[i].velocity.x = -1*Init_V;
        particles[i].velocity.y = 0;
        particles[i].velocity.z = 0;
        // calculate the distance r_xz on the plane xz from the center of mass for INER
        float r_xz = sqrt(pow((particles[i].position.x - CENTER_MASS_X),2.0) + pow((particles[i].position.z - CENTER_MASS_Z),2.0));
        float theta = atan((particles[i].position.z - CENTER_MASS_Z) / (particles[i].position.x - CENTER_MASS_X));
        particles[i].velocity.x +=  OMEGA * r_xz * sin(theta);
        particles[i].velocity.z +=  -OMEGA * r_xz * cos(theta);
        particles[i].velocity.y +=  0;
    }
}
    
}

// interaction forces

__device__ double interactionForce(double radius, double sqr_r, int force_type, bool if_sep_dec){
// radius is the parapmeter 'r' in Table 2
    double M, K1, K2, KRP1, KRP2;
    double accel;
    if (force_type == 0) {  // accel for silicate particle, silicate-iron 
        M = Mfe;
        K1 = Ksi;
        K2 = Kfe;
        KRP1 = KRPsi;
        KRP2 = KRPfe;
    }
    else if (force_type == 1) { // accel for iron particle, iron-silicate
        M = Msi;
        K1 = Ksi;
        K2 = Kfe;
        KRP1 = KRPsi;
        KRP2 = KRPfe;
    }
    else if (force_type == 2) { // accel for silicate particle, silicate-silicate
        M = Msi;
        K1 = Ksi;
        K2 = Ksi;
        KRP1 = KRPsi;
        KRP2 = KRPsi;
    }
    else if (force_type == 3) { // accel for iron particle, iron-iron
        M = Mfe;
        K1 = Kfe;
        K2 = Kfe;
        KRP1 = KRPfe;
        KRP2 = KRPfe;
    }
    else {
      //  fprintf(stderr,"invalid force type");
    }


    if (radius >= D){
        accel = G * M / sqr_r;
    }  
    else if (radius >= D-D*SDPsi){
        accel = G * M / sqr_r - 0.5*(K1+K2)*(D*D-sqr_r) ;
    }

    else if (radius >= D-D*SDPfe && if_sep_dec==true){
        accel = G * M / sqr_r - 0.5*(K1+K2)*(D*D-sqr_r) ;
    }
    else if (radius >= D-D*SDPfe && if_sep_dec==false){
        accel = G * M / sqr_r - 0.5*(K1*KRP1+K2)*(D*D-sqr_r) ;
    }
    else if (radius >= Epsilon  && if_sep_dec==true){
        accel = G * M / sqr_r - 0.5*(K1+K2)*(D*D-sqr_r) ;
    }
    else if (radius >= Epsilon  && if_sep_dec==false){
        accel = G * M / sqr_r - 0.5*(K1*KRP1+K2*KRP2)*(D*D-sqr_r) ;
    }
    else if(radius < Epsilon){
        radius=Epsilon;
        accel = G * M / sqr_r - 0.5*(K1+K2)*(D*D-sqr_r) ;
    }
    else{
        // no case fit
        // fprintf(stderr,"interaction force can't be calculated");
    }
//	if(accel > 100.0){
//		printf("acc calc is: %f\n", accel);
//		printf("radius is: %f\n", radius);
//	}
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
//  printf("BBI radius is %f\n", radius);

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
  
  if(accel > 100.0){
	printf("acc calc is: %f, radius is: %f, pj is: %f, %f, %f, pi is: %f, %f, %f\n", accel, radius, pj.position.x, pj.position.x, pj.position.x, pi.position.x, pi.position.x, pi.position.x);
//	printf("%f, ", radius, pj.position.x, pj.position.x, pj.position.x);
//	printf("%f, ", radius, pi.position.x, pi.position.x, pi.position.x);
  }

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
//	printf("in tile: %d\n", i);
//	printf("before BBI acc is: %f, %f, %f \n", acc.x, acc.y, acc.z);
//	printf("my particle %d is : %f, %f, %f \n", i, myParticle.position.x, myParticle.position.y, myParticle.position.z);
//	printf("shared particle %d is : %f, %f, %f \n", i, shParticles[i].position.x, shParticles[i].position.y, shParticles[i].position.z);
	if (isnan(myParticle.position.x) || isnan(shParticles[i].position.x)) 
		printf("i is %d\n", i);
	else
    	acc = bodyBodyInteraction(myParticle, shParticles[i], acc);
//	printf("acc is: %f\n", acc);
  }
  return acc;
}

//__global__ void calculate_forces(void *devP, void *devA)
//{
//  extern __shared__ struct Particle shParticles[];
//  struct Particle *global_P = (struct Particle *)devP;
//  struct vecfloat3 *global_A = (struct vecfloat3 *)devA;
//  struct Particle myParticle;
//  int i, tile;
//  struct vecfloat3 acc = {0.0f, 0.0f, 0.0f};
//  int gtid = blockIdx.x * blockDim.x + threadIdx.x;
//  myParticle = global_P[gtid];
//  for (i = 0, tile = 0; i < P_NUM; i += TPB, tile++) {
//    int idx = tile * blockDim.x + threadIdx.x;
//    shParticles[threadIdx.x] = global_P[idx];
//    __syncthreads();
//    acc = tile_calculation(myParticle, acc);
//    __syncthreads();
//  }
//  // Save the result in global memory for the integration step.
//  //float4 acc4 = {acc.x, acc.y, acc.z, 0.0f};
//  if(gtid < P_NUM){
//  global_A[gtid] = acc;
//  printf("thread gid is: %d\n", gtid);
//  } 
//
//}

__global__ void calculate_forces(struct Particle *devP, struct vecfloat3 *devA)
{
  extern __shared__ struct Particle shParticles[]; 
//  printf("blockDim.x is %d \n", blockDim.x);
  
//  struct Particle *global_P = (struct Particle *)devP;
//  struct vecfloat3 *global_A = (struct vecfloat3 *)devA;
  struct Particle myParticle;
  int i, tile;
  struct vecfloat3 acc; //= {0.0f, 0.0f, 0.0f}
  acc.x = 0.0f;
  acc.y = 0.0f;
  acc.z = 0.0f;
  int gtid = blockIdx.x * blockDim.x + threadIdx.x;
//  printf("thread gtid is: %d\n", gtid);
  if(gtid < P_NUM){
//	printf("testing!\n");
//	printf("GPU particle position is: %f, %f, %f\n", devP[gtid].position.x, devP[gtid].position.y, devP[gtid].position.z);
    myParticle = devP[gtid];
	if (isnan(myParticle.position.x)) 
		printf("gtid is %d\n", gtid);
  	for (i = 0, tile = 0; i < P_NUM; i += TPB, tile++) {
    	int idx = tile * blockDim.x + threadIdx.x;
//		printf("tile is : %d\n", tile);
		if(idx < P_NUM){
		// printf("thread idx is: %d\n", idx);
//		printf("shParticles thread idx is: %d, test point 1\n", threadIdx.x);
    	shParticles[threadIdx.x] = devP[idx];
//		printf("after shParticles thread idx is: %d, test point 2\n", threadIdx.x);
    	__syncthreads();
    	acc = tile_calculation(myParticle, acc);
//		printf("CF acc is: %f, %f, %f\n", acc.x, acc.y, acc.z);
    	__syncthreads();
		}
  	}
  // Save the result in global memory for the integration step.
  //float4 acc4 = {acc.x, acc.y, acc.z, 0.0f};

  	devA[gtid] = acc;
 // 	printf("thread gid is: %d\n", gtid);
  } 

}

__global__ void update_pos(struct Particle *devP, struct vecfloat3 *devA)
{
//  struct Particle *global_P = (struct Particle *)devP;
//  struct vecfloat3 *global_A = (struct vecfloat3 *)devA;
 
  int i= blockIdx.x*blockDim.x + threadIdx.x;
//  printf("update pos thread id is: %d\n", i); 

//  if(i<P_NUM){
//  global_P[i].position.x = global_P[i].position.x + global_P[i].velocity.x * Time_step + global_A[i].x / 2.0f * pow(Time_step, 2.0); 
//  global_P[i].position.y = global_P[i].position.y + global_P[i].velocity.y * Time_step + global_A[i].y / 2.0f * pow(Time_step, 2.0); 
//  global_P[i].position.z = global_P[i].position.z + global_P[i].velocity.z * Time_step + global_A[i].z / 2.0f * pow(Time_step, 2.0); 
//  }
  if(i<P_NUM){
  	devP[i].position.x = devP[i].position.x + devP[i].velocity.x * Time_step + devA[i].x / 2.0f * pow(Time_step, 2.0); 
  	devP[i].position.y = devP[i].position.y + devP[i].velocity.y * Time_step + devA[i].y / 2.0f * pow(Time_step, 2.0); 
  	devP[i].position.z = devP[i].position.z + devP[i].velocity.z * Time_step + devA[i].z / 2.0f * pow(Time_step, 2.0); 
  }
	
}

__global__ void update_vel(struct Particle *devP, struct vecfloat3 *devA, struct vecfloat3 *devA_next)
{
//  struct Particle *global_P = (struct Particle *)devP;
//  struct vecfloat3 *global_A = (struct vecfloat3 *)devA;
//  struct vecfloat3 *global_A_next = (struct vecfloat3 *)devA_next;

  int i= blockIdx.x*blockDim.x + threadIdx.x;

//  if(i<P_NUM){
//  global_P[i].velocity.x = global_P[i].velocity.x + (global_A_next[i].x + global_A[i].x) / 2.0f * Time_step; 
//  global_P[i].velocity.y = global_P[i].velocity.y + (global_A_next[i].y + global_A[i].y) / 2.0f * Time_step; 
//  global_P[i].velocity.z = global_P[i].velocity.z + (global_A_next[i].z + global_A[i].z) / 2.0f * Time_step; 
//  }

  if(i<P_NUM){
  	devP[i].velocity.x = devP[i].velocity.x + (devA_next[i].x + devA[i].x) / 2.0f * Time_step; 
  	devP[i].velocity.y = devP[i].velocity.y + (devA_next[i].y + devA[i].y) / 2.0f * Time_step; 
  	devP[i].velocity.z = devP[i].velocity.z + (devA_next[i].z + devA[i].z) / 2.0f * Time_step; 
  }


}

__global__ void print_devP(struct Particle *devP, struct vecfloat3 *devA)
{
//  struct Particle *global_P = (struct Particle *)devP;
//  struct vecfloat3 *global_A = (struct vecfloat3 *)devA;
//  struct vecfloat3 *global_A_next = (struct vecfloat3 *)devA_next;

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
//	printf("particle initialized \n");
//    cudaMemcpy(cpuA, devA, P_NUM*sizeof(struct vecfloat3),cudaMemcpyDeviceToHost);
	
//	return devA;
}


void particle_update(struct Particle *cpuP)
{
	struct Particle *devP=0;
	struct vecfloat3 *devA=0;
	struct vecfloat3 *devA_next=0;
	checkCudaErrors(cudaMalloc(&devP, P_NUM*sizeof(struct Particle)));
	checkCudaErrors(cudaMalloc(&devA, P_NUM*sizeof(struct vecfloat3)));
	checkCudaErrors(cudaMalloc(&devA_next, P_NUM*sizeof(struct vecfloat3)));

//	printf("particle updating, test point 1\n");
	
//	printf("before particle position is: %f, %f, %f\n", cpuP[10].position.x, cpuP[10].position.y, cpuP[10].position.z);
 
	checkCudaErrors(cudaMemcpy(devP, cpuP, P_NUM*sizeof(struct Particle),cudaMemcpyHostToDevice));
	cudaDeviceSynchronize();

//	print_devP<<<(P_NUM+TPB-1)/TPB, TPB>>>(devP, devA);
 

//	printf("particle updating, test point 2\n");

    // update position
    calculate_forces<<<(P_NUM+TPB-1)/TPB, TPB, TPB*sizeof(struct Particle)>>>(devP, devA);
    getLastCudaError("Kernel execution failed");  	// check if kernel execution generated and error
    cudaDeviceSynchronize();

//	printf("particle updating, test point 3\n");

    update_pos<<<(P_NUM+TPB-1)/TPB, TPB>>>(devP, devA);
	getLastCudaError("Kernel execution failed");  	// check if kernel execution generated and error
    cudaDeviceSynchronize();
	
//	printf("particle updating, test point 4\n");

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
//	printf("particle updated \n");
//    cudaMemcpy(cpuA, devA, P_NUM*sizeof(vecfloat3),cudaMemcpyDeviceToHost);

//    struct vecfloat3 *devA_temp = (struct vecfloat3 *)devA;
//    devA = devA_next;
//    devA_next = devA_temp;
}

void cuda_release()
{
//    cudaFree(devP);
//    cudaFree(devA);
//    cudaFree(devA_next);
}


// void launch_kernel_initial_position_velocity(unsigned seed, void *particles)
// {
// 	initial_position_velocity<<<(P_NUM+TPB-1)/TPB, TPB>>>(seed, particles);
// }

// void launch_kernel_calculate_forces(void *devP, void *devA)
// {
// 	calculate_forces<<<(P_NUM+TPB-1)/TPB, TPB>>>(devP, devA);
// }

// void launch_kernel_update_pos(void *devP, void *devA)
// {
// 	update_pos<<<(P_NUM+TPB-1)/TPB, TPB>>>(devP, devA);
// }

// void launch_kernel_update_vel(void *devP, void *devA, void *devA_next)
// {
// 	update_vel<<<(P_NUM+TPB-1)/TPB, TPB>>>(devP, devA, devA_next);
// }

