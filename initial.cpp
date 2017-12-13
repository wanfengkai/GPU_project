# include <math.h>
# include <stdio.h>
# include <string.h>
# include <cuda_runtime.h>
# include <cstdlib>

// universal gravitational constant in km
# define G  6.67408E-20
# define Epsilon 47.0975
# define D 376.78
# define Msi 7.4161E19   
# define Mfe 1.9549E20
# define Ksi 2.9114E14
# define Kfe 5.8228E14
# define KRPsi 0.01
# define KRPfe 0.02
# define SDPsi 0.001
# define SDPfe 0.002
# define Time_step 5.8117
# define R 100 // teacher will post it later
# define R1 200 // teacher will post it later
# define R2 300 // teacher will post it later
# define PI 3.14
# define Init_V 3.2416
# define CENTER_MASS_X 3185.5
# define CENTER_MASS_Z 9042.7
# define OMEGA 10 // wait for teacher
struct Particle
{
	float3 position_iner;
	float3 position_outer;
	float3 velocity_iner;
	float3 velocity_outer;
	bool earthmoon;   // True :earth  False :moon

};


// interaction forces

void interactionForce(double* radius, double* Force,bool if_sep_dec){
// radius is the parapmeter 'r' in Table 2

	if (radius >= D){
	Force= G * Msi* Mfe* sqrt(radius);
	}  
	else if (radius >= D-D*SDPsi){
        Force= G * Msi * Mfe*sqrt(radius)- 0.5*(Ksi+Kfe)*(pow(D,2)-pow(radius,2)) ;
	}

        else if (radius >= D-D*SDPfe &&if_sep_dec==true){
        Force= G * Msi * Mfe*sqrt(radius)- 0.5*(Ksi+Kfe)*(pow(D,2)-pow(radius,2)) ;
	}
        else if (radius >= D-D*SDPfe &&if_sep_dec==false){
        Force= G * Msi * Mfe*sqrt(radius)- 0.5*(Ksi*KRPsi+Kfe)*(pow(D,2)-pow(radius,2)) ;
        }
        else if (radius>= Epsilon  &&if_sep_dec==true){
        Force= G * Msi * Mfe*sqrt(radius)- 0.5*(Ksi+Kfe)*(pow(D,2)-pow(radius,2)) ;
        }
        else if (radius>= Epsilon  &&if_sep_dec==false){
        Force= G * Msi * Mfe*sqrt(radius)- 0.5*(Ksi*KRPsi+Kfe*KRPfe)*(pow(D,2)-pow(radius,2)) ;
        }
        else if(radius< Epsilon){
        radius=Epsilon;
        Force= G * Msi * Mfe*sqrt(radius)- 0.5*(Ksi+Kfe)*(pow(D,2)-pow(radius,2)) ;
        }
        else{
            // no case fit
            fprintf(stderr,"interaction force can't be calculated");
        }
}

// particles' position & velocity
void generate_uniform_random_number(double* rho1, double* rho2, double* rho3 ){
	// 0 -1 range
	rho1= rand()%1;
	rho2= rand()%1;
	rho3= rand()%1;
}


__global__ void initial_position_velocity(struct Particle *array)
{

	int i= blockIdx.x*blockDim.x + threadIdx.x;


	generate_uniform_random_number(double* rho1, double* rho2, double* rho3 );
        // position update for inert shell
	miu= 1- 2 * rho2;
        array[i].position_iner.x = R * cbrt(rho1)*sqrt(1-pow(miu,2))*cos(2*PI*rho3);
        array[i].position_iner.y = R * cbrt(rho1)*sqrt(1-pow(miu,2))*sin(2*PI*rho3);
        array[i].position_iner.z = R * cbrt(rho1)*miu;
	// position update for outer shell
	// generate again
	generate_uniform_random_number(double* rho1, double* rho2, double* rho3 );
	miu= 1- 2 * rho2;	
        array[i].position_outer.x = cbrt(pow(R1,3)+(pow(R2,3)-pow(R1,3))*rho1) * sqrt(1-pow(miu,2))*cos(2*PI*rho3);
        array[i].position_outer.y = cbrt(pow(R1,3)+(pow(R2,3)-pow(R1,3))*rho1) * sqrt(1-pow(miu,2))*sin(2*PI*rho3);
        array[i].position_outer.z = cbrt(pow(R1,3)+(pow(R2,3)-pow(R1,3))*rho1) *miu;

	// velocity initialize

	// sign for + or - 1
	if (rand()%1>0.5){
		sign=1	;
	}
	else{
		sign=-1 ;
	}
        array[i].velocity_iner.x =Init_V*sign;
	array[i].velocity_iner.y =0;
	array[i].velocity_iner.z =0;

        array[i].velocity_outer.x =-1*Init_V*sign;
	array[i].velocity_outer.y =0;
	array[i].velocity_outer.z =0;	

	
	// calculate the distance r_xz on the plane xz from the center of mass for INER
	float r_xz = sqrt(pow((array[i].position_iner.x-CENTER_MASS_X),2)+pow((array[i].position_iner.z-CENTER_MASS_Z),2));
	float theta = atan((array[i].position_iner.z-CENTER_MASS_Z)/(array[i].position_iner.x-CENTER_MASS_X));
	array[i].velocity_iner.x +=  OMEGA * r_xa* sin(theta);
	array[i].velocity_iner.z +=  -OMEGA * r_xa* cos(theta);
	array[i].velocity_iner.y +=  0;

	// calculate the distance r_xz on the plane xz from the center of mass for OUTER
	float r_xz = sqrt(pow((array[i].position_outer.x+CENTER_MASS_X),2)+pow((array[i].position_outer.z+CENTER_MASS_Z),2));
	float theta = atan((array[i].position_outer.z+CENTER_MASS_Z)/(array[i].position_outer.x+ CENTER_MASS_X));
	array[i].velocity_outer.x +=  OMEGA * r_xa* sin(theta);
	array[i].velocity_outer.z +=  -OMEGA * r_xa* cos(theta);
	array[i].velocity_outer.y +=  0;
	
}










