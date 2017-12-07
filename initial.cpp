# include <math.h>
# include <stdio.h>
# include <string.h>
# include <cuda_runtime.h>

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
 


// interaction forces

void interactionForce(double radius, double Force,bool if_sep_dec){
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









// particles' velocities








// unit alloc 
