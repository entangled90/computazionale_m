#include "constants.h"
#include "random.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define METROPOLIS_C


inline double delta_Ham_ising ( short int * configuration,int a, int b){
	double ham;
	ham = -J * ( configuration[ ((a+1)%WIDTH)*WIDTH + b]+ configuration[((a-1+WIDTH)%WIDTH)*WIDTH+b]
		 + configuration[a*WIDTH+(b+1)%HEIGHT]+configuration[a*WIDTH + (b-1+HEIGHT)%HEIGHT])*(-2*configuration[a*WIDTH +b]); 
	return ham;
}

inline void metropolis_ising( short int *x ){
	double dH;
	float *tmp = malloc(sizeof(float)*WIDTH*WIDTH) ;
	ranlxs(tmp,WIDTH*WIDTH);
	int i,j;
	for (i = 0; i<WIDTH;i++){
		for(j= 0; j<WIDTH;j++){
			dH = delta_Ham_ising(x,i,j);
			if( tmp[i*WIDTH+j] <  exp(-BETA*dH)){
				x[i*WIDTH+j] = - x[i*WIDTH+j];
			}
		}
	}
}



/*
inline void metropolis( short int *x ){
	double dH;
	float *tmp = malloc(sizeof(float)*N) ;
	unsigned int * new = malloc(sizeof(float)*N);
	ranlxs(tmp,N);
	int i,j;
	for (i = 0; i<WIDTH;i++){
		for(j= 0; j<WIDTH;j++){
			dH = delta_Ham(x,i,j);
			if( tmp[i*WIDTH+j] <  exp(-BETA*dH)){
				new[i*WIDTH+j] = - x[i*WIDTH+j];
			}
			else{
				new[i*WIDTH+j] = x[i*WIDTH+j];
			}
		}
	}
	for (i = 0; i<N;i++){
		x[i] = new[i];
	}
}
*/

inline double hamiltonian_ising ( short int * configuration){
	double ham=0;
	int a,b;
	for (a = 0; a<WIDTH ; a++){
		for ( b= 0; b<HEIGHT; b++){
			ham += J * ( configuration[ ((a+1+WIDTH)%WIDTH)*WIDTH + b]+ configuration[((a-1+WIDTH)%WIDTH)*WIDTH+b]
				+ configuration[a*WIDTH+(b+1+HEIGHT)%HEIGHT]+configuration[a*WIDTH + (b-1+HEIGHT)%HEIGHT])*(configuration[a*WIDTH +b]); 
		}
	}
	return (ham/2.0);
}

