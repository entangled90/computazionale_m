#ifndef SW_H

#define SW_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cluster.h"

#define N 32 
#define BETA 0.4
#define J 1 

void spin_init (spin * s){
	int i,j;
	double tmp;
	for ( i = 0; i< N; i++){
		for ( j = 0; j<N ; j++){
			tmp = rand()/(double)RAND_MAX;
			if ( tmp > 0.5){
				s[i*N+j].s = +1;
			}
			else {
				s[i*N+j].s = -1;
			}
			s[i*N+j].cluster = 0;
			s[i*N+j].i = i;
			s[i*N+j].j = j;
		}
	}
}



/*Calcola la magnetizzazione */
inline double magnetization(spin *x){
	int i,j;
	int mag=0;
	for (i=0;i<N;i++){
		for(j=0;j<N;j++){
			mag += x[i*N+j].s;
		}
	}
	return mag/((double) (N*N));
}


/*Prima volta che si clusterizza*/
void  clusterize (spin * matrix, cluster * head){
	int i,j;
	int dx,dy;
	double prob;
	for(i = 0; i<N;i++){
		for ( j=0;j<N;j++){
			for ( dx = -1; dx<2;dx+=2){
				for ( dy = -1 ; dy<2 ;dy+=2){
					if ( matrix[((i+dx+N)%N)*N + (j+dy+N)%N].s == matrix[ i*N+j].s){
						prob = 1.0 - exp(-BETA*J);
						if( rand()/((double)RAND_MAX) < prob ){
								// Li metto nello stesso cluster
						}
					}
				}
			}
		}
	}
	if(matrix[i*N+j].cluster == 0){
		new_cluster( head , matrix +(i*N+j));
	}
}

#endif