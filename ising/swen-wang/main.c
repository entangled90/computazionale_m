#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "dynamicList.h"

void spin_init (short int * configuration){
	int i,j;
	double tmp;
	for ( i = 0; i< WIDTH; i++){
		for ( j = 0; j<HEIGHT ; j++){
			tmp = rand()/(double)RAND_MAX;
			if ( tmp > 0.5){
				configuration[i*WIDTH+j] = +1;
			}
			else {
				configuration[i*WIDTH+j] = -1;
			}
		}
	}
}

inline void savePPM(short int * x)
{
    unsigned char white[3] = {255,255,255};
    unsigned char black[3] = {0,0,0};
   FILE *f = fopen("image.ppm", "wb");
    fprintf(f, "P6\n%d %d\n255\n", WIDTH, HEIGHT);
    int i,j;
    for(i = WIDTH-1; i >= 0; i--)
        for(j = 0; j < WIDTH; j++){
            if(x[i*WIDTH+j] == 1)
                fwrite(white, sizeof(unsigned char), 3, f);
            if(x[i*WIDTH+j] == -1)
                fwrite(black, sizeof(unsigned char), 3, f);
        }
    fclose(f);
}


/*Calcola la magnetizzazione */
inline double magnetization(short int *x){
	int i,j;
	int mag=0;
	for (i=0;i<WIDTH;i++){
		for(j=0;j<WIDTH;j++){
			mag += x[i*WIDTH+j];
		}
	}
	return mag/((double) (WIDTH*WIDTH));
}


int main (int argc, char *argv[]){	
	int i;
	int iteration = 0;
	double mag;
	double sum = 0;
	double mag_prevista = ( 1 - pow(sinh(2*BETA*J),-4));
	short int  * configuration;
	srand(time(NULL));
//	rlxs_init(1,time(NULL));
/*
	for ( i=0 ; i<5000;i++){
		mersenne_generate();
	}
*/
	FILE *f_mag = fopen("data/magnetization.dat","w+");
	configuration = malloc(WIDTH*WIDTH*sizeof(short int));
	spin_init(configuration);
	while (iteration < ITERATION_MAX){
		// do stuff
	}


	fclose(f_mag);
	printf("Magnetizzazione media: %lf\n",sum/((double) ITERATION_MAX));
	exit(EXIT_SUCCESS);
}