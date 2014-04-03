#include <stdlib.h>
#include <stdio.h>
#include "constants-metro.h"
#include <math.h>
#include "mtwist.h"

#define METROPOLIS_C

// Inizializza gli spin
void spin_init (short int * configuration, int N){
	int i,j;
	double tmp;
	for ( i = 0; i< N; i++){
		for ( j = 0; j<N ; j++){
			tmp = mt_drand()/(double)RAND_MAX;
			if ( tmp > 0.5){
				configuration[i*N+j] = +1;
			}
			else {
				configuration[i*N+j] = -1;
			}
		}
	}
}


//Salva lo stato del sistema come immagine .ppm in bianco/nero
inline void savePPM(short int * x, int N)
{
    unsigned char white[3] = {255,255,255};
    unsigned char black[3] = {0,0,0};
   FILE *f = fopen("image.ppm", "wb");
    fprintf(f, "P6\n%d %d\n255\n", N, N);
    int i,j;
    for(i = N-1; i >= 0; i--)
        for(j = 0; j < N; j++){
            if(x[i*N+j] == 1)
                fwrite(white, sizeof(unsigned char), 3, f);
            if(x[i*N+j] == -1)
                fwrite(black, sizeof(unsigned char), 3, f);
        }
    fclose(f);
}

//calcola la magnetizzazione moltiplicata per il numero di siti
inline double magnetization(short int *x, int N){
	int i,j;
	int mag=0;
	for (i=0;i<N;i++){
		for(j=0;j<N;j++){
			mag += x[i*N+j];
		}
	}
	return mag;
}

// Calcola la variazione dell'hamiltoniana in seguito all'inversione dello spin nel sito (a,b)
inline double delta_Ham_ising ( short int * configuration,int a, int b, int N){
	double ham;
	ham = -J * ( configuration[ ((a+1)%N)*N + b]+ configuration[((a-1+N)%N)*N+b]
		 + configuration[a*N+(b+1)%N]+configuration[a*N + (b-1+N)%N])*(-2*configuration[a*N +b]); 
	return ham;
}


// Esegue uno step dell'algoritmo metropolis per ogni sito
inline void metropolis_ising( short int *x, int N , double BETA ){
	double dH;
	int i,j;
	for (i = 0; i<N;i++){
		for(j= 0; j<N;j++){
			dH = delta_Ham_ising(x,i,j,N);
			if( mt_drand() <  exp(-BETA*dH)){
				x[i*N+j] = - x[i*N+j];
			}
		}
	}
}

//Ritorna l'hamiltoniana. Viene calcolata partendo da in alto a sinistra e vengono usati i link destra e in basso.
inline double hamiltoniana ( short int * configuration, int N){
	double ham=0;
	int a,b;
	for (a = 0; a<N ; a++){
		for ( b= 0; b<N; b++){
			ham += -( configuration[ ((a+1+N)%N)*N + b]	+ configuration[a*N+(b+1+N)%N])*(configuration[a*N +b]); 
		}
	}
	return (ham);
}

/*Somma su righe per la correlazione*/
double sum_row(short int * configuration, int row, int N){
	double sum = 0;
	int j = 0;
	for (j = 0; j<N ;j++){
		sum += configuration[row*N+j];
	}
	return sum /= (double)N ;
}

/*Somma su colonne per la correlazione*/

double sum_col(short int * configuration, int col, int N){
	double sum = 0;
	int j = 0;
	for (j = 0; j<N ;j++){
		sum += configuration[j*N+col];
	}
	return sum /= (double)N ;
}

