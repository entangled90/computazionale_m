
/****************
 *  COMANDO
 gcc -lm -Wall -O3 -funroll-loops  -o main main.c
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>

/*Numero di dimensioni */
#define N 3
#define TERM_TIME 20000
#define MAX_COLLISION 50000

/*Numero particelle */
#define  NUMBER_OF_PARTICLES  250
/* Diametro sfere */
double SIGMA =  0;
/*Tavola delle collisioni */

/*  (i,j,tempo di collisione) */
unsigned int index_collision[2];
double time_collision = 0;
unsigned int numOfCollisions = 0;
double total_time = 0;
double DIST_RET=0;
double K_BOLTZ=1;
double pression = 0;
double D_speed_norm = 0;
double time_prec=0;

typedef struct particle_s {
	double position[N];
	double speed[N];
	} particle_s ;

void print_coordinate ( particle_s *particleList){
	FILE *f = fopen ( "data/pack.dat","w");
	int i ;
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		fprintf(f,"%e\t%e\t%e\n", particleList[i].position[0], particleList[i].position[1],particleList[i].position[2]);
	}
	fclose(f);
	}

void print_speed ( particle_s *particleList ){
	FILE *f = fopen ( "data/speed.dat","w");
	int i = 0;
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		fprintf(f,"%e\t%e\t%e\n", particleList[i].speed[0], particleList[i].speed[1], particleList[i].speed[2]);
	}
	fclose(f);
	}

void particle_init ( particle_s *particleList ){
	int i_part= 0;
	int i,j;
	double r_0[N];
	double speed_cm[3];
	int z_row = 0;
	double offsetX=0;
	double offsetY=0;
	
	for ( i = 0; i<N;i++) {
		r_0[i]=0;
		speed_cm[i]=0.0;
	}
	
	while (i_part <NUMBER_OF_PARTICLES){
		for ( i=0; i<N;i++){
		particleList[i_part].position[i] = r_0[i];
		}
		//particleList[i_part].dist = 0;
		//particleList[i_part].num_collision=0;
	//printf("%d Printed: (%e,%e,%e)\n", i_part,r_0[0],r_0[1],r_0[2]);
	r_0[0] += DIST_RET;
	if ( r_0[0] > 1.0 - SIGMA + offsetX){
		r_0[0] = offsetX;
		r_0[1] += DIST_RET;
		if ( r_0[1] > 1.0 - SIGMA + offsetY){
			z_row++;
			offsetX= (z_row%2)*DIST_RET/2.0 ;
			r_0[0] = offsetX;
			offsetY= (z_row%2)*DIST_RET/2.0;
			r_0[1] = offsetY;
			r_0[2] += DIST_RET/sqrt(2.0);
		}
	}
	if (r_0[2]> 1 - SIGMA){
		printf("%e\n",r_0[2]);
		print_coordinate(particleList);
		printf("Impacchettamento non completato\n");
		exit(1);
	}
	i_part++;
	}
	for ( i = 0; i< NUMBER_OF_PARTICLES; i++){
		for ( j = 0; j<N;j++){
			particleList[i].speed[j] = 2*(rand()/(RAND_MAX*1.0)) - 1.0 ;
			speed_cm[j] += particleList[i].speed[j];
			}
	}
	for (i =0 ; i< NUMBER_OF_PARTICLES; i++){
		for ( j = 0; j<N;j++){
			particleList[i].speed[j] -= (speed_cm[j]/((double) NUMBER_OF_PARTICLES));
		}
	}
}

inline void  sum ( double * v1 , double * v2, double * v_sum){
	int i = 0;
	for (i = 0; i < N ; i++){
		v_sum[i] = v1[i] + v2[i];
	}
	}

inline void diff (double * v1 , double *v2, double * v_diff){
	int i = 0;
	for (i = 0; i < N ; i++){
		v_diff[i] = v1[i] - v2[i];
	}
}

inline double scalar_prod ( double * v1 , double * v2){
	int i = 0;
	double sum  = 0;
	for ( i = 0; i< N ; i++){
		sum += v1[i]* v2[i];
	}
	return (sum);
	}
inline void scalar_mult  ( double scalar , double* vec){
	int i = 0;
	for( i = 0 ; i< N ; i++){
		vec[i] *= scalar;
	}
	}

void check_distance (particle_s *particleList){
	int i,j;
	double distance = 0;
	double diff_v[N];
	int x,y,z;
	particle_s temp_part;
	//particle_s temp_part;
	for (i = 0 ; i< NUMBER_OF_PARTICLES ; i++){
		for(j = i+1;j <NUMBER_OF_PARTICLES ; j++){
			for ( x= -1; x < 2 ; x++){
				for ( y = -1; y<2 ; y++){
					for(z=-1;z<2;z++){
						temp_part = particleList[j];
						temp_part.position[0] += x;
						temp_part.position[1] += y;
						temp_part.position[2] += z;
						diff(particleList[i].position,temp_part.position,diff_v);
						distance = sqrt(scalar_prod(diff_v,diff_v));
						if( distance <SIGMA){
							printf("Sfere (%d,%d) troppo vicine!\n",i,j);
						}
					}
				}
			}			
		}
	}
}


/* Calcola la tabella delle collisioni (solo la diagonale superiore)*/
 void collision_table (particle_s *particleList, double * collTable){
	int i,j;
	double x,y,z;
	double min= DBL_MAX;
	double r_diff[N];
	double v_diff[N];
	double det;
	double temp;
	double rv,rr,vv;
	particle_s temp_part;
	for (i = 0; i < NUMBER_OF_PARTICLES ; i++){
		for ( j= i+1; j < NUMBER_OF_PARTICLES; j++){
			min = DBL_MAX;
			for ( x= -1; x < 2 ; x++){
				for ( y = -1; y<2 ; y++){
					for ( z=-1; z<2; z++){
						temp_part = particleList[j];
						temp_part.position[0] += x;
						temp_part.position[1] += y;
						temp_part.position[2] += z;
						diff(particleList[i].position,temp_part.position, r_diff);
						diff( particleList[i].speed,temp_part.speed, v_diff);
						rv = scalar_prod( r_diff, v_diff);
						if( rv < 0){
							rr = scalar_prod(r_diff,r_diff);
							vv = scalar_prod(v_diff,v_diff);
							det = rv*rv - vv*( rr -SIGMA*SIGMA);
							if (det > 0){
								temp =  ( - rv - sqrt( det ))/ vv;
								if ( temp < min ){
									min = temp;
								}
							}
						}
					}
				}
			}
			collTable[i*NUMBER_OF_PARTICLES + j] = min;
		}
	}
}
inline void search_min_coll ( particle_s * particleList, double *collTable){
 int i,j;
 double min = DBL_MAX;
 for ( i= 0; i<NUMBER_OF_PARTICLES;i++){
	for(j=i+1; j<NUMBER_OF_PARTICLES;j++){
		if(  collTable[i*N+j] < min){
			min = collTable[i*N+j];
			index_collision[0] = i;
			index_collision[1] = j;
		}
	}
 }
time_collision = min;
}
inline void step (particle_s *particleList){
	int i,j;
	for ( i = 0; i < NUMBER_OF_PARTICLES ; i++){
		for (j =0 ; j< N ;j++){
			particleList[i].position[j] += time_collision*particleList[i].speed[j];
		}
	}
}
void substract_t0 (particle_s * particleList, double * collTable){
	int i,j;
	for ( i = 0; i<NUMBER_OF_PARTICLES;i++){
		for ( j = i+1; j<NUMBER_OF_PARTICLES;j++){
			collTable[i*N+j] -= time_collision;
		}
	}
}

void switch_speeds(particle_s *particleList){
	int j;
	int  x,y,z;
	double r_diff[N];
	double v_diff[N];
	double temp_r_diff[N];
	/* r_diff = R0 - R1
	 * v_diff = V0 _ V1
	 */
	double min = DBL_MAX;
	double tmp_dbl;
	particle_s temp_part;
	double v_temp;
	for ( x= -1; x < 2 ; x++){
		for ( y = -1; y<2 ; y++){
			for( z=-1; z<2 ; z++){
				temp_part = particleList[index_collision[1]];
				temp_part.position[0] += x;
				temp_part.position[1] += y;
				temp_part.position[2] += z;
				diff(particleList[index_collision[0]].position,temp_part.position, r_diff);
				tmp_dbl = scalar_prod(r_diff,r_diff) ; 
				if ( tmp_dbl < min){
					min = tmp_dbl;
					for ( j= 0; j<N; j++){
					temp_r_diff[j ] = r_diff[j];
					}
				}
			}
		}
	}
	diff( particleList[index_collision[0]].speed, particleList[index_collision[1]].speed, v_diff);
	scalar_mult( 1/(sqrt(scalar_prod(temp_r_diff,temp_r_diff))), temp_r_diff);
	for ( j = 0 ; j < N ; j++){
		v_temp = scalar_prod(v_diff,temp_r_diff)*temp_r_diff[j];
		particleList[index_collision[0]].speed[j] -= v_temp*temp_r_diff[j];
		particleList[index_collision[1]].speed[j] += v_temp*temp_r_diff[j];
	}
}

/* Aggiorna i tempi di collisioni per le righe e le colonne della matrice della particelle che hanno colliso:
 *
 *
 * NOTA BENE :
 * Modifica le righe associate ad una particella -> LA MATRICE é simmetrica 
 *
 */


/*
 void update_coll_table(particle_s *particleList, double * collTable){
	int i,j;
	double r_diff[N];
	double v_diff[N];
	//double det;
	double x,y,z;
	double min= DBL_MAX;
	double det;
	double temp;
	double rv,vv;
	int a,b,c;
	particle_s temp_part;
	for (i = 0; i < 2 ; i++){
		//a,b indici di riga e colonna -> Matrice simmetrica: tengo solo b>a, ossia j> index_collision[i] 
		for ( j= 0 ; j < NUMBER_OF_PARTICLES; j++){
			min = DBL_MAX;
			a=index_collision[i];
			b=j;
			if( a != b){
				if( a > b){
					c=a;
					a=b;
					b=c;
				}
				for ( x= -1; x < 2 ; x++){
					for ( y = -1; y<2 ; y++){
						for(z=-1; z<2 ; z++){
							temp_part = particleList[b];
							temp_part.position[0] += x;
							temp_part.position[1] += y;
							temp_part.position[2] += z;
							diff(particleList[a].position,temp_part.position, r_diff);
							diff( particleList[a].speed,temp_part.speed, v_diff);
							rv = scalar_prod( r_diff, v_diff);
							if( rv < 0){
								vv = scalar_prod(v_diff,v_diff);
								det = rv*rv - vv*( scalar_prod(r_diff,r_diff) -SIGMA*SIGMA);
								if (det > 0){
									temp = time_prec + ( - rv - sqrt( det ))/ vv ;
									if ( temp < min ){
										min = temp;
									}
								}
							}
						}
					}
				}
			collTable[a*NUMBER_OF_PARTICLES+b]= min;
			}
		}
	}
}
*/

double calc_min ( particle_s * particleList,int i, int j){
	double x,y,z;
	double min= DBL_MAX;
	double r_diff[N];
	double v_diff[N];
	double det;
	double temp;
	particle_s temp_part;
	for ( x= -1; x < 2 ; x++){
		for ( y = -1; y<2 ; y++){
			for ( z = -1; z<2;z++){
				temp_part = particleList[j];
				temp_part.position[0] += x;
				temp_part.position[1] += y;
				temp_part.position[2] += z;
				diff(particleList[i].position,temp_part.position, r_diff);
				diff( particleList[i].speed,temp_part.speed, v_diff);
				if( scalar_prod( r_diff, v_diff) < 0){
					det = scalar_prod(r_diff,v_diff)*scalar_prod(r_diff,v_diff) - scalar_prod(v_diff,v_diff)*( scalar_prod(r_diff,r_diff) -SIGMA*SIGMA);
					if (det > 0){
						//printf("Scalar prod: %e \n",(scalar_prod(v_diff,v_diff)));
						temp = ( - scalar_prod( r_diff, v_diff) - sqrt( det ))/ (scalar_prod(v_diff,v_diff) );
						if ( temp < min ){
							min = temp;
						}	
					}
				}	
			}
		}
	}
	return min;
}

void update_coll_table(particle_s *particleList, double * collTable){
	int i,j;
	int a,b,c;
	for (i = 0; i < 2 ; i++){
		//a,b indici di riga e colonna -> Matrice simmetrica: tengo solo b>a, ossia j> index_collision[i] 
		for ( j= 0 ; j < NUMBER_OF_PARTICLES; j++){
			a=index_collision[i];
			b=j;
			if( a != b){
				if( a > b){
					c=a;
					a=b;
					b=c;
				}
				
			collTable[a*NUMBER_OF_PARTICLES+b]= calc_min(particleList,a,b);
			}
		}
	}
}


inline void fix_boundaries (particle_s *particleList){
	int i = 0;
	int j = 0;
	for (i = 0 ; i< NUMBER_OF_PARTICLES ; i++){
		for( j= 0; j< N ; j++){
			particleList[i].position[j] -= floor(particleList[i].position[j]);
			}
		}
}

double kin_en ( particle_s *particleList) {
	int i = 0;
	double  sum = 0;
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		sum += scalar_prod(particleList[i].speed, particleList[i].speed);
	}
	return sum;
	}
double total_momentum (particle_s *particleList){
	int i,j;
	double  sum = 0;
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		for ( j = 0; j< N ; j++){
		sum += particleList[i].speed[j];
		}
	}
	return sum;
	}




void histo( const char * input , int n){
    double width;
    int* freq;
    int i ;
 //   double mean;
  //  double sigma;
    float tmp;
    double max; //=DBL_MIN;
    double min;//=DBL_MAX;
    freq =  malloc(n*sizeof(int));
 	for(i=0;i<n;i++){
		freq[i]=0;
	}
	FILE * f =fopen(input,"r");;
	max=-1;
	min=1;
	width=(max-min)/(double)n;
	while( fscanf(f,"%e\n",&tmp) == 1){
		for(i= 0;i<n;i++){
		  if( (tmp>min+i*width) && (tmp<=min+(i+1)*width) ){
		 	 freq[i]++;
		  }
		}
	}
	fclose(f);
	/* Necessario per la creazione del file che gnuplot può fittare, diviso negli opportuni intervalli*/
	f=fopen("data/vx-histo.dat","w");
	for(i=0;i<n;i++){
        fprintf(f,"%lf\t%d\t\n",min+(i+0.5)*width,freq[i]);
    }
    fclose(f);
    free(freq);
}

void boltzmann_file_save ( particle_s *particleList){
	int i = 0;
	FILE *f = fopen("data/v2.dat","a");
	FILE *fx = fopen ("data/vx.dat","a");
	FILE *fy = fopen("data/vy.dat","a");
	FILE *fz= fopen ("data/vz.dat","a");
	for (i = 0; i< NUMBER_OF_PARTICLES ; i++){
		fprintf(fx,"%e\n",particleList[i].speed[0]);
		fprintf(fy,"%e\n",particleList[i].speed[1]);
		fprintf(fz,"%e\n",particleList[i].speed[2]);
		fprintf(f,"%e\n",scalar_prod(particleList[i].speed,particleList[i].speed));
	}
	fclose(fx);
	fclose(fy);
	fclose(fz);
	fclose(f);
}
void print_coll_table (particle_s *particleList, double * collTable){
	int i,j;
	for(i = 0; i<NUMBER_OF_PARTICLES; i++){
		for(j=i+1; j<NUMBER_OF_PARTICLES; j++){
		printf("Tempo collisione (%d,%d): %e\n", i,j,collTable[i*NUMBER_OF_PARTICLES+j]);
		}
	}
	}
	
	
inline void copylist ( particle_s * in , particle_s * out, int n){
	int i = 0;
	for ( i= 0; i< n;i++){
		out[i] = in[i];
	}
}
//(sqrt(N)-1)*DIST_RET = L-SIGMA
// DIST_RET = (L-SIGMA)/(Sqrt(N)-1)
int main (int argc, char *argv[]) {
//	int i,j;
particle_s * particleList;
particle_s * List_0;
double fraz_imp=0.1;
double tempoterm=0;
double deltaV_pre[N];
double deltaV_post[N];
double deltaV[N];
double * collTable;
if (argc == 2){
	fraz_imp = atof(argv[1]);
}
/*******************
** INIT VARI
*******************/
//Frazione di impacchettamento
srand(time(NULL));
SIGMA = cbrt(6*fraz_imp/ NUMBER_OF_PARTICLES / M_PI);
//DIST_RET = (1.0-SIGMA)/(sqrt(NUMBER_OF_PARTICLES)-1.0);
//DIST_RET = (1.0-SIGMA)*sqrt(1/sqrt(3)/(double) NUMBER_OF_PARTICLES);
DIST_RET = cbrt(6*0.61/ (NUMBER_OF_PARTICLES *M_PI));
printf("SIGMA = %e\n",SIGMA);
printf("Frazione di impacchettamento: %e\n", fraz_imp);
collTable = malloc (NUMBER_OF_PARTICLES*NUMBER_OF_PARTICLES*sizeof(double));
particleList = malloc ( NUMBER_OF_PARTICLES * sizeof(particle_s));
List_0 = malloc(NUMBER_OF_PARTICLES* sizeof(particle_s));
particle_init (particleList);
fix_boundaries(particleList);
check_distance(particleList);
print_coordinate(particleList);
printf("#Collisions: %d \n", numOfCollisions);
printf(" K = %e \t P= %e \n", kin_en(particleList), total_momentum(particleList));
/**** FILE************/
FILE *f = fopen("data/v2.dat","w+");
FILE *fx = fopen ("data/vx.dat","w+");
FILE *fy = fopen("data/vy.dat","w+");
FILE *fz= fopen ("data/vz.dat","w+");
fclose(fx);
fclose(fy);
fclose(fz);
fclose(f);
/*****************
TERMALIZZAZIONE
*****************/
collision_table(particleList, collTable);
while ( numOfCollisions < TERM_TIME){
	//time_prec = total_time;
	search_min_coll(particleList,collTable);
	step(particleList);
	total_time += time_collision;
	substract_t0(particleList,collTable);
	switch_speeds(particleList);
	fix_boundaries(particleList);
	update_coll_table(particleList,collTable);
	numOfCollisions++;
}
//copylist(particleList,List_0,NUMBER_OF_PARTICLES);
/***************
 * PRODUZIONE
 **************/
tempoterm=total_time;
printf("Termalizzato\n");
while (numOfCollisions < MAX_COLLISION){
	/* Calcola la matrice dei tempi e mette il minimo*/
	search_min_coll(particleList,collTable); 
	/* Muove le palline */
	step(particleList);
	diff(particleList[index_collision[0]].speed,particleList[index_collision[1]].speed,deltaV_pre);
	/*Modifica le velocità*/
	switch_speeds(particleList);
	fix_boundaries(particleList);
	total_time += time_collision;
	diff(particleList[index_collision[0]].speed,particleList[index_collision[1]].speed,deltaV_post);
	diff(deltaV_pre,deltaV_post,deltaV);
	substract_t0(particleList,collTable);
	update_coll_table(particleList,collTable);
	pression += sqrt(scalar_prod(deltaV,deltaV));
	if( numOfCollisions % 500 == 0 ){
			boltzmann_file_save(particleList); 
			/*
			for ( i = 0; i< NUMBER_OF_PARTICLES ;i++){
				for ( j = 0; j< NUMBER_OF_PARTICLES ;j++){
				sum +=  
				}
			}
		*/

	}
	if(numOfCollisions % 5000 == 0){
	printf("#Collisions: %d K=%e\t P=%e\n", numOfCollisions,kin_en(particleList),total_momentum(particleList));
	}
	numOfCollisions +=1;
}
// I termini della somma hanno un fattore comune, M=1
//NOTA CHE VALE IN DUE DIMENSIONI
pression*=SIGMA/((total_time-tempoterm)*2.0*kin_en(particleList));
pression+=1.0;
// DA MODIFICARE LA CORREZIONE!!!
pression*=fraz_imp/0.7405;
//NUMBER_OF_PARTICLES*3.0/2.0*kin_en()*fraz_imp/M_PI*2*sqrt(3.00);
// Stampa il tempo medio di collisione
FILE *f_collision=fopen("data/mean_time_collision.dat","a");
fprintf(f_collision,"%e\t%e\n",fraz_imp,2*numOfCollisions/(double)NUMBER_OF_PARTICLES*(total_time-tempoterm));
FILE *f_pression=fopen("data/pression-eta.dat","a");
fprintf(f_pression,"%e\t%e\n",fraz_imp, pression);
//print_coordinate();
//print_speed();
//histo("data/vx.dat",50);
//histo("data/time_collision.dat",50);
free(particleList);
free(collTable);
free(List_0);
exit(EXIT_SUCCESS);
}
