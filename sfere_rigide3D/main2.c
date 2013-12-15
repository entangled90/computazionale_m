#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "gnuplot_i.h"

/*Numero di dimensioni */
#define N 3
#define TERM_TIME 30000
#define MAX_COLLISION 8e4
/*Numero particelle */
unsigned int NUMBER_OF_PARTICLES = 250;
/* Diametro sfere */
double SIGMA =  0;
/*Tavola delle collisioni */
double * collTable;
/*  (i,j,tempo di collisione) */
unsigned index_collision[2];
double time_collision = 0;
unsigned int numOfCollisions = 0;
double total_time = 0;
double temperature = 0;
double K_BOLTZ=1;
double pression = 0;
double DIST_RET = 0;
double deltaV_pre[N];
double deltaV_post[N];
typedef struct particle_s {
	double position[N];
	double speed[N];
	} particle_s ;
particle_s * particleList;

/************
 * 
 * 
 * DEVI FARE DELTA R^2
 * DEVI FARE DELTA R^2
 * DEVI FARE DELTA R^2
 * DEVI FARE DELTA R^2DEVI FARE DELTA R^2DIST_RETDEVI FARE DELTA R^2
 * DEVI FARE DELTA R^2
 * DEVI FARE DELTA R^2
 * DEVI FARE DELTA R^2DEVI FARE DELTA R^2DEVI FARE DELTA R^2DEVI FARE DELTA R^2DIST_RET
 * DEVI FARE DELTA R^2
 * 
 * 
 * 
 * DEVI FARE DELTA R^2DIST_RETDEVI FARE DELTA R^2
 * DEVI FARE DELTA R^2
 * 
 *******************/






void print_coordinate (){
	FILE *f = fopen ( "data/pack.dat","w");
	int i = 0;
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		fprintf(f,"%e \t %e\t%e\n", particleList[i].position[0], particleList[i].position[1], particleList[i].position[2]);
	}
	fclose(f);
	}


void print_speed (){
	FILE *f = fopen ( "data/speed.dat","w");
	int i = 0;
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		fprintf(f,"%e \t %e\t%e\n", particleList[i].speed[0], particleList[i].speed[1],particleList[i].speed[2]);
	}
	fclose(f);
	}




void particle_init ( particle_s *particleList ){
	int i_part= 0;
	int i,j;
	double r_0[N];
	double speed_cm[N];
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
			particleList[i].speed[j] = 2*(rand()/((double)RAND_MAX)) - 1.0 ;
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
		sum += v1[i]*v2[i];
	}
	return (sum);
}
inline void scalar_mult  ( double scalar , double* vec){
	int i = 0;
	for( i = 0 ; i< N ; i++){
		vec[i] *= scalar;
	}
}

void check_distance (){
	int i,j;
	double distance = 0;
	double diff_v[N];
	int x,y,z;
	particle_s temp_part;
	for (i = 0 ; i< NUMBER_OF_PARTICLES ; i++){
		for(j = i+1;j <NUMBER_OF_PARTICLES ; j++){
			for ( x= -1; x < 2 ; x++){
				for ( y = -1; y<2 ; y++){
					for ( z= -1; z<2;z++){	
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
// Calcola il tempo minimo fra le 9 immagini 
inline double calc_min ( int i , int j){
	double x,y,z;
	double min= DBL_MAX;
	double r_diff[N];
	double v_diff[N];
	double det;
	double temp;
	double rv,vv;
	particle_s temp_part;
	for ( x= -1; x < 2 ; x++){
		for ( y = -1; y<2 ; y++){
			for ( z= -1; z<2;z++){
				temp_part = particleList[j];
				temp_part.position[0] += x;
				temp_part.position[1] += y;
				temp_part.position[2] += z;
				diff(particleList[i].position,temp_part.position, r_diff);
				diff( particleList[i].speed,temp_part.speed, v_diff);
				rv = scalar_prod( r_diff, v_diff);
				if( rv < 0){
					vv = scalar_prod(v_diff,v_diff);
					det = rv*rv - vv*( scalar_prod(r_diff,r_diff) -SIGMA*SIGMA);
					if (det > 0){
						temp = ( - rv - sqrt( det ))/ vv ;
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

/* Riempie la matrice dei tempi delle collisioni per j>i */
void collision_table (){
	int i,j;
	for (i = 0; i < NUMBER_OF_PARTICLES ; i++){
		for ( j= i+1 ; j < NUMBER_OF_PARTICLES; j++){
			collTable[i*NUMBER_OF_PARTICLES + j] = calc_min(i,j);
		}
	}
}

void  search_min_coll (){
	int i,j;
	time_collision = DBL_MAX;
	for (i = 0; i < NUMBER_OF_PARTICLES ; i++){
		for ( j= i+1 ; j < NUMBER_OF_PARTICLES; j++){
			if (collTable[i*NUMBER_OF_PARTICLES+j] < time_collision){
				time_collision = collTable[i*NUMBER_OF_PARTICLES+j];
				index_collision[0] = i;
				index_collision[1] = j;
			}
		}
	}
}

void substract_t0 (){
int i,j;
	for (i = 0 ; i < NUMBER_OF_PARTICLES ; i++){
		for ( j = i+1 ; j<NUMBER_OF_PARTICLES ; j++){
			collTable[i*NUMBER_OF_PARTICLES+ j] -= time_collision;
		}
	}
}


void step (){
	int i,j;
	for ( i = 0; i < NUMBER_OF_PARTICLES ; i++){
		for (j =0 ; j< N ;j ++){
			particleList[i].position[j] += time_collision*particleList[i].speed[j];
		}
	}
}
void switch_speeds(){
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
			for ( z= -1; z<2;z++){			
				temp_part = particleList[index_collision[1]];
				temp_part.position[0] += x;
				temp_part.position[1] += y;
				temp_part.position[2] += z;
				diff(particleList[index_collision[0]].position,temp_part.position, r_diff);
				tmp_dbl = scalar_prod(r_diff,r_diff) ; 
				if ( tmp_dbl < min){
					min = tmp_dbl;
					for ( j= 0; j<N; j++){
						temp_r_diff[j] = r_diff[j];
					}
				}
			}
		}
	}
	diff( particleList[index_collision[0]].speed, particleList[index_collision[1]].speed, v_diff);
	scalar_mult( 1/(sqrt(scalar_prod(temp_r_diff,temp_r_diff))), temp_r_diff);
	v_temp = scalar_prod(v_diff,temp_r_diff);
	for ( j = 0 ; j < N ; j++){
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
 *
 * */
void update_coll_table(){
	int i,j;
	int a,b,c;
	for (i = 0; i < 2 ; i++){
		/* a,b indici di riga e colonna -> Matrice simmetrica: tengo solo b>a, ossia j> index_collision[i] */
		for ( j= 0 ; j < NUMBER_OF_PARTICLES; j++){
			a=index_collision[i];
			b=j;
			if( a != b){
				if( a > b){
					c =a;
					a=b;
					b=c;
				}
			}
			collTable[a*NUMBER_OF_PARTICLES+b]= calc_min(a,b);
		}
	}
}


void fix_boundaries (){
	int i = 0;
	int j = 0;
	for (i = 0 ; i< NUMBER_OF_PARTICLES ; i++){
		for( j= 0; j< N ; j++){
			particleList[i].position[j] -= floor(particleList[i].position[j]);
		}
	}
}

double kin_en ( void) {
	int i = 0;
	double  sum = 0;
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		sum += scalar_prod(particleList[i].speed, particleList[i].speed);
		if ( scalar_prod(particleList[i].speed, particleList[i].speed) < 0){
			printf("Vx = %e Vy = %e V^2 = %e\n",particleList[i].speed[0],particleList[i].speed[1],scalar_prod(particleList[i].speed, particleList[i].speed) );
		}
	}
	return sum;
	}
double total_momentum (){
	int i,j;
	double  sum = 0;
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		for ( j = 0; j< N ; j++){
		sum += particleList[i].speed[j];
		}
	}
	return sum;
	}
	
inline void evolve ( ) {
	double deltaV[N];
	/* Calcola la matrice dei tempi*/
	/* Trova la prima coppia che collide  e salva in index_collision e time_collision */
	search_min_coll();
	//printf("Tempo minimo: %e\n",time_collision);
	/* Muove le palline */
	step();
	fix_boundaries();
	diff(particleList[index_collision[0]].speed,particleList[index_collision[1]].speed,deltaV_pre);
	/*Modifica le velocità*/
	switch_speeds();
	total_time += time_collision;
	// Condizioni periodiche
	diff(particleList[index_collision[0]].speed,particleList[index_collision[1]].speed,deltaV_post);
	diff(deltaV_pre,deltaV_post,deltaV);
	substract_t0();
	update_coll_table();
	numOfCollisions +=1;
	pression+=sqrt(scalar_prod(deltaV,deltaV));
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
  /* EX funzione binning */
	for(i=0;i<n;i++){
		freq[i]=0;
	}
	FILE * f;
/*
	f=fopen(input,"r");
	width=(max-min)/(double)n;
	while( fscanf(f,"%e\n",&tmp) == 1){
		if ( tmp > max)
			max= tmp;
		if ( tmp < min)
			min=tmp;
	}
	fclose(f);
*/
	max=3;
	min=0;
	f=fopen(input,"r");
	width=(max-min)/(double)n;
	while( fscanf(f,"%e\n",&tmp) == 1){
		for(i=0;i<n;i++){
		  if( (tmp>min+i*width) && (tmp<=min+(i+1)*width) ){
		 	 freq[i]++;
		  }
		}
	}
	fclose(f);
	/* Necessario per la creazione del file che gnuplot può fittare, diviso negli opportuni intervalli*/
	f=fopen("data/boltzmann-histo.dat","w");
	for(i=0;i<n;i++){
//      tmp = sqrt( (double) freq[i]);
  //    if (tmp == 0)
	//tmp =1;
        fprintf(f,"%lf\t%d\t\n",min+(i+0.5)*width,freq[i]);
    }
    fclose(f);
    free(freq);
}

void boltzmann_file_save ( ){
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
void print_coll_table (){
	int i,j;
	for(i = 0; i<NUMBER_OF_PARTICLES; i++){
		for(j=i+1; j<NUMBER_OF_PARTICLES; j++){
		printf("Tempo collisione (%d,%d): %e\n", i,j,collTable[i*NUMBER_OF_PARTICLES+j]);
		}
	}
	}
/**************MAIN**************/
/*
 */
int main (int argc, char *argv[]) {
/*******************
** INIT VARI
*******************/
gnuplot_ctrl *h;
h=gnuplot_init();

system("rm data/v2.dat");
system("rm data/vz.dat");
system("rm data/vy.dat");
system("rm data/vx.dat");
//Frazione di impacchettamento
double fraz_imp=0.1;
if (argc == 2){
	fraz_imp = atof(argv[1]);
}
srand(time(NULL));
SIGMA = cbrt(6*fraz_imp/ NUMBER_OF_PARTICLES / M_PI);
DIST_RET = cbrt(6*0.61/ (NUMBER_OF_PARTICLES *M_PI));
printf("SIGMA = %e\n",SIGMA);
printf("Frazione di impacchettamento: %e\n", fraz_imp);
collTable = malloc (NUMBER_OF_PARTICLES*NUMBER_OF_PARTICLES*sizeof(double));
particleList = malloc ( NUMBER_OF_PARTICLES * sizeof(particle_s));
particle_init ( particleList);
fix_boundaries();
check_distance();
print_coordinate();
printf("#Collisions: %d \n", numOfCollisions);
printf(" K = %e \t P= %e \n", kin_en(), total_momentum());
temperature = 2*kin_en()/((double) N)/(double) NUMBER_OF_PARTICLES/K_BOLTZ;
/************, norm = 1*****
EVOLUZIONE
*****************/
collision_table();
while ( numOfCollisions < TERM_TIME){
	evolve();
//	printf("min time: %e\n",time_collision);
}
while (numOfCollisions < MAX_COLLISION){
	evolve();
	gnuplot_cmd(h,"splot '-' with points  ps 1 \n");
	for( i = 0; i<NUMBER_OF_PARTICLES;i++){
		gnuplot_cmd(h,"%e\t%e\t%e\n",particleList[i].position[0],particleList[i].position[1],particleList[i].position[2]);
	}
	gnuplot_cmd(h,"e");
//	check_distance();
//	printf("min time: %e\n",time_collision);
	if( numOfCollisions % 1000 == 0 ){
		printf("#Collisions: %d \n", numOfCollisions);
		boltzmann_file_save();
	}
}
pression*=SIGMA/total_time/2.0/kin_en();
pression+=1.0;
pression*=fraz_imp/0.74;
FILE *f_collision=fopen("data/mean_time_collision.dat","a");
fprintf(f_collision,"%e\t%e\n",fraz_imp,2*numOfCollisions/(double)NUMBER_OF_PARTICLES*total_time);
FILE *f_pression=fopen("data/pression-eta.dat","a");
fprintf(f_pression,"%e\t%e\n",fraz_imp, pression);
histo("data/v2.dat",50);
free(particleList);
free(collTable);
gnuplot_close(h);
fclose(f_collision);
fclose(f_pression);
exit(EXIT_SUCCESS);
}

