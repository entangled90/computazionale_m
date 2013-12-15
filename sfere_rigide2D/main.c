/*
* USO DEL PROGRAMMA: ./main eta 
altrimenti eta viene impostato di default a eta = 0.1 (fraz_imp)
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "vect2d.h"
/*Numero di dimensioni */
#define N 2
#define TERM_TIME 20000
#define MAX_COLLISION 2e5
#define TIME_MAX 30
/*Numero particelle */
int NUMBER_OF_PARTICLES = 256;
/* Diametro sfere */
double SIGMA =  0;
/*Tavola delle collisioni */
double * collTable;
/*  (i,j,tempo di collisione) */
int index_collision[2];
double time_collision = 0;
int numOfCollisions = 0;
double total_time = 0;
double temperature = 0;
double K_BOLTZ=1;
double pression = 0;
double D_speed_norm = 0;
double DIST_RET = 0;
double deltaV_pre[N];
double deltaV_post[N];
double DeltaT= 0.1;
double time_prec;
unsigned int time_counted = 0;
unsigned int NUM_TEMPI_SALVATI;
typedef struct particle_s {
	double position[N];
	double speed[N];
	double last_time_collision;
	unsigned int n_collision;
	double distance;
	} particle_s ;
particle_s * particleList;
particle_s * time_list;


void print_coordinate (){
	FILE *f = fopen ( "data/pack.dat","w");
	int i = 0;
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		fprintf(f,"%e \t %e\n", particleList[i].position[0], particleList[i].position[1]);
	}
	fclose(f);
	}


void print_speed (){
	FILE *f = fopen ( "data/speed.dat","w");
	int i = 0;
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		fprintf(f,"%e \t %e\n", particleList[i].speed[0], particleList[i].speed[1]);
	}
	fclose(f);
	}
	
	
inline void boltzmann_file_save ( void ){
	int i = 0;
	double speed_squared = 0;
	FILE *f = fopen ("data/boltzmann.dat","a");
	for (i = 0; i< NUMBER_OF_PARTICLES ; i++){
		speed_squared = sqrt(scalar_prod(particleList[i].speed,particleList[i].speed));
		fprintf(f,"%e\n",speed_squared);
	}
	fclose(f);
}



void particle_init (){
	int i_part= 0;
	int i,j;
	int row=0;
	double x_cur = - DIST_RET;
	double y_cur = 0;
	double speed_cm[2];
	speed_cm[0] = 0;
	speed_cm[1] = 0;
	/*if ( DIST_RET < SIGMA){
		printf("Troppe particelle!\n");
		exit(EXIT_FAILURE);
	}
	*/
	for ( i_part = 0; i_part< NUMBER_OF_PARTICLES; i_part++){
		x_cur += DIST_RET;
		particleList[i_part].distance = 0;
		particleList[i_part].n_collision=0;
		particleList[i_part].last_time_collision=0;
		if ( x_cur > 1 - SIGMA){
			row++;
			x_cur =  (row%2)*DIST_RET/2.0;
			y_cur += sqrt(3)/2.0*DIST_RET;
			
		}
		if ( y_cur > 1 - SIGMA){
				printf("Impacchettamento non completato: raggiunto il margine superiore\n");
				exit(EXIT_FAILURE);
		}
			particleList[i_part].position[0] = x_cur;
			particleList[i_part].position[1] = y_cur;
			particleList[i_part].speed[0] =2* (rand()/(RAND_MAX*1.0)) -1.0 ;
			particleList[i_part].speed[1] = 2*(rand()/(RAND_MAX*1.0)) -1.0 ;
			speed_cm[0] += particleList[i_part].speed[0];
			speed_cm[1] += particleList[i_part].speed[1];
			//print_coordinate();
		}
	for ( i= 0; i<NUMBER_OF_PARTICLES;i++){
		for(j=0;j<N;j++){
			particleList[i].speed[j] -= speed_cm[j]/((double) NUMBER_OF_PARTICLES);
		}
	}
}

void check_distance (){
	int i,j;
	double distance = 0;
	double diff_v[N];
	int x,y;
	particle_s temp_part;
	for (i = 0 ; i< NUMBER_OF_PARTICLES ; i++){
		for(j = i+1;j <NUMBER_OF_PARTICLES ; j++){
			for ( x= -1; x < 2 ; x++){
				for ( y = -1; y<2 ; y++){
					temp_part = particleList[j];
					temp_part.position[0] += x;
					temp_part.position[1] += y;
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
// Calcola il tempo minimo fra le 9 immagini 
double calc_min ( int i , int j){
	double x,y;
	double min= DBL_MAX;
	double r_diff[N];
	double v_diff[N];
	double det;
	double temp;
	particle_s temp_part;
	for ( x= -1; x < 2 ; x++){
		for ( y = -1; y<2 ; y++){
			temp_part = particleList[j];
			temp_part.position[0] += x;
			temp_part.position[1] += y;
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
	return min;
}
/* Riempie la matrice dei tempi delle collisioni per j>i */
void collision_table (){
	int i,j;
	for (i = 0; i < NUMBER_OF_PARTICLES ; i++){
		for ( j= i+1 ; j < NUMBER_OF_PARTICLES; j++){
			collTable[i*NUMBER_OF_PARTICLES + j] = calc_min ( i, j );
		}
	}
}


/* Ritorna un bivettore e un numero: prime due componenti indici di particella, numero è il tempo di collisione */
void search_min_coll (){
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


void step (double time_step){
	int i,j;
	for ( i = 0; i < NUMBER_OF_PARTICLES ; i++){
		for (j =0 ; j< N ;j ++){
			particleList[i].position[j] += time_step*particleList[i].speed[j];
		}
	}
}
void switch_speeds(){
	int j;
	int  x,y;
	double r_diff[N];
	double v_diff[N];
	double temp_r_diff[N]={0,0}; // 
	/* r_diff = R0 - R1
	 * v_diff = V0 _ V1
	 */
	double min = DBL_MAX;
	double tmp_dbl;
	particle_s temp_part;
	double v_temp;
	for ( x= -1; x < 2 ; x++){
		for ( y = -1; y<2 ; y++){
			temp_part = particleList[index_collision[1]];
			temp_part.position[0] += x;
			temp_part.position[1] += y;
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
				if( a>b){
					c=a;
					a=b;
					b=c;
				}
				collTable[a*NUMBER_OF_PARTICLES+b]= calc_min(a,b);
			}
		}
	}
}
/*
				diff( particleList[a].position, particleList[b].position, r_diff); 
				diff( particleList[a].speed, particleList[b].speed, v_diff);
				if( scalar_prod( r_diff, v_diff) < 0){
					det =  pow(scalar_prod( r_diff, v_diff),2) - scalar_prod(v_diff,v_diff)*( scalar_prod(r_diff,r_diff) -SIGMA*SIGMA);
					if ( det > 0){
						
						collTable[a*NUMBER_OF_PARTICLES + b] = ( - scalar_prod( r_diff, v_diff) - sqrt( det ))/ ( pow(scalar_prod(v_diff,v_diff),2) );
					}
				}
*/

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
	double  sum[2] = {0,0};
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		for ( j = 0; j< N ; j++){
		sum[j] += particleList[i].speed[j];
		}
	}
	return sqrt(scalar_prod(sum,sum));
	}
	
inline void  mean_free_path (){
	unsigned i;
	for ( i = 0; i<2;i++){
		particleList[index_collision[i]].n_collision++;
		particleList[index_collision[i]].distance += (total_time+time_collision-particleList[index_collision[i]].last_time_collision)*sqrt(scalar_prod(particleList[index_collision[i]].speed,particleList[index_collision[i]].speed));
	}
}

inline void evolve ( ) {
	double deltaV[N];
	unsigned int j = 0;
	search_min_coll();
	mean_free_path();
	if ( total_time + time_collision - DeltaT -time_prec < 0){
		step(time_collision);
	}
	else{
		time_counted++;
		step( time_prec + DeltaT - total_time);
		for ( j = 0; j< NUMBER_OF_PARTICLES;j++){
			time_list[time_counted*NUMBER_OF_PARTICLES+j] = particleList[j];
		}
		step( total_time+ time_collision - time_prec - DeltaT);
		time_prec += DeltaT;
	}
	diff(particleList[index_collision[0]].speed,particleList[index_collision[1]].speed,deltaV_pre);
	switch_speeds();
	//calcoli pressione
	diff(particleList[index_collision[0]].speed,particleList[index_collision[1]].speed,deltaV_post);
	diff(deltaV_pre,deltaV_post,deltaV);
	//condizioni al bordo
	fix_boundaries();
	substract_t0();
	update_coll_table();
	numOfCollisions +=1;
	total_time+=time_collision;
	pression+= sqrt(scalar_prod(deltaV,deltaV));
	}

inline void evolve_therm ( ) {
	double deltaV[N];
	search_min_coll();
	step(time_collision);
	diff(particleList[index_collision[0]].speed,particleList[index_collision[1]].speed,deltaV_pre);
	switch_speeds();
	//calcoli pressione
	diff(particleList[index_collision[0]].speed,particleList[index_collision[1]].speed,deltaV_post);
	diff(deltaV_pre,deltaV_post,deltaV);
	//condizioni al bordo
	fix_boundaries();
	substract_t0();
	update_coll_table();
	numOfCollisions +=1;
	total_time+=time_collision;
	pression+= sqrt(scalar_prod(deltaV,deltaV));
	}

void vel_file_save ( ){
	int i = 0;
	FILE *f = fopen("data/v2.dat","a");
	FILE *fx = fopen ("data/vx.dat","a");
	FILE *fy = fopen("data/vy.dat","a");
	for (i = 0; i< NUMBER_OF_PARTICLES ; i++){
		fprintf(fx,"%e\n",particleList[i].speed[0]);
		fprintf(fy,"%e\n",particleList[i].speed[1]);
		fprintf(f,"%e\n",sqrt(scalar_prod(particleList[i].speed,particleList[i].speed)));
	}
	fclose(fx);
	fclose(fy);
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
	

inline void copyList ( particle_s * in , particle_s * out){
	unsigned int i;
	for ( i = 0; i< NUMBER_OF_PARTICLES;i++){
		out[i] = in[i];
	}
}


inline double r_squared_calc ( particle_s * list_0, particle_s * list_1){
	unsigned int i,k;
	double sum = 0;
	double rdiff[N];
	double distance, min;
	double rdiff2[N];
	int x,y;
	particle_s temp_part;
	for ( i = 0; i< NUMBER_OF_PARTICLES;i++){
		min = DBL_MAX;
		for ( x= -1; x < 2 ; x++){
			for ( y = -1; y<2 ; y++){
				temp_part = list_0[i];
				temp_part.position[0] += x;
				temp_part.position[1] += y;
				diff(list_1[i].position,temp_part.position,rdiff);
				distance = sqrt(scalar_prod(rdiff,rdiff));
				if( distance < min ){
					min = distance;
					for ( k = 0; k<N;k++){
						rdiff2[k] = rdiff[k];
					}
				}
			}
		}
		sum += scalar_prod(rdiff2,rdiff2);
	}
	return sum/NUMBER_OF_PARTICLES;
} 

void r_squared_save ( char * filename){
	FILE *f1 = fopen(filename, "w");
	fclose(f1);
	FILE *f = fopen(filename, "a");
	double sum=0;
	unsigned int delta,init;
	unsigned int count ;
	for ( delta = 1; delta  <  time_counted-1; delta++){
		sum = 0;
		count = 0;
		for ( init = 0; init +delta < time_counted; init++){
			sum += r_squared_calc( time_list+(init+delta)*NUMBER_OF_PARTICLES,time_list + init*NUMBER_OF_PARTICLES);
			count++;
		}
		sum /= (double) count;
//		sum=r_squared_calc(time_list,time_list+delta);
		fprintf(f,"%e\t%e\n",delta*DeltaT, sum);
	}
	fclose(f);
}

int main (int argc, char *argv[]) {
/*******************
** INIT VARI
*******************/
NUM_TEMPI_SALVATI = (int) (floor( (double) TIME_MAX / DeltaT)+1);
unsigned int i ;
srand(time(NULL));
double dist_tot=0;
double fraz_imp=0.1;

if (argc > 1){
	fraz_imp = atof(argv[1]);
}
SIGMA = sqrt(4*fraz_imp/ NUMBER_OF_PARTICLES / M_PI);
DIST_RET = sqrt(4*0.76/ NUMBER_OF_PARTICLES / M_PI);
printf("SIGMA = %e\n",SIGMA);
printf("Frazione di impacchettamento: %e\n", fraz_imp);
collTable = malloc (NUMBER_OF_PARTICLES*NUMBER_OF_PARTICLES*sizeof(double));
particleList = malloc ( NUMBER_OF_PARTICLES * sizeof(particle_s));
time_list = malloc (NUM_TEMPI_SALVATI*NUMBER_OF_PARTICLES * sizeof(particle_s));
particle_init ( particleList);
fix_boundaries();
temperature = 2*kin_en()/((double) N)/(double) NUMBER_OF_PARTICLES/K_BOLTZ;
/****** Variabili per nomi dei file */
char r2_file[64] = "";
snprintf(r2_file,64,"data/dr2/%f_%d.dat",fraz_imp, rand()); //Aggiunto un numero casuale per evitare sovrascritture
char * press_file = "data/press.dat";
char * tc_file = "data/tc.dat";
char tcpdf_file[64] = "";
snprintf(tcpdf_file, 64, "data/pdf_tc/%f.dat", fraz_imp);
char * mfp_file = "data/mfp.dat";
char  header_file[256] = "";
snprintf(header_file, 256, "#header: N=%d\t eta=%f\tTIME_MAX=%d\tTERM_TIME=%d\tTEMP=%f\n",N,fraz_imp,TIME_MAX,TERM_TIME,temperature);
//check_distance();
print_coordinate();
printf("#Collisions: %d \n", numOfCollisions);
printf(" K = %e \t P= %e \n", kin_en(), total_momentum());
temperature = 2*kin_en()/((double) N)/(double) NUMBER_OF_PARTICLES/K_BOLTZ;
//FILE *f_time_collision = fopen("data/time_collision.dat","w");
/*****************
EVOLUZIONE
*****************/
collision_table();
while ( numOfCollisions < TERM_TIME){
	evolve_therm();
}
total_time = 0;
printf("Termalizzato: %d urti\n",numOfCollisions);
while (total_time < TIME_MAX){
	evolve();
	if( numOfCollisions % 10000 == 0 ){
		printf("#Collisions: %d  Total Time: %e\n", numOfCollisions, total_time);
	}
}
if (time_counted > NUM_TEMPI_SALVATI){
	printf("ERROR \n");
}
r_squared_save(r2_file);
//fclose(f_time_collision);
pression*=SIGMA/total_time/3.0/kin_en();
pression+=1.0;
pression*=fraz_imp/M_PI*2*sqrt(3.00);
FILE *f_collision=fopen(tc_file,"a");
fprintf(f_collision,"%e\t%e\n",fraz_imp,2*numOfCollisions/(double)NUMBER_OF_PARTICLES*total_time);
FILE *f_pression=fopen(press_file,"a");
fprintf(f_pression,"%e\t%e\n",fraz_imp, pression);
FILE *f_mean_path = fopen(mfp_file,"w");
for ( i = 0; i< NUMBER_OF_PARTICLES;i++){
	fprintf(f_mean_path,"%e\n",particleList[i].distance/((double)particleList[i].n_collision));
}
FILE *f_mean_mfp = fopen( "data/mfp_eta.dat","a");
for ( i = 0; i<NUMBER_OF_PARTICLES;i++){
	dist_tot += particleList[i].distance;
}
dist_tot /= (double) numOfCollisions;
fprintf(f_mean_mfp,"%e\t%e\n",fraz_imp, dist_tot);
fclose(f_mean_mfp);
fclose(f_mean_path);
fclose(f_collision);
fclose(f_pression);
free(particleList);
free(collTable);
free(time_list);
exit(EXIT_SUCCESS);
}
