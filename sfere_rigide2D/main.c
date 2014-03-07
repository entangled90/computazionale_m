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
int number_of_particles = 256;
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


char  header_file[256] = "";
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


/**
**
** Necessarie per controllare stato particelle
*
*/
void print_coordinate (){
	FILE *f = fopen ( "data/pack.dat","w");
	int i = 0;
	for ( i = 0; i< number_of_particles ; i++){
		fprintf(f,"%e \t %e\n", particleList[i].position[0], particleList[i].position[1]);
	}
	fclose(f);
	}


void print_speed (){
	FILE *f = fopen ( "data/speed.dat","w");
	int i = 0;
	for ( i = 0; i< number_of_particles ; i++){
		fprintf(f,"%e \t %e\n", particleList[i].speed[0], particleList[i].speed[1]);
	}
	fclose(f);
	}
	
	
inline void boltzmann_file_save ( void ){
	int i = 0;
	double speed_squared = 0;
	FILE *f = fopen ("data/boltzmann.dat","a");
	for (i = 0; i< number_of_particles ; i++){
		speed_squared = sqrt(scalar_prod(particleList[i].speed,particleList[i].speed));
		fprintf(f,"%e\n",speed_squared);
	}
	fclose(f);
}
/*******************************************************************************************/
// n deve esser maggiore di 0
int min_square ( int n) {
	int i = 1;
	int result = 1;	
	while ( i*i <= n){
		result = i;
		i++;
	}
	return result;
}
/* Inizializzazione delle particelle */
void particle_init (){
	int i_part= 0;
	int i,j;
	int row=0;
	//DIST_RET = 1/((double) min_square(number_of_particles));
	double x_cur = 0;
	double y_cur = DIST_RET/2.0;
	double speed_cm[2];
	speed_cm[0] = 0;
	speed_cm[1] = 0;
	/*if ( DIST_RET < SIGMA){
		printf("Troppe particelle!\n");
		exit(EXIT_FAILURE);
	}
	*/
	for ( i_part = 0; i_part< number_of_particles; i_part++){
		x_cur += DIST_RET;
		particleList[i_part].distance = 0;
		particleList[i_part].n_collision=0;
		particleList[i_part].last_time_collision=0;
		if ( x_cur > 1 - SIGMA ){
			row++;
			x_cur =  (row%2)*DIST_RET/2.0;
			y_cur += sqrt(3)/2.0*DIST_RET;
			
		}
		if ( y_cur > 1){
				printf("Impacchettamento non completato: raggiunto il margine superiore\n");
				exit(EXIT_FAILURE);
		}
			particleList[i_part].position[0] = x_cur;
			particleList[i_part].position[1] = y_cur;
			particleList[i_part].speed[0] =2* (rand()/(RAND_MAX*1.0)) -1.0 ;
			particleList[i_part].speed[1] = 2*(rand()/(RAND_MAX*1.0)) -1.0 ;
			speed_cm[0] += particleList[i_part].speed[0];
			speed_cm[1] += particleList[i_part].speed[1];
			print_coordinate();
		}
	for ( i= 0; i<number_of_particles;i++){
		for(j=0;j<N;j++){
			particleList[i].speed[j] -= speed_cm[j]/((double) number_of_particles);
		}
	}
}


/* Controlla che le sfere non si compenetrino.
*Utilizzata solo all'inizio
EXIT CODE 1 =  ERRORE, si toccano
EXIT CODE 0 =  TUTTO OK
*/
int  check_distance (){
	int i,j;
	double distance = 0;
	double diff_v[N];
	int x,y;
	particle_s temp_part;
	for (i = 0 ; i< number_of_particles ; i++){
		for(j = i+1;j <number_of_particles ; j++){
			for ( x= -1; x < 2 ; x++){
				for ( y = -1; y<2 ; y++){
					temp_part = particleList[j];
					temp_part.position[0] += x;
					temp_part.position[1] += y;
					diff(particleList[i].position,temp_part.position,diff_v);	
					distance = sqrt(scalar_prod(diff_v,diff_v));
					if( distance <SIGMA){
						printf("Sfere (%d,%d) troppo vicine!\n",i,j);
						return (1);
					}
				}
			}			
		}
	}
	return (0);
}
/* Calcola il tempo minimo fra le 9 immagini  */
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
					//uso debug: printf("Scalar prod: %e \n",(scalar_prod(v_diff,v_diff)));
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
	for (i = 0; i < number_of_particles ; i++){
		for ( j= i+1 ; j < number_of_particles; j++){
			collTable[i*number_of_particles + j] = calc_min ( i, j );
		}
	}
}


/*
*	Calcola la prima coppia che colliderà.
Gli indici di particella sono salvati in "index_collision",
mentre il tempo mancante alla collisione in time_collision.
GLOBALI
*/
void search_min_coll (){
	int i,j;
	time_collision = DBL_MAX;
	for (i = 0; i < number_of_particles ; i++){
		for ( j= i+1 ; j < number_of_particles; j++){
			if (collTable[i*number_of_particles+j] < time_collision){
				time_collision = collTable[i*number_of_particles+j];
				index_collision[0] = i;
				index_collision[1] = j;
			}
		}
	}
	}
/* Sottrae il tempo dell'avvenuta collisione a tutta la matrice (parte superiore dx) */
void substract_t0 (){
int i,j;
	for (i = 0 ; i < number_of_particles ; i++){
		for ( j = i+1 ; j<number_of_particles ; j++){
			collTable[i*number_of_particles+ j] -= time_collision;
		}
	}
}

/* Muove le particelle di uno step temporale*/
void step (double time_step){
	int i,j;
	for ( i = 0; i < number_of_particles ; i++){
		for (j =0 ; j< N ;j ++){
			particleList[i].position[j] += time_step*particleList[i].speed[j];
		}
	}
}


void switch_speeds(){
	int j;
	int  x,y;
	double temp_r_diff[N]; /* Vettore differenza temporaneo per le 9 immagini*/
	double v_diff[N];
	double rdiff[2]={0,0}; /*Vero vettore differenza*/
	/* temp_r_diff = R0 - R1
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
			diff(particleList[index_collision[0]].position,temp_part.position, temp_r_diff); /*vettore differenza salvato in temp_r_diff*/
			tmp_dbl = scalar_prod(temp_r_diff,temp_r_diff) ; 
			if ( tmp_dbl < min){
				min = tmp_dbl;
				for ( j= 0; j<N; j++){
				rdiff[j] = temp_r_diff[j];
				}
			}
		}
	}
	diff( particleList[index_collision[0]].speed, particleList[index_collision[1]].speed, v_diff);
	scalar_mult( 1/(sqrt(scalar_prod(rdiff,rdiff))), rdiff);
	v_temp = scalar_prod(v_diff,rdiff);
	for ( j = 0 ; j < N ; j++){
		particleList[index_collision[0]].speed[j] -= v_temp*rdiff[j];
		particleList[index_collision[1]].speed[j] += v_temp*rdiff[j];
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
		for ( j= 0 ; j < number_of_particles; j++){
			a=index_collision[i];
			b=j;
			if( a != b){
				/*considero sempre solo la parte in alto a dx della matrice*/
				if( a>b){
					c=a;
					a=b;
					b=c;
				}
				collTable[a*number_of_particles+b]= calc_min(a,b);
			}
		}
	}
}

/*Rimette le particelle nella scatola*/
void fix_boundaries (){
	int i = 0;
	int j = 0;
	for (i = 0 ; i< number_of_particles ; i++){
		for( j= 0; j< N ; j++){
			particleList[i].position[j] -= floor(particleList[i].position[j]);
		}
	}
}

/*Calcola energia cinetica*/
double kin_en ( void) {
	int i = 0;
	double  sum = 0;
	for ( i = 0; i< number_of_particles ; i++){
		sum += scalar_prod(particleList[i].speed, particleList[i].speed);
		if ( scalar_prod(particleList[i].speed, particleList[i].speed) < 0){
			printf("Vx = %e Vy = %e V^2 = %e\n",particleList[i].speed[0],particleList[i].speed[1],scalar_prod(particleList[i].speed, particleList[i].speed) );
		}
	}
	return sum;
	}

/*Calcola momento totale (in norma)*/
double total_momentum (){
	int i,j;
	double  sum[2] = {0,0};
	for ( i = 0; i< number_of_particles ; i++){
		for ( j = 0; j< N ; j++){
		sum[j] += particleList[i].speed[j];
		}
	}
	return sqrt(scalar_prod(sum,sum));
	}
	
/*Calcola il libero cammino medio:
*Esso viene calcolato come lo spazio percorso dalla particella a partire dall'ultimo urto che ha fatto.
Questo viene salvato in .last_time_collision
*/
inline void  mean_free_path (){
	unsigned i;
	for ( i = 0; i<2;i++){
		particleList[index_collision[i]].n_collision++;
		particleList[index_collision[i]].distance += (total_time+time_collision-particleList[index_collision[i]].last_time_collision)*sqrt(scalar_prod(particleList[index_collision[i]].speed,particleList[index_collision[i]].speed));
	}
}

/*Evolve il sistema di uno step
* Volendo calcolar dr2(t) l'evoluzione non va di step din step, ma di dt in dt.
*/
inline void evolve ( ) {
	double deltaV_pre[N];
	double deltaV_post[N];
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
		for ( j = 0; j< number_of_particles;j++){
			time_list[time_counted*number_of_particles+j] = particleList[j];
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

/* Evolve ma utilizzata solo in fase di termalizzazione, senza alcuna presa dati*/
inline void evolve_therm ( ) {
	double deltaV_pre[N];
	search_min_coll();
	step(time_collision);
	diff(particleList[index_collision[0]].speed,particleList[index_collision[1]].speed,deltaV_pre);
	switch_speeds();
	//condizioni al bordo
	fix_boundaries();
	substract_t0();
	update_coll_table();
	numOfCollisions +=1;
	total_time+=time_collision;
	}

void vel_file_save ( ){
	int i = 0;
	FILE *f = fopen("data/v2.dat","a");
	FILE *fx = fopen ("data/vx.dat","a");
	FILE *fy = fopen("data/vy.dat","a");
	for (i = 0; i< number_of_particles ; i++){
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
	for(i = 0; i<number_of_particles; i++){
		for(j=i+1; j<number_of_particles; j++){
		printf("Tempo collisione (%d,%d): %e\n", i,j,collTable[i*number_of_particles+j]);
		}
	}
	}
	

inline void copyList ( particle_s * in , particle_s * out){
	unsigned int i;
	for ( i = 0; i< number_of_particles;i++){
		out[i] = in[i];
	}
}

/*Calcola il minimo di dr2 fra tutte le immagini*/
inline double r_squared_calc ( particle_s * list_0, particle_s * list_1){
	unsigned int i,k;
	double sum = 0;
	double rdiff[N];
	double distance, min;
	double rdiff2[2]={0,0};
	int x,y;
	particle_s temp_part;
	for ( i = 0; i< number_of_particles;i++){
		min = DBL_MAX;
		for ( x= -1; x < 2 ; x++){
			for ( y = -1; y<2 ; y++){
				temp_part = list_0[i];
				temp_part.position[0] += x;
				temp_part.position[1] += y;
				diff(list_1[i].position,temp_part.position,rdiff);
				distance = scalar_prod(rdiff,rdiff);
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
	return sum/number_of_particles;
} 

/* Fa una media sui tempi dei dr2(delta) per tutti i delta e per tempi tali che sono distanti delta tra di loro */
void r_squared_save ( char * filename){
	FILE *f = fopen(filename, "w");
	double sum=0;
	unsigned int delta,init;
	unsigned int count ;
	fprintf(f,"%s",header_file);
	for ( delta = 1; delta  <  time_counted-1; delta++){
		sum = 0;
		count = 0;
		for ( init = 0; init+delta<time_counted; init++){
			sum += r_squared_calc( time_list+(init+delta)*number_of_particles,time_list + init*number_of_particles);
			count++;
		}
		sum /= (double) count;
		fprintf(f,"%e\t%e\n",delta*DeltaT, sum);
	}
	fclose(f);
}



/****************************************************************************************
*****************************************************************************************
MAIN
*****************************************************************************************
*****************************************************************************************
*/






int main (int argc, char *argv[]) {
/* acquisisce il tempo per salvare i file con nomi sempre diversi*/
time_t rawtime;
struct tm * timeinfo;
char date_buffer [30];
time (&rawtime);
timeinfo = localtime (&rawtime);
strftime (date_buffer,30,"%F--%T",timeinfo);

/*******************
** INIT VARI
*******************/
/*Calcola il numero di istanti temporali che verranno salvati*/
NUM_TEMPI_SALVATI = (int) (floor( (double) TIME_MAX / DeltaT)+1);
unsigned int i ;
srand(time(NULL));
double dist_tot=0;
double fraz_imp=0.1;

if (argc > 1){
	fraz_imp = atof(argv[1]);
}
SIGMA = sqrt(4*fraz_imp/ number_of_particles / M_PI);
/* DA dove salta fuori?*/
DIST_RET = sqrt(4*0.74/ number_of_particles / M_PI);
printf("\n\n*****************************************************\n");
printf("Starting simulation with:");
printf("SIGMA = %e\t",SIGMA);
printf("Frazione di impacchettamento: %e\n", fraz_imp);
collTable = malloc (number_of_particles*number_of_particles*sizeof(double));
particleList = malloc ( number_of_particles * sizeof(particle_s));
time_list = malloc (NUM_TEMPI_SALVATI*number_of_particles * sizeof(particle_s));
particle_init ( particleList);
fix_boundaries();
temperature = 2*kin_en()/((double) N)/(double) number_of_particles/K_BOLTZ;
printf(" K = %e \t P= %e \t", kin_en(), total_momentum());
printf("Temperature is: %f \n",temperature );

/****** GESTIONE FILE  ******/
char r2_file[64] = "";
snprintf(r2_file,64,"data/dr2/%.2f__%s.dat",fraz_imp,date_buffer); 
char * press_file = "data/press.dat";
char * tc_file = "data/tc.dat";
char tcpdf_file[64] = "";
snprintf(tcpdf_file, 64, "data/pdf_tc/%2f__%s.dat", fraz_imp,date_buffer);
char  mfp_file[64] = "";
snprintf(mfp_file,64,"data/mfp/mfp%.6lf.dat",fraz_imp);
snprintf(header_file, 256, "#header: N=%d\t eta=%f\tTIME_MAX=%d\tTERM_TIME=%d\tTEMP=%f\n",N,fraz_imp,TIME_MAX,TERM_TIME,temperature);
/****FINE GESTIONE FILE***/


if ( check_distance() != 0){
	printf("Sfere troppo vicine tra loro. Avvio annullato\n");
	exit(EXIT_FAILURE);
}
print_coordinate();
printf("#Collisions: %d \n", numOfCollisions);

FILE *pdf_time_coll_fileindex = fopen(tcpdf_file,"w");
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
	/*if( numOfCollisions % 10000 == 0 ){
		printf("#Collisions: %d  Total Time: %e\n", numOfCollisions, total_time);
	}
	*/
	fprintf(pdf_time_coll_fileindex,"%f\n",time_collision);
}
fclose(pdf_time_coll_fileindex);
if (time_counted > NUM_TEMPI_SALVATI){
	printf("ERROR \n");
}
r_squared_save(r2_file);
pression*=SIGMA/total_time/3.0/kin_en();
pression+=1.0;
pression*=fraz_imp/M_PI*2*sqrt(3.00);
pression *= number_of_particles*temperature;
FILE *f_collision=fopen(tc_file,"a");
fprintf(f_collision,"%e\t%e\n",fraz_imp,total_time/(2*numOfCollisions/(double)number_of_particles));
FILE *f_pression=fopen(press_file,"a");
fprintf(f_pression, "%s\n",header_file);
fprintf(f_pression,"%e\t%e\n",fraz_imp, pression);
FILE *f_mean_path = fopen(mfp_file,"w");
for ( i = 0; i< number_of_particles;i++){
	fprintf(f_mean_path,"%e\n",particleList[i].distance/((double)particleList[i].n_collision));
}
FILE *f_mean_mfp = fopen( "data/mfp_eta.dat","a");
for ( i = 0; i<number_of_particles;i++){
	dist_tot += particleList[i].distance;
}
dist_tot /= (double) (numOfCollisions*number_of_particles);
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
