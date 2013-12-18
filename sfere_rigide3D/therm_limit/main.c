/*
* USO DEL PROGRAMMA: ./main eta 
altrimenti eta viene impostato di default a eta = 0.1 (fraz_imp)
*/



/****************************
VERSIONE DEL PROGRAMMA PER LIMITE TERMODINAMICO!!!
TUTTTO IL RESTO E' ELIMINATO

**********************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "vect3d.h"
/*Numero di dimensioni */
#define N 3
#define TERM_TIME 20000
#define MAX_COLLISION 2e5
#define TIME_MAX 30
/*Numero particelle da tenere fissato a 256 */
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
	
	while (i_part <number_of_particles){
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
	for ( i = 0; i< number_of_particles; i++){
		for ( j = 0; j<N;j++){
			particleList[i].speed[j] = 2*(rand()/(RAND_MAX*1.0)) - 1.0 ;
			speed_cm[j] += particleList[i].speed[j];
			}
	}
	for (i =0 ; i< number_of_particles; i++){
		for ( j = 0; j<N;j++){
			particleList[i].speed[j] -= (speed_cm[j]/((double) number_of_particles));
		}
	}
}
/* Controlla che le sfere non si compenetrino.
*Utilizzata solo all'inizio
*/
void check_distance (){
	int i,j;
	double distance = 0;
	double diff_v[N];
	int x,y,z;
	particle_s temp_part;
	for (i = 0 ; i< number_of_particles ; i++){
		for(j = i+1;j <number_of_particles ; j++){
			for ( x= -1; x < 2 ; x++){
				for ( y = -1; y<2 ; y++){
					for ( z = -1 ; z<2 ; z++){
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
/* Calcola il tempo minimo fra le 9 immagini  */
double calc_min ( int i , int j){
	double x,y,z;
	double min= DBL_MAX;
	double r_diff[N];
	double v_diff[N];
	double det;
	double temp;
	particle_s temp_part;
	for ( x= -1; x < 2 ; x++){
		for ( y = -1; y<2 ; y++){
			for ( z = -1 ; z<2 ; z++){
				temp_part = particleList[j];
				temp_part.position[0] += x;
				temp_part.position[1] += y;
				temp_part.position[2] += z;
				diff(particleList[i].position,temp_part.position, r_diff);
				diff( particleList[i].speed,temp_part.speed, v_diff);
				if(	 scalar_prod( r_diff, v_diff) < 0){
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
	int  x,y,z;
	double temp_r_diff[N]; /* Vettore differenza temporaneo per le 9 immagini*/
	double v_diff[N];
	double rdiff[N]={0,0,0}; /*Vero vettore differenza*/
	/* temp_r_diff = R0 - R1
	 * v_diff = V0 _ V1
	 */
	double min = DBL_MAX;
	double tmp_dbl;
	particle_s temp_part;
	double v_temp;
	for ( x= -1; x < 2 ; x++){
		for ( y = -1; y<2 ; y++){
			for ( z = -1 ; z<2 ; z++){
				temp_part = particleList[index_collision[1]];
				temp_part.position[0] += x;
				temp_part.position[1] += y;
				temp_part.position[2] += z;
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
	double  sum[N] = {0,0,0};
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
	search_min_coll();
//	mean_free_path();
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
inline void copyList ( particle_s * in , particle_s * out){
	unsigned int i;
	for ( i = 0; i< number_of_particles;i++){
		out[i] = in[i];
	}
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
srand(time(NULL));
double dist_tot=0;
// DA NON CAMBIARE!!!!!!!!!!!
double fraz_imp=0.3;
// ORA L'ARGOMENTO E' IL NUMERO DI PARTICELLE!!!!
if (argc > 1){
	number_of_particles = atof(argv[1]);
}

SIGMA = cbrt(6*fraz_imp/ NUMBER_OF_PARTICLES / M_PI);
DIST_RET = cbrt(6*0.61/ (NUMBER_OF_PARTICLES *M_PI));printf("\n\n*****************************************************\n");
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
//char tcpdf_file[64] = "";
//snprintf(tcpdf_file, 64, "data/pdf_tc/%2f__%s.dat", fraz_imp,date_buffer);
///char * mfp_file = "data/mfp.dat";
snprintf(header_file, 256, "#header: N=%d\t eta=%f\tTIME_MAX=%d\tTERM_TIME=%d\tTEMP=%f\n",N,fraz_imp,TIME_MAX,TERM_TIME,temperature);
/****FINE GESTIONE FILE***/


//check_distance();
print_coordinate();
printf("#Collisions: %d \n", numOfCollisions);

//FILE *pdf_time_coll_fileindex = fopen(tcpdf_file,"w");
/*****************
EVOLUZIONE
*****************/
collision_table();
while ( numOfCollisions < TERM_TIME){
	evolve();
}
total_time = 0;
printf("Termalizzato: %d urti\n",numOfCollisions);
while (total_time < TIME_MAX){
	evolve();
	/*if( numOfCollisions % 10000 == 0 ){
		printf("#Collisions: %d  Total Time: %e\n", numOfCollisions, total_time);
	}
	*/
}
if (time_counted > NUM_TEMPI_SALVATI){
	printf("ERROR \n");
}
pression*=SIGMA/total_time/3.0/kin_en();
pression+=1.0;
pression*=fraz_imp/0.7405;
FILE *f_collision=fopen(tc_file,"a");
fprintf(f_collision,"%e\t%e\n",fraz_imp,total_time/(2*numOfCollisions/(double)number_of_particles));
FILE *f_pression=fopen(press_file,"a");
fprintf(f_pression, "%s\n",header_file);
fprintf(f_pression,"%e\t%e\n\n",fraz_imp, pression);
/*
FILE *f_mean_path = fopen(mfp_file,"w");
for ( i = 0; i< number_of_particles;i++){
	fprintf(f_mean_path,"%e\n",particleList[i].distance/((double)particleList[i].n_collision));
}
FILE *f_mean_mfp = fopen( "data/mfp_eta.dat","a");
for ( i = 0; i<number_of_particles;i++){
	dist_tot += particleList[i].distance;
}
*/
dist_tot /= (double) numOfCollisions;
//fprintf(f_mean_mfp,"%e\t%e\n",fraz_imp, dist_tot);
//fclose(f_mean_mfp);
//fclose(f_mean_path);
fclose(f_collision);
fclose(f_pression);
free(particleList);
free(collTable);
free(time_list);
exit(EXIT_SUCCESS);
}
