/*
* USO DEL PROGRAMMA: ./main eta 
altrimenti eta viene impostato di default a eta = 0.1 (fraz_imp)
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "vect3d.h"
#include "raccolta_dati.h"
/*Numero di dimensioni */
#define N 3
#define TERM_TIME 20000
#define MAX_COLLISION 2e5
#define TIME_MAX 150
/*Numero particelle */
int NUMBER_OF_PARTICLES = 128;
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
double T_D = 1;


/**
**
** Necessarie per controllare stato particelle
*
*/
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
		fprintf(f,"%e \t %e \t%e\n", particleList[i].speed[0], particleList[i].speed[1],particleList[i].speed[2]);
	}
	fclose(f);
	}
	/*Calcola energia cinetica*/
double kin_en ( void) {
	int i = 0;
	double  sum = 0;
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		sum += scalar_prod(particleList[i].speed, particleList[i].speed);
	}
	return sum/2.0;
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
/*******************************************************************************************/

void genera_sottoreticolo(double rx_in, double ry_in,double rz_in,int q,int start, double passo){
	int p = start; 
	int c=0;
	int d=0;
	int e=0;
	double rx;
	double ry;
	double rz;
	rx=rx_in;
	ry=ry_in;
	rz=rz_in;
	while(p<NUMBER_OF_PARTICLES && e<q-1){
		while (p < NUMBER_OF_PARTICLES && d<q-1) { 
			while (p < NUMBER_OF_PARTICLES && c<q-1) {
				particleList[p].distance=0;
				particleList[p].n_collision=0;
				particleList[p].last_time_collision=0;
				particleList[p].speed[0] =2*(rand()/(RAND_MAX*1.0)) -1.0 ;
				particleList[p].speed[1] = 2*(rand()/(RAND_MAX*1.0)) -1.0 ;
				particleList[p].speed[2] = 2*(rand()/(RAND_MAX*1.0)) -1.0 ;
				particleList[p].position[0]=rx;
				particleList[p].position[1]=ry;		
				particleList[p].position[2]=rz;
	//			printf("P %d %lf \t %lf\n", p,particleList[p].position[0],particleList[p].position[1]);
				speed_cm[0] += particleList[p].speed[0];
				speed_cm[1] += particleList[p].speed[1];
				speed_cm[2] += particleList[p].speed[2];
				rx = rx + passo;
				p++;
				c++; 
			} 
			rx = rx_in;
			c=0;
			ry = ry + passo;
			d++;
		}
		c=0;
		d=0;
		rz=rz+passo;
		rx=rx_in;
		ry=ry_in;
		e++;
	}
}

//genera un reticolo BCC
void reticolo () { 
	double passo = 0.0;//passo del reticolo 
	//contatori 
	int q=0;
    int m=0;
    int i,j;
    double speed_cm[3]={0.0,0.0,0.0};
     //Definisco il passo del reticolo cercando il minimo doppio di un cubo: m >= n.
      //Questa procedur  	a permette di sfruttare l'intero spazio a disposizione per la creazione del reticolo.
     for (q = 0; m < NUMBER_OF_PARTICLES; q++){
    	m = 2*q*q*q;
    }
    passo = cbrt(2/(double)(m));
	printf("passo %lf\n", passo);
      	  //creazione reticolo

	printf("Primo reticolo\n");
  	genera_sottoreticolo(0,0,0,q,0,passo);

	printf("Secondo reticolo\n");
  	genera_sottoreticolo(passo/2.0,passo/2.0,passo/2.0,q,NUMBER_OF_PARTICLES/2, passo);
	for (i =0 ; i< NUMBER_OF_PARTICLES; i++){
		for ( j = 0; j<N;j++){	
				speed_cm[j] += particleList[i].speed[j];
		}
	}

	for (i =0 ; i< NUMBER_OF_PARTICLES; i++){
		for ( j = 0; j<N;j++){
			particleList[i].speed[j] -= speed_cm[j]/((double) NUMBER_OF_PARTICLES);
		}
	}

	print_coordinate();
	print_speed();

}
void reset_mfp_variables(){
	int i;
	for(i = 0; i<NUMBER_OF_PARTICLES;i++){
			particleList[i].distance=0;
			particleList[i].n_collision=0;
			particleList[i].last_time_collision=0;
	
	}
}

/* Controlla che le sfere non si compenetrino.
*Utilizzata solo all'inizio
*/
int check_distance (){
	int i,j;
	double distance = 0;
	double diff_v[N];
	int x,y,z;
	particle_s temp_part;
	for (i = 0 ; i< NUMBER_OF_PARTICLES ; i++){
		for(j = i+1;j <NUMBER_OF_PARTICLES ; j++){
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
							printf("Sfere (%d,%d) troppo vicine: distanza %lf\n",i,j,distance);
							return 1;
						}
					}
				}
			}			
		}
	}
	return 0;
}
/* Calcola il tempo minimo fra le 27 immagini  */
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
	for (i = 0; i < NUMBER_OF_PARTICLES ; i++){
		for ( j= i+1 ; j < NUMBER_OF_PARTICLES; j++){
			collTable[i*NUMBER_OF_PARTICLES + j] = calc_min ( i, j );
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
/* Sottrae il tempo dell'avvenuta collisione a tutta la matrice (parte superiore dx) */
void substract_t0 (){
int i,j;
	for (i = 0 ; i < NUMBER_OF_PARTICLES ; i++){
		for ( j = i+1 ; j<NUMBER_OF_PARTICLES ; j++){
			collTable[i*NUMBER_OF_PARTICLES+ j] -= time_collision;
		}
	}
}

/* Muove le particelle di uno step temporale*/
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
		for ( j= 0 ; j < NUMBER_OF_PARTICLES; j++){
			a=index_collision[i];
			b=j;
			if( a != b){
				/*considero sempre solo la parte in alto a dx della matrice*/
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

/*Rimette le particelle nella scatola*/
void fix_boundaries (){
	int i = 0;
	int j = 0;
	for (i = 0 ; i< NUMBER_OF_PARTICLES ; i++){
		for( j= 0; j< N ; j++){
			particleList[i].position[j] -= floor(particleList[i].position[j]);
		}
	}
}


/*Calcola momento totale (in norma)*/
double total_momentum (){
	int i,j;
	double  sum[N] = {0,0,0};
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
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
	//Calcola la prossima coppia che si scontra e mette il tempo di collisione in time_collision
	search_min_coll();
	//Mfp da calcolare prima che si siano scambiate le velocità
	mean_free_path();
	//
	/*
	Ossia:
	time_prec è l'ultimo tempo in cui si son salvati i dati
	total_time è il tempo corrente 
	DeltaT è la larghezza di step temporale a cui si vuole calcolare dr2
	if ( total_time + time_collision <time_prec+DeltaT ){
	*/
	// Se non ha superato lo step temporale, muovi sempre prendere dati
	if ( total_time + time_collision - DeltaT -time_prec < 0){
		step(time_collision);
	}
	/* Supererebbe lo step:
	* ~ muovi del tempo necessario per arrivare allo step
	* ~ prende dati
	* ~ muove del tempo necessario per arrivare a time_collision
	*/
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
	particleList[index_collision[0]].last_time_collision=total_time;
	particleList[index_collision[1]].last_time_collision=total_time;
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
	}

void vel_file_save ( ){
	int i = 0;
	FILE *f = fopen("data/v2.dat","a");
	FILE *fx = fopen ("data/vx.dat","a");
	FILE *fy = fopen("data/vy.dat","a");
	FILE *fz = fopen("data/vz.dat","a");
	for (i = 0; i< NUMBER_OF_PARTICLES ; i++){
		fprintf(fx,"%e\n",particleList[i].speed[0]);
		fprintf(fy,"%e\n",particleList[i].speed[1]);
		fprintf(fz,"%e\n",particleList[i].speed[2]);
		fprintf(f,"%e\n",sqrt(scalar_prod(particleList[i].speed,particleList[i].speed)));
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
	

inline void copyList ( particle_s * in , particle_s * out){
	unsigned int i;
	for ( i = 0; i< NUMBER_OF_PARTICLES;i++){
		out[i] = in[i];
	}
}

/*Calcola il minimo di dr2 fra tutte le immagini*/
inline double r_squared_calc ( particle_s * list_0, particle_s * list_1){
	unsigned int i,k;
	double sum = 0;
	double rdiff[N];
	double distance, min;
	double rdiff2[N]={0,0,0};
	int x,y,z;
	particle_s temp_part;
	for ( i = 0; i< NUMBER_OF_PARTICLES;i++){
		min = DBL_MAX;
		for ( x= -1; x < 2 ; x++){
			for ( y = -1; y<2 ; y++){
				for ( z = -1 ; z<2 ; z++){
					temp_part = list_0[i];
					temp_part.position[0] += x;
					temp_part.position[1] += y;
					temp_part.position[2] += z;
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
		}
		sum += scalar_prod(rdiff2,rdiff2);
	}
	return sum/(double)NUMBER_OF_PARTICLES;
} 

/* Fa una media sui tempi dei dr2(delta) per tutti i delta e per tempi tali che sono distanti delta tra di loro */
void r_squared_save ( char * filename){
	FILE *f = fopen(filename, "w");
	double sum=0;
	unsigned int delta,init;
	unsigned int count ;
//	fprintf(f,"%s",header_file);
	double tmp;
	double var =0;
	for ( delta = 1; delta  <  time_counted-1; delta++){
		sum = 0;
		count = 0;
		for ( init = 0; init+delta<time_counted; init++){
			tmp=r_squared_calc( time_list+(init+delta)*NUMBER_OF_PARTICLES,time_list + init*NUMBER_OF_PARTICLES);
			sum += tmp;
			var +=tmp*tmp;
			count++;
		}
		sum /= (double) count;
		var /= (double) count;
		var -= sum*sum;
		fprintf(f,"%.14e\t%.14e\t%.14e\n",delta*DeltaT, sum, sqrt(var/(double)count));
	}
	fclose(f);
}

////// CANNATA é per 2D -> inizializza a T*(2/3)
inline void riscala_vel_temp (){
	int i,j;
	double k_en = kin_en();
	for ( i = 0; i<NUMBER_OF_PARTICLES;i++){
		for (j = 0; j<N;j++){
			particleList[i].speed[j] *= sqrt( NUMBER_OF_PARTICLES* T_D/k_en);
		}
	}
}


/****************************************************************************************
*****************************************************************************************
MAIN
*****************************************************************************************
*****************************************************************************************
*/






int main (int argc, char *argv[]) {
/*******************
** INIT VARI
*******************/
/*Calcola il numero di istanti temporali che verranno salvati*/
NUM_TEMPI_SALVATI = (int) (floor( (double) TIME_MAX / DeltaT)+1);
//unsigned int i ;
srand(time(NULL));
//double dist_tot=0;
double fraz_imp=0.1;

if (argc > 1){
	fraz_imp = atof(argv[1]);
}
SIGMA = cbrt(6*fraz_imp/ NUMBER_OF_PARTICLES / M_PI);
//DIST_RET = sqrt(4*0.76/ NUMBER_OF_PARTICLES / M_PI);
printf("\n\n*****************************************************\n");
printf("Starting simulation with:");
printf("SIGMA = %e\t",SIGMA);
printf("Frazione di impacchettamento: %e\n", fraz_imp);
collTable = malloc (NUMBER_OF_PARTICLES*NUMBER_OF_PARTICLES*sizeof(double));
particleList = malloc ( NUMBER_OF_PARTICLES * sizeof(particle_s));
time_list = malloc (NUM_TEMPI_SALVATI*NUMBER_OF_PARTICLES * sizeof(particle_s));
//particle_init ( particleList);
reticolo();
fix_boundaries();
temperature = 2.0*kin_en()/((double) N)/(double) NUMBER_OF_PARTICLES/K_BOLTZ;
if ( check_distance() != 0){
	printf("Sfere troppo vicine tra loro. Avvio annullato\n");
	exit(EXIT_FAILURE);
}
printf(" K = %e \t P= %e \t", kin_en(), total_momentum());
printf("Temperature is: %f \n",temperature );
riscala_vel_temp();

/****** GESTIONE FILE  ******/
char r2_file[64] = "";
snprintf(r2_file,64,"data/dr2/%d/%.6lf.dat",NUMBER_OF_PARTICLES, fraz_imp); 
//char * press_file = "data/press.dat";
char  pression_filename[128] = "";
snprintf(pression_filename,128,"data/pression/%d/pression%.6lf.dat",NUMBER_OF_PARTICLES,fraz_imp);


print_coordinate();
printf("#Collisions: %d \n", numOfCollisions);

//FILE *pdf_time_coll_fileindex = fopen(tcpdf_file,"w");
/*****************
EVOLUZIONE
*****************/
collision_table();
while ( numOfCollisions < TERM_TIME){
	evolve_therm();
}
numOfCollisions=0;
reset_mfp_variables();
boltzmann_file_save();
total_time = 0;
pression=0;

printf("Termalizzato: %d urti ---- kin_en = %lf\n",numOfCollisions,kin_en());
while (total_time < TIME_MAX){
	evolve();
}
if (time_counted > NUM_TEMPI_SALVATI){
	printf("ERROR \n");
}
r_squared_save(r2_file);
pression /= (double) (2*(total_time)*kin_en());
pression *= SIGMA;
pression +=1.0;
pression *= (fraz_imp/0.7405); //dovuto a PV_0/NKT

FILE * file_pression = fopen(pression_filename,"a");
fprintf(file_pression,"%e\n",pression);
fclose(file_pression);

free(particleList);
free(collTable);
free(time_list);
exit(EXIT_SUCCESS);
}
