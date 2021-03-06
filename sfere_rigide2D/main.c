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
#include "raccolta_dati.h"
/*Numero di dimensioni */
#define N 2
//Numero collisioni di termalizzazione
#define TERM_TIME 10000
#define MAX_COLLISION 2e5
//tempo adimensionato in cui si simula il sistema
#define TIME_MAX 300
/*Numero particelle */
int number_of_particles = 128;
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
double DeltaT= 0.03;
double time_prec;
unsigned int time_counted = 0;
unsigned int NUM_TEMPI_SALVATI;

//Struttura per la particella. COntine velocità, posizione, 
//il tempo in cui ha effettuato l'ultima collisione, il numero di collisioni che ha fatto e la distanza percorsa.
typedef struct particle_s {
	double position[N];
	double speed[N];
	double last_time_collision;
	unsigned int n_collision;
	double distance;
	} particle_s ;
particle_s * particleList;
//Necessaria per il calcolo di Delta r^2 per memorizzare la "storia" di tutta la simulazione
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
	/*Calcola energia cinetica*/
double kin_en ( void) {
	int i = 0;
	double  sum = 0;
	for ( i = 0; i< number_of_particles ; i++){
		sum += scalar_prod(particleList[i].speed, particleList[i].speed);
	}
	return sum/2.0;
	}
	/*Salva il modulo della velocità in un file per poter fare un istogramma */
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

/* Inizializzazione delle particelle */
//Genera un reticolo quadrato a partire dal punto rx,ry con passo p e q particelle per lato
void genera_sottoreticolo(double rx_in, double ry_in, int q,int start, double passo){
	int p = start; 
	int c=0;
	int d=0;
	double rx;
	double ry;
	rx=rx_in;
	ry=ry_in;
	while (p < number_of_particles && d<q-1) { 
		while (p < number_of_particles && c<q-1) {
			particleList[p].distance=0;
			particleList[p].n_collision=0;
			particleList[p].last_time_collision=0;
			particleList[p].speed[0] =2*(rand()/(RAND_MAX*1.0)) -1.0 ;
			particleList[p].speed[1] = 2*(rand()/(RAND_MAX*1.0)) -1.0 ;
			particleList[p].position[0]=rx;
			particleList[p].position[1]=ry;
//			printf("P %d %lf \t %lf\n", p,particleList[p].position[0],particleList[p].position[1]);
			rx = rx + passo;
			p++;
			c++; 
		} 
		rx = rx_in;
		c=0;
		ry = ry + passo;
		d++;
	} 

}
void reset_mfp_variables(){
	int i;
	for(i = 0; i<number_of_particles;i++){
			particleList[i].distance=0;
			particleList[i].n_collision=0;
			particleList[i].last_time_collision=0;
	
	}
}
void reticolo () { 
	double passo = 0.0;//passo del reticolo 
	//contatori 
	double rx = 0.0;
	double ry = 0.0; 
	int q=0;
    int m=0;
    int i,j;
    double speed_cm[2]={0.0,0.0};
     //Definisco il passo del reticolo cercando il minimo doppio di un quadrato: m >= n.
      //Questa procedur  	a permette di sfruttare l'intero spazio a disposizione per la creazione del reticolo.
     for (q = 0; m < number_of_particles; q++){
    	m = 2*q*q;
    }
    passo = sqrt(2/(double)m);
	printf("passo %lf\n", passo);
      	  //creazione reticolo
	rx=0;
	ry=0;
    genera_sottoreticolo(rx,ry,q,0,passo);
	rx = passo/2.0;
	ry = passo/2.0;
    genera_sottoreticolo(rx,ry,q,number_of_particles/2, passo);
	for (i =0 ; i< number_of_particles; i++){
		for ( j = 0; j<N;j++){	
				speed_cm[j] += particleList[i].speed[j];
		}
	}

	for ( i= 0; i<number_of_particles;i++){
		for(j=0;j<N;j++){
			particleList[i].speed[j] -= speed_cm[j]/((double) number_of_particles);
		}
	}
	print_coordinate();
	print_speed();

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


/*Riscala la velocità in modo da avere la temperatura desiderata T_D*/
inline void riscala_vel_temp (){
	int i,j;
	double k_en = kin_en();
	for ( i = 0; i<number_of_particles;i++){
		for (j = 0; j<N;j++){
			particleList[i].speed[j] *= sqrt( number_of_particles* T_D/k_en);
		}
	}
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

/*Modifica le velocità delle particelle coinvolte nella collisione*/
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
* Volendo calcolar dr2(t) l'evoluzioneva di step din step e tiene conto del fatto se facendo lo step temporale che porta alla collisione successiva si supera lo step temporale
* fissato per il dr2(t)
*/
void evolve ( ) {
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
	particleList[index_collision[0]].last_time_collision=total_time;
	particleList[index_collision[1]].last_time_collision=total_time;
	fix_boundaries();
	substract_t0();
	update_coll_table();
	numOfCollisions ++;
	total_time+=time_collision;
	pression += sqrt(scalar_prod(deltaV,deltaV));
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
	//total_time+=time_collision;
	}

void vel_file_save ( ){
	int i = 0;
	FILE *f = fopen("data/v2.dat","a");
	FILE *fx = fopen ("data/vx.dat","a");
	FILE *fy = fopen("data/vy.dat","a");
	for (i = 0; i< number_of_particles ; i++){
		fprintf(fx,"%e\n",particleList[i].speed[0]);
		fprintf(fy,"%e\n",particleList[i].speed[1]);
		fprintf(f,"%e\n",scalar_prod(particleList[i].speed,particleList[i].speed));
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

/*Calcola il minimo di dr2 fra tutte le immagini
	Viene calcolato per tutte le particelle. Le due liste passate sono le liste di particelle a istanti di tempo diversi
	Deve essere chiamata da r_squared_save
*/
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
	return sum/(double)number_of_particles;
} 

/* Fa una media sui tempi dei dr2(delta) per tutti i delta e per tempi tali che sono distanti delta tra di loro */
void r_squared_save ( char * filename){
	FILE *f = fopen(filename, "w");
	double sum=0;
	unsigned int delta,init;
	unsigned int count ;
	fprintf(f,"%s",header_file);
	double tmp;
	double var =0;
	for ( delta = 1; delta  <  time_counted-1; delta++){
		sum = 0;
		count = 0;
		for ( init = 0; init+delta<time_counted; init++){
			tmp=r_squared_calc( time_list+(init+delta)*number_of_particles,time_list + init*number_of_particles);
			sum += tmp;
			var += tmp*tmp;
			count++;
		}
		sum /= (double) count;
		var /= (double) count;
		var -= sum*sum;
		fprintf(f,"%.14e\t%.14e\t%.14e\n",delta*DeltaT, sum, sqrt(var/(double)count));
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
/*******************
** INIT VARI
*******************/
/*Calcola il numero di istanti temporali che verranno salvati*/
NUM_TEMPI_SALVATI = (int) (floor( (double) TIME_MAX / DeltaT)+1);
unsigned int i ;
srand(time(NULL));
double dist_tot=0;
double fraz_imp=0.1;
//double pression_bin[N_BIN_PRESS];
if (argc > 1){
	fraz_imp = atof(argv[1]);
}
SIGMA = sqrt(4*fraz_imp/ number_of_particles / M_PI);
/* DA dove salta fuori?*/
DIST_RET = sqrt(4*0.74/ number_of_particles / M_PI);
printf("\n\n*****************************************************\n");
printf("Starting simulation with:\n");
printf("SIGMA = %e\t",SIGMA);
printf("Frazione di impacchettamento: %e\n", fraz_imp);
collTable = malloc (number_of_particles*number_of_particles*sizeof(double));
particleList = malloc ( number_of_particles * sizeof(particle_s));
time_list = malloc (NUM_TEMPI_SALVATI*number_of_particles * sizeof(particle_s));
reticolo ();
fix_boundaries();

if ( check_distance() != 0){
	printf("Sfere troppo vicine tra loro. Avvio annullato\n");
	exit(EXIT_FAILURE);
}

temperature = 2*kin_en()/((double) N)/(double) number_of_particles/K_BOLTZ;
printf(" K = %e \t P= %e \t", kin_en(), total_momentum());
printf("Temperature is: %f \n",temperature );
riscala_vel_temp();
temperature = 2*kin_en()/((double) N)/(double) number_of_particles/K_BOLTZ;
printf(" K = %e \t P= %e \t", kin_en(), total_momentum());
printf("Temperature is: %f \n",temperature );

/****** GESTIONE FILE  ******/
char r2_file[64] = "";
snprintf(r2_file,64,"data/dr2/dr2_%d_%.6lf.dat",(int)time(NULL),fraz_imp); 
//char * press_file = "data/press.dat";
char  tc_filename[64] = "";
snprintf(tc_filename, 64, "data/tc/%d/tc%6f.dat",number_of_particles, fraz_imp);

char tcpdf_filename[64] = "";
snprintf(tcpdf_filename, 64, "data/pdf_tc/%d/%6f.dat",number_of_particles, fraz_imp);
char  mfp_filename[64] = "";
snprintf(mfp_filename,64,"data/mfp/%d/mfp%.6lf.dat",number_of_particles,fraz_imp);
/****FINE GESTIONE FILE***/
char  pression_filename[128] = "";
snprintf(pression_filename,128,"data/pression/%d/pression%.6lf.dat",number_of_particles,fraz_imp);


print_coordinate();
printf("#Collisions: %d \n", numOfCollisions);

FILE *pdf_tc_file = fopen(tcpdf_filename,"w");
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
	fprintf(pdf_tc_file,"%lf\n",time_collision*number_of_particles/2.0);
}
printf("Num collisioni: %d\n",numOfCollisions);
fclose(pdf_tc_file);
if (time_counted > NUM_TEMPI_SALVATI){
	printf("ERROR \n");
}
r_squared_save(r2_file);
/****** CALCOLO PV/NKT = 1 + 1/(3*number_of_particles*k_boltz*temp)*massa*diametro*Somma collisioni******/
pression /= (double) (3*(total_time)*kin_en());
pression *= SIGMA;
pression +=1.0;
pression *= (fraz_imp/0.9069); //dovuto a PV_0/NKT
FILE *f_collision=fopen(tc_filename,"a");
fprintf(f_collision,"%e\t%e\n",fraz_imp,total_time/(2*numOfCollisions)*(number_of_particles));

FILE * file_pression = fopen(pression_filename,"a");
fprintf(file_pression,"%e\n",pression);
fclose(file_pression);

FILE *f_mean_mfp = fopen( mfp_filename,"a");
for ( i = 0; i<number_of_particles;i++){
	dist_tot += particleList[i].distance/(double)(particleList[i].n_collision); 
}
dist_tot /= (double) (number_of_particles);
fprintf(f_mean_mfp,"%.14e\t%.14e\t\n",fraz_imp, dist_tot);
fclose(f_mean_mfp);
fclose(f_collision);
//fclose(f_pression);
free(particleList);
free(collTable);
free(time_list);
exit(EXIT_SUCCESS);
}
