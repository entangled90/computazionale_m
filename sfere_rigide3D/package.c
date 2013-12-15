#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>

/*Numero di dimensioni */
#define N 3

/*Numero particelle */
int  NUMBER_OF_PARTICLES = 100;
/* Diametro sfere */
double SIGMA =  0;
/*Tavola delle collisioni */
double * collTable;
/*  (i,j,tempo di collisione) */
int index_collision[2];
double time_collision = 0;
int numOfCollisions = 0;
double total_time = 0;
double fraz_imp =0.4;
typedef struct particle_s {
	double position[N];
	double speed[N];
	} particle_s ;
particle_s * particleList;

void particle_init (){
	int i_part= 0;
	int i,j;
	//double r_0= SIGMA/2.0;
	double x_cur = 0;
	double y_cur = 0;
	double z_cur = 0;
	double speed_cm[3];
	int z_row = 0;
	int y_row = 0;
//	Inizializzo il vettore velocità del centro di massa 
	speed_cm[0] = 0;
	speed_cm[1] = 0;
	speed_cm[2] = 0;
	for ( i_part = 0; i_part< NUMBER_OF_PARTICLES ; i_part++){
		x_cur += SIGMA;
		if ( x_cur > 1.0 - SIGMA){
			y_row++;
			x_cur = (double)(y_row%2)*SIGMA/2.0;
			y_cur += sqrt(3)/2.0*SIGMA;
			if ( y_cur > 1.0 - SIGMA){
				y_row = 0;
				y_cur = ((z_row+1)%2)*sqrt(3)*SIGMA/2.0;
				z_row++;
				x_cur = ((z_row)%2)*SIGMA/2.0;
				z_cur += SIGMA/2.0;
			}
		}
		if (z_cur > 1 - SIGMA){
			printf("Impacchettamento non completato\n");
			exit(1);
		}
		particleList[i_part].position[0] = x_cur;
		particleList[i_part].position[1] = y_cur;
		particleList[i_part].position[2] = z_cur;
		//Init velocità random fra [-1,1] */
		particleList[i_part].speed[0] =2* (rand()/(RAND_MAX*1.0)) -1.0 ;
		particleList[i_part].speed[1] = 2*(rand()/(RAND_MAX*1.0)) -1.0 ;
		particleList[i_part].speed[2] = 2*(rand()/(RAND_MAX*1.0)) -1.0 ;
		//Calcolo velocità centro di massa 
		speed_cm[0] += particleList[i_part].speed[0];
		speed_cm[1] += particleList[i_part].speed[1];
		speed_cm[2] += particleList[i_part].speed[2];
		for (i =0 ; i< NUMBER_OF_PARTICLES; i++){
			for ( j = 0; j<N;j++){
			particleList[i].speed[j] -= speed_cm[j]/((double) NUMBER_OF_PARTICLES);
			}
		}
	}
	
	}


void  sum ( double * v1 , double * v2, double * v_sum){
	int i = 0;
	for (i = 0; i < N ; i++){
		v_sum[i] = v1[i] + v2[i];
	}
	}


void diff (double * v1 , double *v2, double * v_diff){
	int i = 0;
	for (i = 0; i < N ; i++){
		v_diff[i] = v1[i] - v2[i];
	}
}

double scalar_prod ( double * v1 , double * v2){
	int i = 0;
	int sum  = 0;
	for ( i = 0; i< N ; i++){
		sum += v1[i]* v2[i];
	}
	return (sum);
	}
void scalar_mult  ( double scalar , double* vec){
	int i = 0;
	for( i = 0 ; i< N ; i++){
		vec[i] *= scalar;
	}
	}

void check_distance (){
	int i,j;
	double distance = 0;
	double diff_v[N];
	for (i = 0 ; i< NUMBER_OF_PARTICLES ; i++){
		for(j = 0;j <NUMBER_OF_PARTICLES ; j++){
		diff(particleList[i].position,particleList[j].position,diff_v);
		if( sqrt(scalar_prod(diff_v,diff_v))<SIGMA){
			printf("Sfere (%d,%d) troppo vicine!\n",i,j);
			}
		}
	}
}

/* Riempie la matrice dei tempi delle collisioni per j>i 
void collision_table (){
	int i,j;
	int x,y;
	double r_diff[N];
	double v_diff[N];
	double det;
	double temp;
	double min ;
	particle_s temp_part;
	for (i = 0; i < NUMBER_OF_PARTICLES ; i++){
		for ( j= i+1 ; j < NUMBER_OF_PARTICLES; j++){
			//Tempo minimo per ognuna della 9 caselle (8 immagini)
			min = DBL_MAX;
			for ( x= -1; x < 2 ; x++){
				for ( y = -1; y<2 ; y++){
					temp_part = particleList[j];
					temp_part.position[0] += x*L;
					temp_part.position[1] += y*L;
					diff(particleList[i].position,temp_part.position, r_diff);
					diff( particleList[i].speed,temp_part.speed, v_diff);
					if( scalar_prod( r_diff, v_diff) < 0){
						det = scalar_prod(r_diff,v_diff)*scalar_prod(r_diff,v_diff) - scalar_prod(v_diff,v_diff)*( scalar_prod(r_diff,r_diff) -SIGMA*SIGMA);
						if (det > 0){
							//printf("Scalar prod: %e \n",(scalar_prod(v_diff,v_diff)));
							temp = ( - scalar_prod( r_diff, v_diff) - sqrt( det ))/ (scalar_prod(v_diff,v_diff) );
							//if ( temp < 0){
								//printf("Tempi negativi: %e\n",temp);
							//}
							if ( temp < min ){
								min = temp;
								//printf("MIN = %e\n",min);
							}
						}
					}
				}
			}
			collTable[i*NUMBER_OF_PARTICLES + j] = min;
			//if (min<0)
				//printf("Tempo minimo per la coppia(%d,%d): %e\n",i,j,min);
		}
	}
}
* */
/* Ritorna un bivettore e un numero: prime due componenti indici di particella, numero è il tempo di collisione 
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
	//printf("Time Collision: %e\n",time_collision);
	}
*/
/*
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
void update_speed_coll (){
	int j;
	double r_diff[N];
	double v_diff[N];
	// r_diff = R0 - R1
	// * v_diff = V0 _ V1
	 
	diff( particleList[index_collision[0]].position, particleList[index_collision[1]].position, r_diff); 
	diff( particleList[index_collision[0]].speed, particleList[index_collision[1]].speed, v_diff);
	scalar_mult( 1/(sqrt(scalar_prod(r_diff,r_diff))), r_diff);
		for ( j = 0 ; j < N ; j++){
			particleList[index_collision[0]].speed[j] -= scalar_prod(v_diff,r_diff) * scalar_prod(particleList[index_collision[0]].speed,r_diff);
			particleList[index_collision[1]].speed[j] += scalar_prod(v_diff,r_diff) * scalar_prod(particleList[index_collision[0]].speed,r_diff);
		}
}
*/
	
/* Aggiorna i tempi di collisioni per le righe e le colonne della matrice della particelle che hanno colliso:
 *
 *
 * NOTA BENE :
 * Modifica le righe associate ad una particella -> LA MATRICE é simmetrica 
 *
 *
 *
void update_coll_table(){
	int i,j;
	double r_diff[N];
	double v_diff[N];
	double det;
	int a,b;
	for (i = 0; i < 2 ; i++){
		// a,b indici di riga e colonna -> Matrice simmetrica: tengo solo b>a, ossia j> index_collision[i] 
		for ( j= 0 ; j < NUMBER_OF_PARTICLES; j++){
			if ( j> index_collision[i] ){
				a = index_collision[i];
				b = j;
				}
			else if ( index_collision[i] > j){
				b = index_collision[i];
				a = j;
			}
			else {
				collTable[index_collision[i]*NUMBER_OF_PARTICLES + j] = 0;
				a=0;
				b=0;
			}
			if ( a != b ){
				diff( particleList[a].position, particleList[b].position, r_diff); 
				diff( particleList[a].speed, particleList[b].speed, v_diff);
				if( scalar_prod( r_diff, v_diff) < 0){
					det =  pow(scalar_prod( r_diff, v_diff),2) - scalar_prod(v_diff,v_diff)*( scalar_prod(r_diff,r_diff) -SIGMA*SIGMA);
					if ( det > 0){
						collTable[a*NUMBER_OF_PARTICLES + b] = ( - scalar_prod( r_diff, v_diff) - sqrt( det ))/ ( pow(scalar_prod(v_diff,v_diff),2) );
					}
				}
			}
		}
	}
}

void fix_boundaries (){
	int i = 0;
	int j = 0;
	for (i = 0 ; i< NUMBER_OF_PARTICLES ; i++){
		for( j= 0; j< N ; j++){
			particleList[i].position[j]= fmod(particleList[i].position[j],1 );
			if ( particleList[i].position[j] < 0 ){
				particleList[i].position[j] += 1;
			}
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


void evolve ( ) {
	// Calcola la matrice dei tempi
	collision_table();
	// Trova la prima coppia che collide  e salva in index_collision e time_collision 
	search_min_coll();
	// Muove le palline 
	step();
	// Modifica le velocità
	update_speed_coll();
	//Riporta la palline nei confini
	fix_boundaries();
	numOfCollisions +=1;
	total_time += time_collision;
	if ( numOfCollisions %100 == 0){
		printf("#Collisions: %d \n", numOfCollisions);
		printf(" K = %e \t P= %e \n", kin_en(), total_momentum());
	}
	}	
*/


void print_coordinate (){
	FILE *f = fopen ( "pack.dat","w");
	int i = 0;
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		fprintf(f,"%e\t%e\t%e\n", particleList[i].position[0], particleList[i].position[1],particleList[i].position[2]);
	}
	fclose(f);
	}


void print_speed (){
	FILE *f = fopen ( "speed.dat","w");
	int i = 0;
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		fprintf(f,"%e\t%e\t%e\n", particleList[i].speed[0], particleList[i].speed[1], particleList[i].speed[2]);
	}
	fclose(f);
	}


void fit( const char * input ){
    double width;
    int* freq;
    int n = 50;
    int i ;
 //   double mean;
  //  double sigma;
    float tmp;
    double max;
    double min ;
    freq =  malloc(n*sizeof(int));
  /* EX funzione binning */
    for(i=0;i<n;i++)
      freq[i]=0;
	max = 2.5;
	min = 0.0;
	FILE * f;
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
	f=fopen("boltzmann-histo.dat","w");
	for(i=0;i<n;i++){
//      tmp = sqrt( (double) freq[i]);
  //    if (tmp == 0)
	//tmp =1;
        fprintf(f,"%lf\t%d\t\n",min+(i+0.5)*width,freq[i]);
    }
    fclose(f);
    free(freq);
}

void boltzmann_file_save ( void ){
	int i = 0;
	double speed_squared = 0;
	FILE *f = fopen ("boltzmann.dat","a");
	for (i = 0; i< NUMBER_OF_PARTICLES ; i++){
		speed_squared = pow(particleList[i].speed[0],2) + pow(particleList[i].speed[1],2);
		fprintf(f,"%e\n",speed_squared);
	}
	fclose(f);
}

int main () {
int maxCollision = 1e4;

SIGMA = cbrt(6*fraz_imp/ ((double)NUMBER_OF_PARTICLES) / M_PI);
printf("SIGMA = %e\n",SIGMA);
printf("Frazione di impacchettamento: %e\n", fraz_imp);
collTable = malloc (NUMBER_OF_PARTICLES*NUMBER_OF_PARTICLES*sizeof(double));
particleList = malloc ( NUMBER_OF_PARTICLES * sizeof(particle_s));
particle_init ();
check_distance();
/*printf("Scalar prod: %e \n", scalar_prod(vec1,vec1));
//print_coordinate();
while ( numOfCollisions < maxCollision){
	evolve();
	if ( numOfCollisions % 200 == 0 ){
		boltzmann_file_save();
	}
}
print_speed();
	fit("boltzmann.dat");


system("rm boltzmann.dat");
*/
print_coordinate();
print_speed();
free(particleList);
free(collTable);
exit(EXIT_SUCCESS);
}
