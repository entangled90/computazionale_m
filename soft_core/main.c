#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>

#define NUMBER_OF_PARTICLES 108
#define N 3
#define ITERATION_MAX 1e4
#define ITERATION_THERM 10000

double SIGMA=1;
double DIST_RET = 1;
double EPS = 1;
double acceleration[N];
double R_LIM;
double r;
double D_T = 1e-3;
double total_time=0;
double DeltaF;
double u_R_LIM;
double L;
double V_MAX = 1;
double T_D = 1.19;
typedef struct particle_s {
	double position[N];
	double speed[N];
	double acc[N];
	} particle_s ;

particle_s * particleList;
void print_coordinate (){
	FILE *f = fopen ( "data/pack.dat","w");
	int i = 0;
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		fprintf(f,"%e\t%e\t%e\n", particleList[i].position[0], particleList[i].position[1],particleList[i].position[2]);
	}
	fclose(f);
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

inline double kin_en ( void) {
	int i = 0;
	double  sum = 0;
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		sum += 0.5*scalar_prod(particleList[i].speed, particleList[i].speed);
	}
	return (sum);
	}
inline double total_momentum (){
	int i,j;
	double  sum = 0;
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		for ( j = 0; j< N ; j++){
		sum += particleList[i].speed[j];
		}
	}
	return (sum);
	}

void fix_boundaries ( particle_s * particleList){
	int i = 0;
	int j = 0;
	for (i = 0 ; i< NUMBER_OF_PARTICLES ; i++){
		for( j= 0; j< N ; j++){
			particleList[i].position[j]= fmod(particleList[i].position[j],L);
			if ( particleList[i].position[j] < 0 ){
				particleList[i].position[j] += L;
			}
		}
	}
}

void particle_init (){
	int i_part= 0;
	double speed_cm[3];
	int i,j;
	int x,y,z;
	for ( i = 0; i<N;i++) {
		speed_cm[i]=0.0;
	}
	int cube_max=0;
	do{
	 cube_max++;
	}while( cube_max*cube_max*cube_max<NUMBER_OF_PARTICLES);
	printf("CUBE MAX = %d\n",cube_max);
	DIST_RET = L/(double)cube_max;
	printf("Dist_ret = %e\n",DIST_RET);
	//printf("%d Printed: (%e,%e,%e)\n", i_part,r_0[0],r_0[1],r_0[2]);
	for ( x=0; x<cube_max ; x++){
		for ( y = 0;  y<cube_max; y++){
			for ( z=0; z<cube_max; z++){
				if(i_part<NUMBER_OF_PARTICLES){
					particleList[i_part].position[0] = DIST_RET*(x+sqrt(2)*0.5);
					particleList[i_part].position[1] = DIST_RET*(y+sqrt(2)*0.5);
					particleList[i_part].position[2] = DIST_RET*(z+sqrt(2)*0.5);
				}
				i_part++;
			}
		}
	}
	for ( i = 0; i< NUMBER_OF_PARTICLES; i++){
		for ( j = 0; j<N;j++){
			particleList[i].speed[j] = (2*(rand()/(RAND_MAX*1.0)) - 1.0 )*V_MAX;
			speed_cm[j] += particleList[i].speed[j];
			particleList[i].acc[j]=0;
			}
		
	}
	for (i =0 ; i< NUMBER_OF_PARTICLES; i++){
		for ( j = 0; j<N;j++){
			particleList[i].speed[j] -= speed_cm[j]/((double) NUMBER_OF_PARTICLES);
		}
	}
}
/*Lennard-Jones 6-12 -> U* = U/epsilon */
inline double potential ( double r ){
	if (r<R_LIM){
		return (( 4*(1/(pow(r,12))-1/(pow(r,6)))) - u_R_LIM - DeltaF*(r-R_LIM));
	}
	else{
		return 0;
	}
}

/*Si intende forza* (adimensionata)*/
inline double force (double r){
		return (24*(2*pow(1/r,13) - pow(1/r,7))+ DeltaF);
}
// Senza liste

inline void calc_acc (){
	double r;
	double F;
	int i,j,l;
	int x,y,z;
	int found;
	particle_s temp_part;
	double r_versore[N];
	for (i=0;i<NUMBER_OF_PARTICLES;i++){
		for ( l=0;l<N;l++){
			particleList[i].acc[l] = 0;
		}
	}
	for (i=0;i<NUMBER_OF_PARTICLES;i++){
		for ( j = i+1; j < NUMBER_OF_PARTICLES;j++){
			found = 0;
			for ( x= -1; x < 2 ; x++){
				for ( y = -1; y<2 ; y++){
					for( z=-1 ; z<2 ; z++){
						if ( found == 0){
							temp_part = particleList[j];
							temp_part.position[0] += x*L;
							temp_part.position[1] += y*L;
							temp_part.position[2] += z*L;
							diff(particleList[i].position,temp_part.position ,r_versore );
							r = scalar_prod(r_versore,r_versore);
							if ( r < R_LIM*R_LIM){
								found++;
								F= force(r);
								scalar_mult(1/r,r_versore);
								for ( l=0;l<N;l++){
									particleList[i].acc[l] += F*r_versore[l];
									particleList[j].acc[l] -= F*r_versore[l];
								}
							}
						}
					}
				}
			}
		}
	}
}


inline void verlet( particle_s * partList){
	int i,j;
	for ( i = 0 ; i < NUMBER_OF_PARTICLES ; i++){
		for (j = 0; j<N;j++){
			partList[i].position[j] += partList[i].speed[j]*D_T + partList[i].acc[j]*0.5*D_T*D_T;
			partList[i].speed[j] += partList[i].acc[j]*0.5*D_T;
			
		}
	}
	fix_boundaries(particleList);
	calc_acc();
	for ( i = 0 ; i < NUMBER_OF_PARTICLES ; i++){
		for (j = 0; j<N;j++){
			partList[i].speed[j] += partList[i].acc[j]*0.5*D_T;
		}
	}
}



inline double total_energy(){
	int i,j;
	int x,y,z;
	double k_en=0;
	double u=0;
	double v[N];
	particle_s temp_part;
	for ( i = 0;i<NUMBER_OF_PARTICLES;i++){
		k_en += 0.5*scalar_prod(particleList[i].speed,particleList[i].speed);
	}
	for ( i = 0; i<NUMBER_OF_PARTICLES;i++){
		for ( j=i+1; j<NUMBER_OF_PARTICLES;j++){
			for ( x= -1; x < 2 ; x++){
				for ( y = -1; y < 2 ; y++){
					for( z=-1 ; z < 2 ; z++){
						temp_part = particleList[j];
						temp_part.position[0] += x*L;
						temp_part.position[1] += y*L;
						temp_part.position[2] += z*L;
						diff(particleList[i].position,temp_part.position,v);
						u+= potential( sqrt(scalar_prod(v,v)));
					}
				}
			}
		}
	}
	return (( k_en+u)/(double)NUMBER_OF_PARTICLES);
}

inline void riscala_vel_temp (){
	int i,j;
	double k_en = kin_en();
	for ( i = 0; i<NUMBER_OF_PARTICLES;i++){
		for (j = 0; j<N;j++){
			particleList[i].speed[j] *= sqrt( T_D/k_en);
		}
	}
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
	max=3;
	min=0;
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
	f=fopen("data/boltzmann-histo.dat","w");
	for(i=0;i<n;i++){
        fprintf(f,"%lf\t%d\t\n",min+(i+0.5)*width,freq[i]);
    }
    fclose(f);
    free(freq);
} //Giusti dixit "io voglio apparire come un nulla tenente hahahah" 


void boltzmann_file_save ( void ){
	int i = 0;
	double speed_squared = 0;
	FILE *f = fopen ("data/boltzmann.dat","a");
	for (i = 0; i< NUMBER_OF_PARTICLES ; i++){
		speed_squared = sqrt(scalar_prod(particleList[i].speed,particleList[i].speed));
		fprintf(f,"%e\n",speed_squared);
	}
	fclose(f);
}


void create_box_file(){
	FILE *f = fopen("data/box.tcl","w");
	fprintf(f,"set minx 0\n");
	fprintf(f,"set miny 0\n");
	fprintf(f,"set minz 0\n");
	fprintf(f,"set maxx %e\n",L);
	fprintf(f,"set maxy %e\n",L);
	fprintf(f,"set maxz %e\n",L);
	fprintf(f,"draw materials off \n draw color yellow \n ");
	fprintf(f,"draw line \" $minx $miny $minz \" \" $maxx $miny $minz\"\n ");
	fprintf(f,"draw line \"$minx $miny $minz\" \"$minx $maxy $minz\"\n ");
	fprintf(f,"draw line \"$minx $miny $minz\" \"$minx $miny $maxz\"\n ");
	fprintf(f,"draw line \"$maxx $miny $minz\" \"$maxx $maxy $minz\"\n ");
	fprintf(f,"draw line \"$maxx $miny $minz\" \"$maxx $miny $maxz\"\n ");
	fprintf(f,"draw line \"$minx $maxy $minz\" \"$maxx $maxy $minz\"\n" );
	fprintf(f,"draw line \"$minx $maxy $minz\" \"$minx $maxy $maxz\"\n ");
	fprintf(f,"draw line \"$minx $miny $maxz\" \"$maxx $miny $maxz\"\n ");
	fprintf(f,"draw line \"$minx $miny $maxz\" \"$minx $maxy $maxz\"\n ");
	fprintf(f,"draw line \"$maxx $maxy $maxz\" \"$maxx $maxy $minz\"\n ");
	fprintf(f,"draw line \"$maxx $maxy $maxz\" \"$minx $maxy $maxz\"\n ");
	fprintf(f,"draw line \"$maxx $maxy $maxz\" \"$maxx $miny $maxz\" \n ");
	fprintf(f,"mol default style VDW\n");
	fprintf(f,"mol addfile vmd.xyz\n");
	fprintf(f,"set sel [atomselect top all]\n$sel set radius 0.1\n");
	fclose(f);
}


int main (int argc, char *argv[]){
double rho=0.7;
int i ;
double x;
R_LIM = 2.5*SIGMA;
DeltaF = -0.039;
int iteration = 0 ;
u_R_LIM = 4*(1/(pow(R_LIM,12))-1/(pow(R_LIM,6)));
if (argc == 2){
	rho = atof(argv[1]);
}
srand(112);
L = cbrt(NUMBER_OF_PARTICLES/rho);
printf("L = %e\n",L);
printf("Frazione di impacchettamento: %e\n", rho);
fflush(stdout);
if ( R_LIM > L/2){
	printf("R_LIM > L/2\n");
	exit(1);
}
FILE *f = fopen ("data/boltzmann.dat","w+");
fclose(f);
particleList = malloc( sizeof(particle_s)*NUMBER_OF_PARTICLES);
create_box_file();
FILE *f_vmd=fopen("data/vmd.xyz","w");
particle_init();
fix_boundaries(particleList);
print_coordinate();
fprintf(f_vmd,"%d\n\n",NUMBER_OF_PARTICLES);
for( i = 0; i<NUMBER_OF_PARTICLES;i++){
	fprintf(f_vmd,"He\t%14.10e\t%14.10e\t%14.10e\n",particleList[i].position[0],particleList[i].position[1],particleList[i].position[2]);
}
riscala_vel_temp();
FILE *f_energy = fopen("data/energy.dat","w");
FILE *f_mom = fopen("data/momentum.dat","w");

while ( iteration < ITERATION_THERM){
	if (iteration %500== 0){
		printf("Iterazione %d\n",iteration);
		printf("E_TOT = %e\t P = %e\n", total_energy() ,total_momentum());
	}
	if ( iteration %10== 0){
		riscala_vel_temp();
	}
	verlet(particleList);
	total_time+=D_T;
	iteration++;
}

iteration = 0;
total_time=0;



while ( iteration < ITERATION_MAX){
	verlet(particleList);
	total_time+=D_T;
//	if (iteration % 30== 0){
		fprintf(f_vmd,"%d\n\n",NUMBER_OF_PARTICLES);
		for( i = 0; i<NUMBER_OF_PARTICLES;i++){
			fprintf(f_vmd,"He\t%14.10e\t%14.10e\t%14.10e\n",particleList[i].position[0],particleList[i].position[1],particleList[i].position[2]);
//		}
	}
//	if (iteration % 100 == 0){
		fprintf(f_energy,"%e\t%e\n",total_time,total_energy());
//		fprintf(f_mom,"%e\t%e\n",total_time,total_momentum());
//	}
//	riscala_vel_temp();
	if (iteration %1000== 0){
		printf("Iterazione %d\n",iteration);
		printf("E_TOT = %e\t P = %e\n", total_energy() ,total_momentum());
		boltzmann_file_save();
	}
	iteration++;
}
histo("data/boltzmann.dat",50);
/************* FREEE DELLA MEMORIA*/
free(particleList);
fclose(f_vmd);
return (EXIT_SUCCESS);
}
