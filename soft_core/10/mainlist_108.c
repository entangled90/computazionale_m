#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>

#define NUMBER_OF_PARTICLES 108
#define N 3
#define ITERATION_MAX 2e4
#define ITERATION_THERM 5000

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
double R_LIST;
int last_index;
double press=0;
int iteration = 0 ;
//unsigned int NUM_TEMPI_SALVATI ;
//unsigned int PERIOD_R2 = 100;
//unsigned int time_counted;
typedef struct particle_s {
	double position[N];
	double speed[N];
	double acc[N];
	int list_start;
	} particle_s ;

particle_s * particleList;
particle_s * * neighboursList;
//particle_s * time_list;

void print_coordinate (){
	FILE *f = fopen ( "data/pack.dat","w");
	int i = 0;
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		fprintf(f,"%e\t%e\t%e\n", particleList[i].position[0], particleList[i].position[1],particleList[i].position[2]);
	}
	fclose(f);
	}
	
	void print_speed (){
	FILE *f = fopen ( "data/speed.dat","a");
	int i = 0;
	for ( i = 0; i< NUMBER_OF_PARTICLES ; i++){
		fprintf(f,"%e \t %e \t%e\n", particleList[i].speed[0], particleList[i].speed[1], particleList[i].speed[2]);
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
	return (sum/(double)NUMBER_OF_PARTICLES);
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

/*Lennard-Jones 6-12 -> U* = U/epsilon */
inline double potential ( double r ){
	if (r<R_LIM){
		return (( 4*(1/(pow(r,12))-1/(pow(r,6)))) - u_R_LIM - DeltaF*(r-R_LIM));
	}
	else{
		return 0;
	}
}
	
/* Crea liste*/
void create_list (){
	int i,j;
	int x,y,z;
	int found;
	double r_versore[N];
	particle_s temp_part;
	int list_index = 0;
	int has_neighbour=0;
	for ( i = 0; i<NUMBER_OF_PARTICLES-1;i++){
		has_neighbour=0;
		for ( j = i+1; j<NUMBER_OF_PARTICLES;j++){
			found = 0;
				for ( x= -1; x < 2 ; x++){
					for ( y = -1; y<2 ; y++){
						for( z=-1 ; z<2 ; z++){
							if (found == 0){
								temp_part = particleList[j];
								temp_part.position[0] += x*L;
								temp_part.position[1] += y*L;
								temp_part.position[2] += z*L;
								diff(temp_part.position , particleList[i].position,r_versore );
								r = scalar_prod(r_versore,r_versore);
								if(r <  R_LIST*R_LIST ){
									found++;
									neighboursList[list_index] = particleList+j;
									if(has_neighbour ==0){
										particleList[i].list_start = list_index;
									}
									list_index++;
									has_neighbour++;
								}
							}
						}
					}
				}
		}
		if ( has_neighbour == 0){
			particleList[i].list_start =list_index;
		}
	}
	last_index = list_index;
	particleList[NUMBER_OF_PARTICLES].list_start=last_index;
}
/*Si intende forza* (adimensionata)*/
inline double force (double r){
		return (24*(2*pow(1/r,13) - pow(1/r,7))+ DeltaF);
}

inline void calc_acc (){
	double r;
	double F;
	int i,j,l;
	double r_versore[N];
	int x,y,z;
	int found;
	particle_s temp_part;
	for (i=0;i<NUMBER_OF_PARTICLES;i++){
		for ( l=0;l<N;l++){
			particleList[i].acc[l] = 0;
		}
	}
	press=0;
	for (i=0;i<NUMBER_OF_PARTICLES-1;i++){
		for ( j = particleList[i].list_start; ((j < particleList[i+1].list_start) && ( j<last_index)) ;j++){
			found =0;
			for ( x= -1; x < 2 ; x++){
				for ( y = -1; y<2 ; y++){
					for( z=-1 ; z<2 ; z++){
						if (found ==0){
							temp_part = *(neighboursList[j]) ;
							temp_part.position[0] += x*L;
							temp_part.position[1] += y*L;
							temp_part.position[2] += z*L;
							diff(particleList[i].position,temp_part.position ,r_versore );
							r = scalar_prod(r_versore,r_versore);
							if ( r < R_LIM*R_LIM){
								r = sqrt(r);
								found++;
								F= force(r);
								scalar_mult(1/r,r_versore);
								for ( l=0;l<N;l++){
									particleList[i].acc[l] += F*r_versore[l];
									neighboursList[j]->acc[l] -= F*r_versore[l];
								}
								if ( iteration > ITERATION_THERM){
									press -= F*r;
									
								}
							}
						}
					}
				}
			}
		}
	}
	if( iteration >ITERATION_THERM){
		p /= (3*NUMBER_OF_PARTICLES*kin_en());
		pression_save();
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

inline double potential_energy(){
	int i,j;
	double u=0;
	double v[N];
	int found = 0;
	int x,y,z;
	particle_s temp_part;
	for ( i = 0; i<NUMBER_OF_PARTICLES-1;i++){
		for ( j = particleList[i].list_start; ((j < particleList[i+1].list_start) && ( j<last_index)) ;j++){
			found =0;
			for ( x= -1; x < 2 ; x++){
				for ( y = -1; y<2 ; y++){
					for( z=-1 ; z<2 ; z++){
						if (found ==0){
							temp_part = *(neighboursList[j]) ;
							temp_part.position[0] += x*L;
							temp_part.position[1] += y*L;
							temp_part.position[2] += z*L;
							diff(particleList[i].position,temp_part.position ,v );
							r = scalar_prod(v,v);
							if ( r < R_LIM*R_LIM){
								r = sqrt(r);
								u+= potential( sqrt(scalar_prod(v,v)));
								found++;
							}
						}
					}
				}
			}
		}
	}
	u /= (double) NUMBER_OF_PARTICLES;
	return u;
}


inline double total_energy(){
	double k_en,u;
	k_en = kin_en();
	u=potential_energy();
	return ( k_en+u);
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






inline void vmd_file_save(){
	int i ;
	FILE *f_vmd = fopen("data/vmd.xyz","a");
	fprintf(f_vmd,"%d\n\n",NUMBER_OF_PARTICLES);
	for( i = 0; i<NUMBER_OF_PARTICLES;i++){
			fprintf(f_vmd,"He\t%14.10e\t%14.10e\t%14.10e\n",particleList[i].position[0],particleList[i].position[1],particleList[i].position[2]);
	}
	fclose(f_vmd);
}

inline void energy_file_save(){
	FILE *f_energy = fopen("data/energy_108.dat","a");
	fprintf(f_energy,"%e\t%e\n",total_time,total_energy());
	fclose(f_energy);
}

inline void momentum_file_save(){
	FILE *f_mom = fopen("data/momentum.dat","a");
	fprintf(f_mom,"%e\t%e\n",total_time,total_momentum());
	fclose(f_mom);
}

inline void pression_save(){
	FILE *f_pres = fopen("data/press.dat","a");
	fprintf(f_pres,"%e\t%e\n",total_time,press);
	fclose(f_pres);
}

int main (int argc, char *argv[]){
double rho=0.7;
R_LIM = 2.5*SIGMA;
R_LIST = 2.8*SIGMA;
DeltaF = -0.039;
int i;
u_R_LIM = 4*(1/(pow(R_LIM,12))-1/(pow(R_LIM,6)));
if (argc == 2){
	rho = atof(argv[1]);
}
srand(time(NULL));
L = cbrt(NUMBER_OF_PARTICLES/rho);
if (R_LIM > L/2){
	printf("R_LIM > L mezzi\n");
	exit(1);
}
//NUM_TEMPI_SALVATI = (ITERATION_MAX)/( (double) PERIOD_R2) +1;
printf("L = %e\n",L);
printf("Frazione di impacchettamento: %e\n DIST_RET = %e\n", rho,DIST_RET);
fflush(stdout);
FILE *f_pres = fopen("data/press.dat","w");
fclose(f_pres);
particleList = malloc( sizeof(particle_s)*NUMBER_OF_PARTICLES);
neighboursList = malloc(sizeof(particle_s)*30*NUMBER_OF_PARTICLES);
//time_list = malloc(sizeof(particle_s)*NUM_TEMPI_SALVATI*NUMBER_OF_PARTICLES);
particle_init();
fix_boundaries(particleList);
//print_coordinate();
create_list();
riscala_vel_temp();
create_box_file();
FILE *f_energy = fopen("data/energy_500.dat","w");
fclose(f_energy);
/*
 * FILE *f_mom = fopen("data/momentum.dat","w+");
FILE *f_vmd=fopen("data/vmd.xyz","w+");
fclose(f_vmd);
fclose(f_mom);
*/
while ( iteration < ITERATION_THERM){
	if (iteration %500== 0){
		printf("Iterazione %d\n",iteration);
		printf("E_TOT = %e\t P = %e\n", total_energy() ,total_momentum());
	}
	riscala_vel_temp();
	if ( iteration %10== 0){
		create_list();
	}
	verlet(particleList);
	energy_file_save();
	total_time+=D_T;
	iteration++;
}

//iteration = 0;
//total_time=0;

while ( iteration < ITERATION_MAX+ITERATION_THERM){
	if (iteration %1000== 0){
		printf("Iterazione %d\n",iteration);
		printf("E_TOT = %e\t P = %e\n", total_energy() ,total_momentum());
//		boltzmann_file_save();
	}

	if ( iteration %10== 0){
		create_list();
	}
//	riscala_vel_temp();
	verlet(particleList);
	total_time+=D_T;
//	vmd_file_save();
//	momentum_file_save();
	energy_file_save();
	iteration++;
}
printf("Calcolo r2\n");
//r_squared_save("data/r2.dat");

free(neighboursList);
free(particleList);


return (EXIT_SUCCESS);
}
