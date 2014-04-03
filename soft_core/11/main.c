#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>

#define NUMBER_OF_PARTICLES 256
#define N 3
#define ITERATION_MAX 1e5
#define ITERATION_THERM 5000
#define skip_times 25
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
double T_D = 0.8;
double R_LIST;
int last_index;
double DeltaT;
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
particle_s * time_list;
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
/*Lennard-Jones 6-12 -> U* = U/epsilon */
inline double potential ( double r ){
	if (r<R_LIM){
		double tmp=pow(r,6);
		return (( 4*(1/(tmp*tmp)-1/(tmp))) - u_R_LIM - DeltaF*(r-R_LIM));
	}
	else{
		return 0;
	}
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
/*******************************************************************************************/

void genera_sottoreticolo(double rx_in, double ry_in,double rz_in,int q,int start, double passo){
	int p = start; 
	int c=0;
	int d=0;
	int e=0;
	int j;
	double rx;
	double ry;
	double rz;
	rx=rx_in;
	ry=ry_in;
	rz=rz_in;
	while(p<NUMBER_OF_PARTICLES && e<q-1){
		while (p < NUMBER_OF_PARTICLES && d<q-1) { 
			while (p < NUMBER_OF_PARTICLES && c<q-1) {
				particleList[p].position[0]=rx;
				particleList[p].position[1]=ry;		
				particleList[p].position[2]=rz;
			for ( j = 0; j<N;j++){
				particleList[p].speed[j] = (2*(rand()/(RAND_MAX*1.0)) - 1.0 )*V_MAX;
				particleList[p].acc[j]=0;
				}
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
    	m = q*q*q;
    }
    passo = cbrt(pow(L,3)/(double)(m));
	printf("passo %lf\n", passo);
 	//creazione reticolo
	//printf("Primo reticolo\n");
  	genera_sottoreticolo(0,0,0,q,0,passo);
//	printf("Secondo reticolo\n");
//	genera_sottoreticolo(passo/2.0,passo/2.0,passo/2.0,q,NUMBER_OF_PARTICLES/2, passo);
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
	double temp = 2/3.0*kin_en();
	for ( i = 0; i<NUMBER_OF_PARTICLES;i++){
		for (j = 0; j<N;j++){
			particleList[i].speed[j] *= sqrt(T_D/temp);
			//printf("Correction %lf\n",sqrt(T_D/temp));
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


void print_vec(char * file, double * vec, int Len){
	int i=0;
	FILE *f=fopen(file,"w");
	for(i=0;i<Len;i++){
		fprintf(f,"%e\t%e\n",i*D_T,vec[i]);
	}
	fclose(f);
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
	for ( delta = 1; delta  <  ITERATION_MAX/skip_times; delta++){
		sum = 0;
		count = 0;
		for ( init = 0; init+delta<ITERATION_MAX/skip_times; init++){
			tmp=r_squared_calc( time_list+(init+delta)*NUMBER_OF_PARTICLES,time_list + init*NUMBER_OF_PARTICLES);
			sum += tmp;
			var +=tmp*tmp;
			count++;
		}
		sum /= (double) count;
		var /= (double) count;
		var -= sum*sum;
		fprintf(f,"%.14e\t%.14e\t%.14e\n",delta*DeltaT, sum/(L*L), sqrt(var/(double)count)/L);
	}
	fclose(f);
}



inline void copyList ( particle_s * in , particle_s * out){
	unsigned int i;
	for ( i = 0; i< NUMBER_OF_PARTICLES;i++){
		out[i] = in[i];
	}
}




int main (int argc, char *argv[]){
double rho=0.9;
DeltaT=skip_times*D_T;
R_LIM = 2.5*SIGMA;
R_LIST = 2.8*SIGMA;
DeltaF = 0.0389994774528;
int iteration = 0 ;
u_R_LIM = 4*(1/(pow(R_LIM,12))-1/(pow(R_LIM,6)));
srand(time(NULL));
L = cbrt(NUMBER_OF_PARTICLES/rho);
if(argc==2){
	T_D=atof(argv[1]);
}
if (R_LIM > L/2){
	printf("R_LIM > L mezzi\n");
	exit(1);
}
//NUM_TEMPI_SALVATI = (ITERATION_MAX)/( (double) PERIOD_R2) +1;
printf("L = %e\n",L);
printf("Frazione di impacchettamento: %e\n\n", rho);
particleList = malloc( sizeof(particle_s)*NUMBER_OF_PARTICLES);
time_list=malloc(sizeof(particle_s)*NUMBER_OF_PARTICLES*ITERATION_MAX/skip_times);
neighboursList = malloc(sizeof(particle_s)*30*NUMBER_OF_PARTICLES);
//time_list = malloc(sizeof(particle_s)*NUM_TEMPI_SALVATI*NUMBER_OF_PARTICLES);
reticolo();
fix_boundaries(particleList);
print_coordinate();
create_list();
riscala_vel_temp();
create_box_file();
unsigned int i;


/**********************FILES ***************************/
char r2_filename[64] = "";
snprintf(r2_filename,64,"data/dr2/dr2_%d_%.5lf.dat",(int)time(NULL),T_D); 

/*
 FILE *f_mom = fopen("data/momentum.dat","w+");
fclose(f_vmd);
fclose(f_mom);
*/
//FILE *f_vmd=fopen("data/vmd.xyz","w+");

total_time=0;
//char  energy_therm_filename[128] = "data/energy_therm.dat";
//FILE * f_energy_therm = fopen(energy_therm_filename,"a");
//double * energy_vec = malloc(sizeof(double)*ITERATION_THERM);
for(iteration=0;iteration<ITERATION_THERM;iteration++){
	if (iteration %2000== 0){
		printf("Iterazione %d\n",iteration);
		printf("TEMP = %e \t E_TOT = %e\t P = %e\n",2/3.0*kin_en(), total_energy() ,total_momentum());
	}
	riscala_vel_temp();
	if ( iteration %10== 0){
		create_list();
	}
	//energy_vec[iteration]=total_energy()/EPS;
	verlet(particleList);
	total_time+=D_T;
}
//print_vec(energy_therm_filename,energy_vec,ITERATION_THERM);
//free(energy_vec);
//fclose(f_energy_therm);
iteration = 0;
total_time=0;
char  energy_filename[128] = "";
snprintf(energy_filename,128,"data/energy/energy%d.dat",NUMBER_OF_PARTICLES);
//energy_vec = malloc(sizeof(double)*ITERATION_MAX);
while ( iteration < ITERATION_MAX){
	if (iteration %4000== 0){
		printf("Iterazione %d\n",iteration);
		printf("TEMP = %e \t E_TOT = %e\t P = %e\n",2/3.0*kin_en(), total_energy() ,total_momentum());
		boltzmann_file_save();
	}
	if (iteration%skip_times==0){
		for ( i = 0; i< NUMBER_OF_PARTICLES;i++){
			time_list[iteration/skip_times*NUMBER_OF_PARTICLES+i] = particleList[i];
		}	
	}
	
	if ( iteration %10== 0){
		create_list();
	//	vmd_file_save();
	}
	verlet(particleList);
	total_time+=D_T;
//	vmd_file_save();
	//energy_vec[iteration] = total_energy()/EPS;
//	fprintf(f_energy,"%e\t%e\n",total_time,tmp);
	iteration++;
}
//print_vec(energy_filename,energy_vec,ITERATION_MAX);
printf("Calcolo r2\n");
r_squared_save(r2_filename);
//free(energy_vec);
free(neighboursList);
free(particleList);
free(time_list);


return (EXIT_SUCCESS);
}
