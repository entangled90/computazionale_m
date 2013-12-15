#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define N 50000
#define omega 1e6
int main (){
  double time[N];
  int i=0;
  for(i = 0; i<N;i++){
    time[i] = (rand()/(double)RAND_MAX)*20-10;
  }
  FILE *f =fopen("tempi.dat","w");
  for ( i = 0; i<N;i++){
    fprintf(f,"%e\n",sin(omega*time[i]));
  }
  fclose(f);
  exit(0);
}
