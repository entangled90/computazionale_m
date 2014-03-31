#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>



int main (){
	FILE * f=fopen("data/prova.dat","w");
	fprintf(f,"prova\n");
	fclose(f);


	return 0;
}