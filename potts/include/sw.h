#ifndef SW_H
#define SW_H

//#ifndef SW_C
#include "list.h"
#include "cnum.h"

void spin_init ( Spin * matrix, Node * n, int N);
void reset_cluster (Spin * matrix, Node * n, int N);
int set_bond (Spin * s1, Spin * s2,float BETA);
void fillCluster( Spin * matrix, Node * nodes, List * l, int N, float BETA);
void startClustering (Spin * matrix, Node * nodes, int N,float BETA);
inline double magnetization(Spin *x, int N);
inline void savePPM(Spin * s, int N, char *);
void drawCluster(Spin * s, int N);
void drawSpin(Spin * s, int N);	
void flip_spin ( Spin * m, int N);
void print_data (Spin * m, int N);
void evolve_therm (Spin * matrix, Node * nodes,int , float);
void evolve( Spin * matrix, Node * nodes, int, float);
double hamiltoniana(Spin * configuration, int N);
cNum sum_row(Spin * s, int row, int N);
cNum sum_col(Spin * s, int col, int N);

//#endif
#endif