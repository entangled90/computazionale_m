#ifndef SW_H
#define SW_H

//#ifndef SW_C
#include "list.h"
void spin_init ( Spin * matrix, Node * n);
void reset_cluster (Spin * matrix, Node * n);
int set_bond (Spin * s1, Spin * s2);
void fillCluster( Spin * matrix, Node * nodes, List * l);
void startClustering (Spin * matrix, Node * nodes);
inline double magnetization(Spin *x);
inline void savePPM(Spin * s);
void drawCluster(Spin * s);
void drawSpin(Spin * s);	
void flip_spin ( Spin * m);
void print_data (Spin * m);
void evolve_therm (Spin * matrix, Node * nodes);
void evolve( Spin * matrix, Node * nodes);

//#endif
#endif
