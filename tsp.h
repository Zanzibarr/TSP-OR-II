#ifndef TSP_H
#define TSP_H

#include "utils.h"

typedef struct {    //node expressed as coordinates
    double x, y;
} tsp_pair;

typedef struct {    //isntance

    int nnodes; //number of nodes
    tsp_pair* coords;   //list of pair of nodes
    double** costs; //cost matrix

} tsp_instance;

// SOLVING
char tsp_alg_type[10];  //store the type of algorithm using
int* tsp_best_solution; //store the best solution found so far
double tsp_best_cost;   //store the best cost found so far
int tsp_best_time;  //store the time of the best solution found so far

double tsp_compute_distance(const tsp_instance* inst, int i, int j) { return sqrt(pow(inst -> coords[i].x - inst -> coords[j].x, 2) + pow(inst -> coords[i].y - inst -> coords[j].y, 2)); }

#endif