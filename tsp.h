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
char tsp_alg_type[10];
int* tsp_best_solution;
double tsp_best_cost;
int tsp_best_time;

double tsp_compute_distance(const tsp_instance* inst, int i, int j) { return sqrt(pow(inst -> coords[i].x - inst -> coords[j].x, 2) + pow(inst -> coords[i].y - inst -> coords[j].y, 2)); }

#endif