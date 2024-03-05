#ifndef TSP_H
#define TSP_H

#include "utils.h"

typedef struct {    //node expressed as coordinates
    double x, y;
} tsp_pair;

typedef struct {    //isntance

    int nnodes; //number of nodes
    tsp_pair* coords;   //list of pair of nodes

} tsp_instance;

#endif