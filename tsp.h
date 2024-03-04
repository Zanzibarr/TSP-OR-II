#ifndef TSP_H
#define TSP_H

#include "utils.h"

typedef struct {    //node expressed as coordinates
    double x, y;
} pair;

typedef struct {    //isntance

    int nnodes; //number of nodes
    pair* coords;   //list of pair of nodes

} instance;

#endif