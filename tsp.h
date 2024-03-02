#ifndef TSP_H
#define TSP_H

#include "utils.h"

typedef struct {
    double x, y;
} pair;

typedef struct {

    //data
    int nnodes;
    pair* coords;

} instance;

#endif