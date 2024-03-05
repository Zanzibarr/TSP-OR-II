#ifndef TSP_H
#define TSP_H

#include "utils.h"


typedef struct {    //node expressed as coordinates
    double x, y;
} tsp_pair;

typedef struct {    //instance

    int nnodes;             //number of nodes
    tsp_pair* coords;       //list of nodes
    double** costs;         //cost matrix

    int* tsp_best_solution; //store the best solution found so far
    double tsp_best_cost;   //store the best cost found so far
    int tsp_best_time;      //store the time of the best solution found so far

} tsp_instance;

// SOLVING
char tsp_alg_type[10];  //store the type of algorithm using



void tsp_allocate_coords_space(tsp_instance* inst) {    //dinamically allocate the space for the coords list

    if (!inst -> nnodes) { printf("The nnodes variable hasn't been assigned yet."); exit(1); }

    inst -> coords = calloc(inst -> nnodes, sizeof(tsp_pair));

}

void tsp_allocate_costs_space(tsp_instance* inst) { //dinamicallu allocate the space for the costs matrix

    if (!inst -> nnodes) { printf("The nnodes variable hasn't been assigned yet."); exit(1); }

    inst -> costs = calloc(inst -> nnodes, sizeof(double*));    //dinamically allocate the rows

    for (int i = 0; i < inst -> nnodes; i++)
        inst -> costs[i] = calloc(inst -> nnodes, sizeof(double)); //dinamically allocate the columns

}

void tsp_allocate_best_sol_space(tsp_instance* inst) {  //dinamically allocate the space for the best solution list

    if (!inst -> nnodes) { printf("The nnodes variable hasn't been assigned yet."); exit(1); }

    inst -> tsp_best_solution = calloc(inst -> nnodes, sizeof(int));
    
}

double tsp_compute_distance(const tsp_instance* inst, int i, int j) {   //euclidian distance between two points in the instance
    
    return sqrt(pow(inst -> coords[i].x - inst -> coords[j].x, 2) + pow(inst -> coords[i].y - inst -> coords[j].y, 2));
    
}

void tsp_precompute_costs(tsp_instance* inst) { //precompute the costs of the edges

    tsp_allocate_costs_space(inst);  //dimanically allocate the space for the cost matrix;  

    for (int i = 0; i < inst -> nnodes; i++) for (int j = 0; j < inst -> nnodes; j++)
        inst -> costs[i][j] = tsp_compute_distance(inst, i, j); //compute the cost of each edge

}

void tsp_init_solution(tsp_instance* inst) {    //initialize the best solution

    tsp_allocate_best_sol_space(inst);
    inst -> tsp_best_cost = INFINITY;
    inst -> tsp_best_time = 0;

}

void tsp_free_instance(tsp_instance* inst) {    //frees the dinamically allocated memory

    free (inst -> coords);
    free (inst -> costs);
    free (inst -> tsp_best_solution);

}

#endif