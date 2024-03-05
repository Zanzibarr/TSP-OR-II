#ifndef TSP_H
#define TSP_H

#include "utils.h"



typedef struct {    //node expressed as coordinates
    double x, y;
} tsp_pair;

typedef struct {    //isntance

    int nnodes;             //number of nodes
    tsp_pair* coords;       //list of pair of nodes
    double** costs;         //cost matrix

    int* tsp_best_solution; //store the best solution found so far
    double tsp_best_cost;   //store the best cost found so far
    int tsp_best_time;      //store the time of the best solution found so far

} tsp_instance;

// SOLVING
char tsp_alg_type[10];  //store the type of algorithm using



//dinamically allocate the space for the coords list
void tsp_allocate_coords_space(tsp_instance* inst) {

    if (!inst -> nnodes) { printf("The nnodes variable hasn't been assigned yet."); exit(1); }

    inst -> coords = calloc(inst -> nnodes, sizeof(tsp_pair));

}

//dinamicallu allocate the space for the costs matrix
void tsp_allocate_costs_space(tsp_instance* inst) {

    if (!inst -> nnodes) { printf("The nnodes variable hasn't been assigned yet."); exit(1); }

    inst -> costs = calloc(inst -> nnodes, sizeof(double*));    //dinamically allocate the rows

    for (int i = 0; i < inst -> nnodes; i++)
        inst -> costs[i] = calloc(inst -> nnodes, sizeof(double)); //dinamically allocate the columns

}

//dinamically allocate the space for the best solution list
void tsp_allocate_best_sol_space(tsp_instance* inst) {

    if (!inst -> nnodes) { printf("The nnodes variable hasn't been assigned yet."); exit(1); }

    inst -> tsp_best_solution = calloc(inst -> nnodes, sizeof(int));
    
}

//euclidian distance between two points in the instance
double tsp_compute_distance(const tsp_instance* inst, int i, int j) {
    
    return sqrt(pow(inst -> coords[i].x - inst -> coords[j].x, 2) + pow(inst -> coords[i].y - inst -> coords[j].y, 2));
    
}

//precompute the costs of the edges
void tsp_precompute_costs(tsp_instance* inst) {

    tsp_allocate_costs_space(inst);  //dimanically allocate the space for the cost matrix;  

    for (int i = 0; i < inst -> nnodes; i++) for (int j = 0; j < inst -> nnodes; j++)
        inst -> costs[i][j] = tsp_compute_distance(inst, i, j); //compute the cost of each edge

}

//initialize the best solution
void tsp_init_solution(tsp_instance* inst) {

    tsp_allocate_best_sol_space(inst);
    inst -> tsp_best_cost = INFINITY;
    inst -> tsp_best_time = 0;

}

//frees the dinamically allocated memory
void tsp_free_instance(tsp_instance* inst) {

    free (inst -> coords);
    free (inst -> costs);
    free (inst -> tsp_best_solution);

}

#endif