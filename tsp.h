#ifndef TSP_H
#define TSP_H

#include "utils.h"

typedef struct {    //instance

    int nnodes;             //number of nodes
    tsp_pair* coords;       //list of nodes

    double** costs;         //cost matrix
    int** min_edges;        //min edges matrix

    int* tsp_best_solution; //store the best solution found so far
    double tsp_best_cost;   //store the best cost found so far
    double tsp_best_time;   //store the time of the best solution found so far (in seconds)

} tsp_instance;

void tsp_update_best_sol(tsp_instance* inst, int* path) {

    for (int i = 0; i < inst -> nnodes; i++) inst -> tsp_best_solution[i] = path[i];

}

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

void tsp_merge(tsp_entry* list, int p, int q, int r) {

    int i = p, j = q+1, k=0;
    tsp_entry b[r-p+1];

    while (i <= q && j <= r) {
        if (list[i].value < list[j].value) {
            b[k] = list[i];
            i++;
        } else {
            b[k] = list[j];
            j++;
        }
        k++;
    }

    while (i <= q) {
        b[k] = list[i];
        i++; k++;
    }

    while (j <= r) {
        b[k] = list[j];
        j++; k++;
    }

    for (k=p; k<=r; k++)
        list[k] = b[k-p];

    return;

}

void tsp_sort(tsp_entry* list, int p, int r) {
    
    if (p < r) {
        int q = (p+r)/2;
        tsp_sort(list, p, q);
        tsp_sort(list, q+1, r);
        tsp_merge(list, p, q, r);
    }

    return;

}

void tsp_precompute_min_edges(tsp_instance* inst) {

    inst -> min_edges = calloc(inst -> nnodes, sizeof(tsp_entry*));

    for (int i = 0; i < inst -> nnodes; i++) {

        tsp_entry list[inst -> nnodes];

        for (int j = 0; j < inst -> nnodes; j++) {
            list[j].key = j;
            list[j].value = inst -> costs[i][j];
        }

        tsp_sort(list, 0, inst -> nnodes - 1);
        
        inst -> min_edges[i] = calloc(inst -> nnodes - 1, sizeof(tsp_entry));

        for (int j = 1; j < inst -> nnodes; j++)
            inst -> min_edges[i][j-1] = list[j].key;

    }

}

void tsp_init_solution(tsp_instance* inst) {    //initialize the best solution

    tsp_allocate_best_sol_space(inst);
    inst -> tsp_best_cost = INFINITY;
    inst -> tsp_best_time = 0;

    tsp_initial_time = clock();

}

void tsp_free_instance(tsp_instance* inst) {    //frees the dinamically allocated memory

    free (inst -> coords);
    free (inst -> costs);
    free (inst -> min_edges);
    free (inst -> tsp_best_solution);

}

#endif