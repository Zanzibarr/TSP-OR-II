#ifndef TSP_H
#define TSP_H

#include "utils.h"

typedef struct {    //instance

    int nnodes;             //number of nodes
    tsp_pair* coords;       //list of nodes

    double* costs;         //cost "matrix"
    int* sort_edges;       //min edges "matrix": row i contains a permutation of the nodes, ordered by increasing distance from node i

    int* best_solution; //store the best solution found so far
    double best_cost;   //store the best cost found so far
    double best_time;   //store the time of the best solution found so far (in seconds)

} tsp_instance;

// GLOBAL VARIABLES
int tsp_over_time;
uint64_t tsp_seed;

double tsp_total_time;

char tsp_alg_type[10];
char tsp_file_name[100];

time_t tsp_time_limit;


// USEFUL METHODS

void tsp_update_best_sol(tsp_instance* inst, int* path, double cost, double time) {

    for (int i = 0; i < inst -> nnodes; i++) inst -> best_solution[i] = path[i];
    inst -> best_cost = cost;
    inst -> best_time = time;

}

void tsp_allocate_coords_space(tsp_instance* inst) {    //dinamically allocate the space for the coords list

    if (!inst -> nnodes) { printf("The nnodes variable hasn't been assigned yet."); exit(1); }

    inst -> coords = (tsp_pair*)calloc(inst -> nnodes, sizeof(tsp_pair));

}

void tsp_allocate_costs_space(tsp_instance* inst) { //dinamicallu allocate the space for the costs matrix

    if (!inst -> nnodes) { printf("The nnodes variable hasn't been assigned yet."); exit(1); }

    inst -> costs = (double*)calloc(inst -> nnodes * inst -> nnodes, sizeof(double));

}

void tsp_allocate_best_sol_space(tsp_instance* inst) {  //dinamically allocate the space for the best solution list

    if (!inst -> nnodes) { printf("The nnodes variable hasn't been assigned yet."); exit(1); }

    inst -> best_solution = (int*)calloc(inst -> nnodes, sizeof(int));
    
}

double tsp_compute_distance(const tsp_instance* inst, int i, int j) {   //euclidian distance between two points in the instance
    
    return sqrt(pow(inst -> coords[i].x - inst -> coords[j].x, 2) + pow(inst -> coords[i].y - inst -> coords[j].y, 2));
    
}

void tsp_precompute_costs(tsp_instance* inst) { //precompute the costs of the edges

    tsp_allocate_costs_space(inst);

    for (int i = 0; i < inst -> nnodes; i++) for (int j = 0; j < inst -> nnodes; j++)
        inst -> costs[i * inst -> nnodes + j] = tsp_compute_distance(inst, i, j);

}

void tsp_merge(tsp_entry* list, int p, int q, int r) {  //sort merge part of the mergesort

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

}

void tsp_sort(tsp_entry* list, int p, int r) {  //mergesort comparing the values of the entries (cost of the edge)
    
    if (p < r) {
        int q = (p+r)/2;
        tsp_sort(list, p, q);
        tsp_sort(list, q+1, r);
        tsp_merge(list, p, q, r);
    }

}

void tsp_check_sort_edges(const tsp_instance* inst) {
    
    for (int i = 0; i < inst -> nnodes; i++) {
        double min_cost = 0;
        for (int j = 0; j < inst -> nnodes - 1; j++) {
            double checked_cost = inst -> costs[i * inst -> nnodes + inst -> sort_edges[i * (inst -> nnodes - 1) + j]];
            if (checked_cost < min_cost - TSP_EPSYLON) {
                printf("SORT_EDGES INTEGRITY COMPROMISED\n");
                exit(1);
            }
            min_cost = checked_cost;
        }
    }

    #if TSP_VERBOSE >= 200
    printf("sort_edges integrity check passed.\n");
    #endif

}

int compare( const void* arg1, const void* arg2) {
    return ((tsp_entry*)arg1)->value - ((tsp_entry*)arg2)->value;
}

void tsp_precompute_sort_edges(tsp_instance* inst) {

    inst -> sort_edges = (int*)calloc(inst -> nnodes * (inst -> nnodes - 1), sizeof(int));
    tsp_entry* list = (tsp_entry*)calloc(inst -> nnodes, sizeof(tsp_entry));

    for (int i = 0; i < inst -> nnodes; i++) {

        for (int j = 0; j < inst -> nnodes; j++) {  //saving the entries to be sorted
            list[j].key = j;
            list[j].value = inst -> costs[i * inst -> nnodes + j];    //considering only the costs of the edges leaving node i
        }

        qsort((void*)list, (size_t)inst -> nnodes, sizeof(tsp_entry), compare);
        //tsp_sort(list, 0, inst -> nnodes - 1);  //sort by cost of the edge (mergesort)

        for (int j = 1; j < inst -> nnodes; j++)
            inst -> sort_edges[i * (inst -> nnodes - 1) + j-1] = list[j].key;    //populate the ith row with the nodes ordered by increasing distance

    }

    free(list);

    #if TSP_VERBOSE >= 100
    tsp_check_sort_edges(inst);
    #endif

}

void tsp_init_solution(tsp_instance* inst) {    //initialize the best solution

    tsp_allocate_best_sol_space(inst);
    inst -> best_cost = INFINITY;
    inst -> best_time = 0;

    tsp_over_time = 0;

    tsp_initial_time = clock();

}

void tsp_free_instance(tsp_instance* inst) {    //frees the dinamically allocated memory

    free (inst -> coords);
    free (inst -> costs);
    free (inst -> sort_edges);
    free (inst -> best_solution);

}

#endif