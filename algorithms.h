#include "tsp.h"

//ALGORITHMS
int     tsp_solve_greedy(tsp_instance* inst, const char g2opt);

void    tsp_check_integrity(const tsp_instance* inst, double cost, int* path);

double  tsp_greedy_from_node(const tsp_instance* inst, int* path, int start_node);
int     tsp_find_min_edge(const tsp_instance* inst, int frontier, int* path);
void    tsp_2opt(const tsp_instance* inst, int* path, double* cost);
double  tsp_find_2optswap(const tsp_instance* inst, int* path, double* cost);
void    tsp_reverse(int* path, int start, int end);


int tsp_solve_greedy(tsp_instance* inst, const char g2opt) {    //solve using greedy (+ 2opt algorithm if g2opt==1)

    int* path = (int*)calloc(inst -> nnodes, sizeof(int));

    for (int i = 0; i < inst -> nnodes; i++) {  //for each starting node
        
        double cost = tsp_greedy_from_node(inst, path, i);

        #if TSP_VERBOSE >= 10
        printf("Greedy solution from %4d:%15.4f", i, cost);
        #endif

        if (g2opt) {   //if I want to use the 2opt

            tsp_2opt(inst, path, &cost);  //fix solution using 2opt

            #if TSP_VERBOSE >= 10
            printf("\t-(2opt)->\t%15.4f", cost);
            #endif

            #if TSP_VERBOSE >= 100
            tsp_check_integrity(inst, cost, path);
            #endif

        }

        double time = tsp_time_elapsed();
        printf("\t\t/\tTime elapsed: %15.4f\n", time);
        
        if (cost < inst -> best_cost - TSP_EPSYLON) {   //if this solution is better than the best seen so far update it

            #if TSP_VERBOSE >= 10
            printf("New best solution\n");
            #endif

            tsp_update_best_sol(inst, path, cost, time);
        }

        if (time > tsp_time_limit) { free(path); return -1; } //if I exceeded the time limit

    }

    free(path);

    return 0;

}


void tsp_check_integrity(const tsp_instance* inst, double cost, int* path) {  //debugging tool

    double c_cost = 0;
    int first = path[0], error = 0;
    int* visited = (int*)calloc(inst -> nnodes, sizeof(int));

    for (int i = 1; i < inst -> nnodes; i++) {
        if (path[i] < 0 || path[i] >= inst -> nnodes || visited[path[i]]) { error = 1; break; }
        visited[path[i]] += 1;
        c_cost += inst -> costs[path[i-1] * inst -> nnodes + path[i]];
    }

    free(visited);

    c_cost += inst -> costs[path[inst -> nnodes - 1] * inst -> nnodes + path[0]];

    if (fabs(c_cost - cost) > TSP_EPSYLON) error = 1;

    if (error == 1) {
        printf("\nINTEGRITY COMPROMISED\n");
        printf("\nCost: %15.4f, Checked cost: %15.4f\n", cost, c_cost);
        exit(1);
    } else {
        #if TSP_VERBOSE >= 200
        printf("Integrity check passed.\n");
        #endif
    }

}


double tsp_greedy_from_node(const tsp_instance* inst, int* path, int start_node) {    //greedy solution starting from specific node

    int frontier = start_node, next = -1;
    double cost = 0;
    int* visited = (int*)calloc(inst -> nnodes, sizeof(int));

    path[0] = start_node;
    visited[start_node] = 1;

    for (int i = 1; i < inst -> nnodes; i++) {  //repeat nnodes times

        for (int j = 0; j < inst -> nnodes - 1; j++) {

            next = inst -> sort_edges[frontier * (inst -> nnodes - 1) + j]; //check the min_edges in the tsp_instance struct for more info
            if (!visited[next]) break;    //if I didn't explore that node yet then it's the closest (new) node

        }
        
        path[i] = next;
        visited[next] = 1;
        cost += inst -> costs[frontier * inst -> nnodes + next];
        frontier = next;

    }

    free(visited);

    cost += inst -> costs[frontier * inst -> nnodes + start_node];    //add the cost of the last edge

    return cost;

}

void tsp_2opt(const tsp_instance* inst, int* path, double* cost) {   //2opt algorithm

    while (tsp_find_2optswap(inst, path, cost) > 0);   //repeat until I can't find new swaps that improve the cost of my solution

}

double tsp_find_2optswap(const tsp_instance* inst, int* path, double* cost) {  //locating a swap for 2opt algorithm

    for (int i = 0; i < inst -> nnodes - 1; i++) {
        for (int j = i + 1; j < inst -> nnodes; j++) {
            int k = (j+1 == inst -> nnodes) ? 0 : j+1;  //allow for the loop over the edge

            if ((inst -> costs[path[i] * inst -> nnodes + path[i+1]] + inst -> costs[path[j] * inst -> nnodes + path[k]]) > (inst -> costs[path[i] * inst -> nnodes + path[j]] + inst -> costs[path[i+1] * inst -> nnodes + path[k]]) + TSP_EPSYLON){   //if c_i,i+1 + c_j,j+1 > c_i,j + c_i+1,j+1 then swap

                *cost = *cost - (inst -> costs[path[i] * inst -> nnodes + path[i+1]] + inst -> costs[path[j] * inst -> nnodes + path[k]]) + (inst -> costs[path[i] * inst -> nnodes + path[j]] + inst -> costs[path[i+1] * inst -> nnodes + path[k]]);    //update the cost

                tsp_reverse(path, i+1, j);  //reverse the part in the middle of the swap

                return 1;    //found swap

            }
        }
    }

    return -1;  //no swap found

}

void tsp_reverse(int* path, int start, int end) {   //reverse the array specified between two specified indexes

    int c;

    while (start < end) {
        c = path[start];
        path[start] = path[end];
        path[end] = c;
        start++; end--;
    }

}
