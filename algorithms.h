#include "tsp.h"

//ALGORITHMS
int     tsp_solve_greedy(tsp_instance* inst, const char g2opt);

void    tsp_check_integrity(const tsp_instance* inst, double cost, int* path);

double  tsp_greedy_from_node(const tsp_instance* inst, int* path, int start_node);
int     tsp_find_min_edge(const tsp_instance* inst, int frontier, int* path);
double  tsp_2opt(const tsp_instance* inst, int* path, double cost);
double  tsp_find_2optswap(const tsp_instance* inst, int* path, double cost);
void    tsp_reverse(int* path, int start, int end);


int tsp_solve_greedy(tsp_instance* inst, const char g2opt) {    //solve using greedy (+ 2opt algorithm if g2opt==1)

    double time, cost;
    int* path = calloc(inst -> nnodes, sizeof(int));

    for (int i = 0; i < inst -> nnodes; i++) {  //for each starting node
        
        cost = tsp_greedy_from_node(inst, path, i);
        if (TSP_VERBOSE > 0) tsp_check_integrity(inst, cost, path);

        if (TSP_VERBOSE > 0) {
            printf("\nIntermediate solution from %d: ", i);
            for (int j = 0; j < inst -> nnodes; j++) printf("%d ", path[j]);
            printf("\nCost: %f\n", cost);
        }

        if (g2opt == 1) {   //if I want to use the 2opt

            cost = tsp_2opt(inst, path, cost);  //fix solution using 2opt
            if (TSP_VERBOSE > 0) tsp_check_integrity(inst, cost, path);

            if (TSP_VERBOSE > 0) {
                printf("Intermediate solution with 2opt from %d : ", i);
                for (int j = 0; j < inst -> nnodes; j++) printf("%d", path[j]);
                printf("\nCost: %f\n", cost);
            }

        }

        time = tsp_time_elapsed();
        
        if (cost < inst -> tsp_best_cost - TSP_EPSYLON) {   //if this solution is better than the best seen so far update it
            if (TSP_VERBOSE > 0) printf("New best solution\n");
            inst -> tsp_best_cost = cost;
            tsp_update_best_sol(inst, path);
            inst -> tsp_best_time = time;
        }

        if (tsp_time_elapsed() > tsp_time_limit) { free(path); return -1; } //if I exceeded the time limit

    }

    free(path);

    return 0;

}


void tsp_check_integrity(const tsp_instance* inst, double cost, int* path) {  //debugging tool

    double c_cost = 0;
    int first = path[0], error = 0;
    int check_feasibility[inst -> nnodes];

    for (int i = 0; i < inst -> nnodes; i++) check_feasibility[i] = 0;

    for (int i = 1; i < inst -> nnodes; i++) {
        if (check_feasibility[path[i]] == 1) { error = 1; break; }
        check_feasibility[path[i]] += 1;
        c_cost += inst -> costs[path[i-1]][path[i]];
    }

    c_cost += inst -> costs[path[inst -> nnodes - 1]][path[0]];

    if (fabs(c_cost - cost) > TSP_EPSYLON) error = 1;

    if (error == 1) {
        printf("INTEGRITY COMPROMISED\nSolution considered: ");
        for(int i = 0; i < inst -> nnodes; i++) {
            printf("%d ", path[i]);
        }
        printf("\nCost: %f, Checked cost: %f\n", cost, c_cost);
        exit(1);
    }

}


double tsp_greedy_from_node(const tsp_instance* inst, int* path, int start_node) {    //greedy solution starting from specific node

    int frontier = start_node, next = -1;
    double cost = 0;
    int visited[inst -> nnodes];

    for (int i = 0; i < inst -> nnodes; i++) visited[i] = 0;

    path[0] = start_node;
    visited[start_node] = 1;

    for (int i = 1; i < inst -> nnodes; i++) {  //repeat nnodes times

        for (int j = 0; j < inst -> nnodes - 1; j++) {

            next = inst -> min_edges[frontier][j]; //check the min_edges in the tsp_instance struct for more info
            if (!visited[next]) break;    //if I didn't explore that node yet then it's the closest (new) node

        }
        
        path[i] = next;
        visited[next] = 1;
        cost += inst -> costs[frontier][next];
        frontier = next;

    }

    cost += inst -> costs[frontier][start_node];    //add the cost of the last edge

    return cost;

}

double tsp_2opt(const tsp_instance* inst, int* path, double cost) {   //2opt algorithm

    double new_cost=cost, check_cost = 1;

    while (check_cost > 0) {   //repeat until I can't find new swaps that improve the cost of my solution

        check_cost = tsp_find_2optswap(inst, path, new_cost);
        if (check_cost > 0) new_cost = check_cost;

    }

    return new_cost;

}

double tsp_find_2optswap(const tsp_instance* inst, int* path, double cost) {  //locating a swap for 2opt algorithm

    double new_cost = cost;
    int k;

    for (int i = 0; i < inst -> nnodes - 1; i++) {
        for (int j = i + 1; j < inst -> nnodes; j++) {
            k = (j+1 == inst -> nnodes) ? 0 : j+1;  //allow for the loop over the edge

            if ((inst -> costs[path[i]][path[i+1]] + inst -> costs[path[j]][path[k]]) > (inst -> costs[path[i]][path[j]] + inst -> costs[path[i+1]][path[k]]) + TSP_EPSYLON){   //if c_i,i+1 + c_j,j+1 > c_i,j + c_i+1,j+1 then swap

                new_cost = new_cost - (inst -> costs[path[i]][path[i+1]] + inst -> costs[path[j]][path[k]]) + (inst -> costs[path[i]][path[j]] + inst -> costs[path[i+1]][path[k]]);    //update the cost

                tsp_reverse(path, i+1, j);  //reverse the part in the middle of the swap

                return new_cost;    //found swap

            }
        }
    }

    return -1;  //no swap found

}

void tsp_reverse(int* path, int start, int end) {   //reverse the array specified betweeb two specified indexes

    int i = start, j = end, c;

    while (i < j) {
        c = path[i];
        path[i] = path[j];
        path[j] = c;
        i++; j--;
    }

}
