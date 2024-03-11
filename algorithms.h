#include "tsp.h"

//ALGORITHMS
//int     tsp_solve_greedy(tsp_instance* inst);
int     tsp_solve_g2opt(tsp_instance* inst, const char g2opt);

void    tsp_check_integrity(const tsp_instance* inst, double cost, int* path);

double  tsp_greedy_from_node(const tsp_instance* inst, int* path, int start_node);
int     tsp_find_min_edge(const tsp_instance* inst, int frontier, int* path);
double  tsp_2opt(const tsp_instance* inst, int* path, double cost);
double  tsp_find_2optswap(const tsp_instance* inst, int* path, double cost);
void    tsp_reverse(int* path, int start, int end);


/*int tsp_solve_greedy(tsp_instance* inst) {   //solve using greedy algorithm

    int* path = calloc(inst -> nnodes, sizeof(int));

    for (int i = 0; i < inst -> nnodes; i++) {  //for each starting node
        
        double cost = tsp_greedy_from_node(inst, path, i);
        if (tsp_verbose > 0) tsp_check_integrity(inst, cost, path);

        if (tsp_verbose > 0) {
            printf("\nIntermediate solution from %d: ", i);
            for (int j = 0; j < inst -> nnodes; j++) printf("%d", path[j]);
            printf("\nCost: %f\n", cost);
        }

        double time = tsp_time_elapsed();
        
        if (cost < inst -> tsp_best_cost - TSP_EPSYLON) {   //if this solution is better than the best seen so far update it
            if (tsp_verbose > 0) printf("New best solution\n");
            inst -> tsp_best_cost = cost;
            tsp_update_best_sol(inst, path);
            inst -> tsp_best_time = time;
        }

        if (tsp_time_elapsed() > tsp_time_limit) { free(path); return -1; } //if I exceeded the time limit

    }

    free(path);

    return 0;

}*/

int tsp_solve_g2opt(tsp_instance* inst, const char g2opt) {    //solve using greedy (+ 2opt algorithm if g2opt==1)

    double time;
    int* path = calloc(inst -> nnodes, sizeof(int));

    for (int i = 0; i < inst -> nnodes; i++) {  //for each starting node
        
        double cost = tsp_greedy_from_node(inst, path, i);
        if (tsp_verbose > 0) tsp_check_integrity(inst, cost, path);

        if (tsp_verbose > 0) {
            printf("\nIntermediate solution from %d: ", i);
            for (int j = 0; j < inst -> nnodes; j++) printf("%d", path[j]);
            printf("\nCost: %f\n", cost);
        }

        if (g2opt == 1) {
            cost = tsp_2opt(inst, path, cost);  //fix solution using 2opt
            if (tsp_verbose > 0) tsp_check_integrity(inst, cost, path);

            if (tsp_verbose > 0) {
                printf("Intermediate solution with 2opt from %d: ", i);
                for (int j = 0; j < inst -> nnodes; j++) printf("%d", path[j]);
                printf("\nCost: %f\n", cost);
                printf("Time: %f\n\n", time);
            }
        }

        time = tsp_time_elapsed();
        
        if (cost < inst -> tsp_best_cost - TSP_EPSYLON) {   //if this solution is better than the best seen so far update it
            if (tsp_verbose > 0) printf("New best solution\n");
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

    for (int i = 0; i < inst -> nnodes; i++) path[i] = -1;

    int frontier = start_node, next = -1;
    double cost = 0;

    for (int i = 0; i < inst -> nnodes - 1; i++) {  //repeat nnodes times

        next = tsp_find_min_edge(inst, frontier, path); //get nearest node
        path[frontier] = next;  //path[3] = 1 means that after the node 3, we must read the node in the position 1 (each element points to where to find the next node)
        cost += inst -> costs[frontier][next];
        frontier = next;

    }

    path[frontier] = start_node;    //close the loop
    cost += inst -> costs[frontier][start_node];    //add the cost of the last edge

    int t_path[inst -> nnodes];     //go from the path[node] = successor's_location representation to the path[i] = ith node
    frontier = start_node, next = path[frontier];
    for (int i = 0; i < inst -> nnodes; i++) {
        t_path[i] = frontier;
        frontier = next;
        next = path[next];
    }

    for (int i = 0; i < inst -> nnodes; i++) path[i] = t_path[i];   //save the new representation in the original array

    return cost;

}

int tsp_find_min_edge(const tsp_instance* inst, int node, int* path) {    //closest node to specific node

    int min_index = -1;

    for (int i = 0; i < inst -> nnodes -1; i++) {

        min_index = inst -> min_edges[node][i]; //check the min_edges in the tsp_instance struct
        if (path[min_index] == -1) return min_index;    //if I didn't explore that node yet then it's the closest (new) node

    }

    printf("\nMIN_EDGE ERROR\n");
    exit(1);

}

double tsp_2opt(const tsp_instance* inst, int* path, double cost) {   //2opt algorithm

    double new_cost=cost, tmp = 1;

    while (tmp > 0) {   //till I can't find new paths keep checking for possible swaps

        tmp = tsp_find_2optswap(inst, path, new_cost);
        if (tmp > 0) new_cost = tmp;

    }

    return new_cost;

}

double tsp_find_2optswap(const tsp_instance* inst, int* path, double cost) {  //locating a swap for 2opt algorithm

    double new_cost = cost;
    int k;

    for (int i = 0; i < inst -> nnodes - 1; i++) {
        for (int j = i + 1; j < inst -> nnodes; j++) {
            k = (j+1 == inst -> nnodes) ? 0 : j+1;

            if ((inst -> costs[path[i]][path[i+1]] + inst -> costs[path[j]][path[k]]) > (inst -> costs[path[i]][path[j]] + inst -> costs[path[i+1]][path[k]]) + TSP_EPSYLON){   //if c_i,i+1 + c_j,j+1 > c_i,j + c_i+1,j+1 then swap

                new_cost = new_cost - (inst -> costs[path[i]][path[i+1]] + inst -> costs[path[j]][path[k]]) + (inst -> costs[path[i]][path[j]] + inst -> costs[path[i+1]][path[k]]);    //update the cost

                tsp_reverse(path, i+1, j);  //reverse the part in the middle of the swap

                return new_cost;

            }
        }
    }

    return -1;

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