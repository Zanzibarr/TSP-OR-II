#include "tsp.h"

//ALGORITHMS
int     tsp_solve_greedy(tsp_instance* inst);
int     tsp_solve_g2opt(tsp_instance* inst);

void    tsp_check_integrity(const tsp_instance* inst, double cost, int* path);

double  tsp_greedy_from_node(const tsp_instance* inst, int* path, int start_node);
int     tsp_find_min_edge(const tsp_instance* inst, int frontier, int* path);
double  tsp_2opt(const tsp_instance* inst, int* path, double cost);
double  tsp_find_2optswap(const tsp_instance* inst, int* path, double cost);
void    tsp_reverse(int* path, int start, int end);


int tsp_solve_greedy(tsp_instance* inst) {   //solve using greedy algorithm

    int i, j, starting_node, current_node, next_node;
    double best_edge_cost, current_edge_cost, current_path_cost, best_path_cost;
    time_t *starting_time, *current_iteration_time, best_time, total_time;
    int *touched_nodes, *current_sol;

    //starting_time = malloc(sizeof(time_t));
    //current_iteration_time = malloc(sizeof(time_t));
    //time(starting_time);
    //time(current_iteration_time);
    touched_nodes = calloc(inst->nnodes, sizeof(int));  //Ricordati di liberare la memoria.
    current_sol = calloc(inst->nnodes, sizeof(int));    //Ricordati di liberare la memoria.
    total_time = 0;
    starting_node = 0;
    best_path_cost = INFINITY;

    while (starting_node<inst->nnodes) {
        printf("-------%d-------\n", starting_node);
        for (i=0; i<inst->nnodes; i++) touched_nodes[i]=0;
        current_path_cost = 0;
        current_sol[0] = starting_node;
        touched_nodes[starting_node] = 1;
        //time(current_iteration_time);
        current_node = starting_node;

        for (i=1; i<inst->nnodes; i++) {
            next_node = -1;
            best_edge_cost = INFINITY;
            for (j=0; j<inst->nnodes; j++) {
                //if (tsp_get_elapsed_time(starting_time)==-1) return;
                if (j==current_node) continue;
                if (touched_nodes[j]==1) continue;
                current_edge_cost = inst->costs[current_node][j];
                if (current_edge_cost<best_edge_cost) {
                    best_edge_cost = current_edge_cost;
                    next_node = j;
                }
            }
            current_sol[i] = next_node;
            touched_nodes[next_node] = 1;
            current_path_cost += inst->costs[current_node][next_node];
            current_node = next_node; 
        }

        for (int i = 0; i < inst -> nnodes; i++) printf("%d", current_sol[i]);
        printf("\n      Cost: %f\n", current_path_cost);

        current_path_cost += inst -> costs[current_sol[inst -> nnodes - 1]][current_sol[0]];

        printf("Fixed Cost: %f\n", current_path_cost);

        if (current_path_cost<best_path_cost) {
            best_path_cost = current_path_cost;
            inst->tsp_best_time = tsp_time_elapsed();
            for (j=0; j<inst->nnodes; j++) inst->tsp_best_solution[j] = current_sol[j];
            inst->tsp_best_cost = best_path_cost;
        }

        starting_node++;

        if (tsp_time_elapsed() > tsp_time_limit) { return -1; }

    }

    return 0;

}

int tsp_solve_g2opt(tsp_instance* inst) {    //solve using greedy + 2opt algorithm

    int* path = calloc(inst -> nnodes, sizeof(int));

    for (int i = 0; i < inst -> nnodes; i++) {
        
        double cost = tsp_greedy_from_node(inst, path, i);
        if (tsp_verbose > 0) tsp_check_integrity(inst, cost, path);

        if (tsp_verbose > 0) {
            printf("\nIntermediate solution from %d: ", i);
            for (int j = 0; j < inst -> nnodes; j++) printf("%d", path[j]);
            printf("\nCost: %f\n", cost);
        }

        cost = tsp_2opt(inst, path, cost);
        if (tsp_verbose > 0) tsp_check_integrity(inst, cost, path);

        double time = tsp_time_elapsed();

        if (tsp_verbose > 0) {
            printf("Intermediate solution with 2opt from %d: ", i);
            for (int j = 0; j < inst -> nnodes; j++) printf("%d", path[j]);
            printf("\nCost: %f\n", cost);
            printf("Time: %f\n\n", time);
        }
        
        if (cost < inst -> tsp_best_cost - TSP_EPSYLON) {
            if (tsp_verbose > 0) printf("New best solution\n");
            inst -> tsp_best_cost = cost;
            tsp_update_best_sol(inst, path);
            inst -> tsp_best_time = time;
        }

        if (tsp_time_elapsed() > tsp_time_limit) { free(path); return -1; }

    }

    free(path);

    return 0;

}


void tsp_check_integrity(const tsp_instance* inst, double cost, int* path) {  //debugging tool

    double c_cost = 0;
    int first = path[0], error = 0;

    for (int i = 1; i < inst -> nnodes; i++) {
        if (path[i] == first) { error = 1; break; }
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

    for (int i = 0; i < inst -> nnodes - 1; i++) {

        next = tsp_find_min_edge(inst, frontier, path);
        path[frontier] = next;
        cost += inst -> costs[frontier][next];
        frontier = next;

    }

    path[frontier] = start_node;
    cost += inst -> costs[frontier][start_node];

    int t_path[inst -> nnodes];
    frontier = start_node, next = path[frontier];
    for (int i = 0; i < inst -> nnodes; i++) {
        t_path[i] = frontier;
        frontier = next;
        next = path[next];
    }

    for (int i = 0; i < inst -> nnodes; i++) path[i] = t_path[i];

    return cost;

}

int tsp_find_min_edge(const tsp_instance* inst, int node, int* path) {    //closest node to specific node

    int min_index = -1;

    for (int i = 0; i < inst -> nnodes -1; i++) {

        min_index = inst -> min_edges[node][i];
        if (path[min_index] == -1) return min_index;

    }

    printf("\nMIN_EDGE ERROR\n");
    exit(1);

}

double tsp_2opt(const tsp_instance* inst, int* path, double cost) {   //2opt algorithm

    double new_cost=cost, tmp = 1;

    while (tmp > 0) {

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

            if ((inst -> costs[path[i]][path[i+1]] + inst -> costs[path[j]][path[k]]) - (inst -> costs[path[i]][path[j]] + inst -> costs[path[i+1]][path[k]]) > TSP_EPSYLON){

                new_cost = new_cost - (inst -> costs[path[i]][path[i+1]] + inst -> costs[path[j]][path[k]]) + (inst -> costs[path[i]][path[j]] + inst -> costs[path[i+1]][path[k]]);

                tsp_reverse(path, i+1, j);

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