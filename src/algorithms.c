#include "../include/algorithms.h"

// TABU
int     tsp_solve_tabu(tsp_instance* inst) {

    int start_node = rand() % inst -> nnodes;
    int* path = calloc(inst -> nnodes, sizeof(int));

    double cost = tsp_greedy_from_node(inst, path, start_node);

    while (tsp_time_elapsed() < tsp_time_limit) tsp_find_2opt_best_swap_tabu(inst, path, &cost);

    if (path != NULL) { free(path); path = NULL; }

    return -1;

}

void tsp_find_2opt_best_swap_tabu(tsp_instance* inst, int* path, double* cost) {

    double best_swap_cost = -INFINITY;
    int best_start = -1, best_end = -1;

    for (int i = 0; i < inst -> nnodes - 1; i++) {

        for (int j = i + 2; j < inst -> nnodes - 1; j++) {
            
            if (tsp_check_tabu(i, j)) continue;

            double difference = (inst -> costs[path[i] * inst -> nnodes + path[i+1]] + inst -> costs[path[j] * inst -> nnodes + path[j+1]]) - (inst -> costs[path[i] * inst -> nnodes + path[j]] + inst -> costs[path[i+1] * inst -> nnodes + path[j+1]]);

            if (difference > best_swap_cost + TSP_EPSILON) {
                best_swap_cost = difference;
                best_start = i; best_end = j;
            }

        }
    }

    printf("Swap: %15.4f\t\t%4d - %4d\t\t", best_swap_cost, best_start, best_end);
    
    if (best_swap_cost < -TSP_EPSILON)
        tsp_add_tabu(best_start, best_end);
    
    *cost = *cost - best_swap_cost;    //update the cost

    tsp_reverse(path, best_start+1, best_end);  //reverse the part in the middle of the swap

    tsp_check_best_sol(inst, path, *cost, tsp_time_elapsed());

    printf("COST: %15.4f\t\tBEST COST: %15.4f\n", *cost, inst -> best_cost);

}

//GREEDY
int tsp_solve_greedy(tsp_instance* inst, const char g2opt) { //solve using greedy (+ 2opt algorithm if g2opt==1)

    int* path = (int*)calloc(inst -> nnodes, sizeof(int));

    for (int i = 0; i < inst -> nnodes; i++) {  //for each starting node
        
        double cost = tsp_greedy_from_node(inst, path, i);

        #if TSP_VERBOSE >= 100
        tsp_check_integrity(inst, cost, path);
        #endif

        #if TSP_VERBOSE >= 10
        printf("Greedy solution from %4d:\t%15.4f\t\t", i, cost);
        #endif

        if (g2opt) {   //if I want to use the 2opt

            int (*swap_function)(const tsp_instance*, int*, double*);
            if (!strcmp(tsp_alg_type, "g2opt")) swap_function = tsp_find_2opt_swap;
            if (!strcmp(tsp_alg_type, "g2opt-best")) swap_function = tsp_find_2opt_best_swap;

            tsp_2opt(inst, path, &cost, swap_function);  //fix solution using 2opt

            #if TSP_VERBOSE >= 10
            printf("(%s) ->\t%15.4f", tsp_alg_type+1, cost);
            #endif

            #if TSP_VERBOSE >= 100
            tsp_check_integrity(inst, cost, path);
            #endif

        }

        double elapsed_time = tsp_time_elapsed();

        #if TSP_VERBOSE >= 10
        printf("\t\t/\tTime elapsed: %15.4f\n", elapsed_time);
        #endif

        #if TSP_VERBOSE >= 100
        printf("Integrity check passed.\n");    //if I get here I've passed the integrity check
        #endif
        
        tsp_check_best_sol(inst, path, cost, elapsed_time); //if this solution is better than the best seen so far update it

        if (elapsed_time > tsp_time_limit) { if(path != NULL) free(path); path = NULL; return -1; } //if I exceeded the time limit

    }

    if (path != NULL) { free(path); path = NULL; }

    return 0;

}

int tsp_solve_greedy_mt(tsp_instance* inst, const char g2opt) { //solve using multithread greedy (+ 2opt if g2opt == 1)

    int (*swap_function)(const tsp_instance*, int*, double*);
    if (!strcmp(tsp_alg_type, "g2opt")) swap_function = tsp_find_2opt_swap;
    if (!strcmp(tsp_alg_type, "g2opt-best")) swap_function = tsp_find_2opt_best_swap;

    tsp_mt_parameters params[N_THREADS];

    for (int i = 0; i < inst -> nnodes; i++) {

        if (tsp_time_elapsed() > tsp_time_limit) {
            tsp_wait_all_threads();
            return -1;
        }

        int thread = tsp_wait_for_thread();

        params[thread].t_index = thread;
        params[thread].inst = inst;
        params[thread].alg = g2opt;
        params[thread].s_node = i;
        params[thread].swap_function = swap_function;

        pthread_create(&tsp_threads[thread], NULL, tsp_greedy_from_node_mt, (void*)&params[thread]);

    }

    tsp_wait_all_threads();

    return 0;

}

void* tsp_greedy_from_node_mt(void* params) {

    int* path = (int*)calloc(((tsp_mt_parameters*)params)->inst -> nnodes, sizeof(int));
        
    double cost = tsp_greedy_from_node(((tsp_mt_parameters*)params)->inst, path, ((tsp_mt_parameters*)params)->s_node); //greedy starting from specified starting node

    #if TSP_VERBOSE >= 100
    tsp_check_integrity(((tsp_mt_parameters*)params)->inst, cost, path);
    #endif

    #if TSP_VERBOSE >= 10
    printf("Greedy solution from %4d:\t%15.4f\t\t", ((tsp_mt_parameters*)params)->s_node, cost);
    #endif

    if (((tsp_mt_parameters*)params)->alg) {   //if I want to use the 2opt

        tsp_2opt(((tsp_mt_parameters*)params)->inst, path, &cost, ((tsp_mt_parameters*)params)->swap_function);  //fix solution using 2opt

        #if TSP_VERBOSE >= 10
        printf("(%s) ->\t%15.4f", tsp_alg_type+1, cost);
        #endif

        #if TSP_VERBOSE >= 100
        tsp_check_integrity(((tsp_mt_parameters*)params)->inst, cost, path);
        #endif

    }

    double elapsed_time = tsp_time_elapsed();

    #if TSP_VERBOSE >= 10
    printf("\t\t/\tTime elapsed: %15.4f\n", elapsed_time);
    #endif

    #if TSP_VERBOSE >= 100
    printf("Integrity check passed.\n");    //if I get here I've passed the integrity check
    #endif

    tsp_check_best_sol(((tsp_mt_parameters*)params)->inst, path, cost, elapsed_time); //if this solution is better than the best seen so far update it

    if(path != NULL) { free(path); path = NULL; }

    tsp_free_thread(((tsp_mt_parameters*)params) -> t_index);

}

double tsp_greedy_from_node(const tsp_instance* inst, int* path, int start_node) { //greedy solution starting from specific node

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

    if (visited != NULL) { free(visited); visited = NULL; }

    cost += inst -> costs[frontier * inst -> nnodes + start_node];    //add the cost of the last edge

    return cost;

}


//_2OPT
void tsp_2opt(const tsp_instance* inst, int* path, double* cost, int (*swap_function)(const tsp_instance*, int*, double*)) { //2opt algorithm

    #if TSP_VERBOSE >= 100

    int counter = 0;
    while ((*swap_function)(inst, path, cost) > 0) { //repeat until I can't find new swaps that improve the cost of my solution
        counter++;
    }
    printf("- %d swaps", counter);

    #else

    while ((*swap_function)(inst, path, cost) > 0); //repeat until I can't find new swaps that improve the cost of my solution

    #endif

}

int tsp_find_2opt_swap(const tsp_instance* inst, int* path, double* cost) { //locating a swap for 2opt algorithm

    for (int i = 0; i < inst -> nnodes - 1; i++) {
        for (int j = i + 1; j < inst -> nnodes; j++) {
            int k = (j+1 == inst -> nnodes) ? 0 : j+1;  //allow for the loop over the edge

            double difference = (inst -> costs[path[i] * inst -> nnodes + path[i+1]] + inst -> costs[path[j] * inst -> nnodes + path[k]]) - (inst -> costs[path[i] * inst -> nnodes + path[j]] + inst -> costs[path[i+1] * inst -> nnodes + path[k]]);

            if (difference > TSP_EPSILON) {
                *cost = *cost - difference;    //update the cost
                tsp_reverse(path, i+1, j);  //reverse the part in the middle of the swap
                return 1;   //found swap
            }

        }
    }

    return -1; //no swap found

}

int tsp_find_2opt_best_swap(const tsp_instance* inst, int* path, double* cost) { //locating the best swap for 2opt algorithm

    double best_swap_cost = 0;
    int best_start = -1, best_end = -1;

    for (int i = 0; i < inst -> nnodes - 1; i++) {
        for (int j = i + 1; j < inst -> nnodes; j++) {
            int k = (j+1 == inst -> nnodes) ? 0 : j+1;  //allow for the loop over the edge

            double difference = (inst -> costs[path[i] * inst -> nnodes + path[i+1]] + inst -> costs[path[j] * inst -> nnodes + path[k]]) - (inst -> costs[path[i] * inst -> nnodes + path[j]] + inst -> costs[path[i+1] * inst -> nnodes + path[k]]);

            if (difference > TSP_EPSILON && difference > best_swap_cost + TSP_EPSILON) {
                best_swap_cost = difference;
                best_start = i; best_end = j;
            }

        }
    }

    if (best_swap_cost < -TSP_EPSILON) return -1; //no swap found

    *cost = *cost - best_swap_cost;    //update the cost

    tsp_reverse(path, best_start+1, best_end);  //reverse the part in the middle of the swap

    return 1;    //found swap

}
