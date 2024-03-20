#include "../include/algorithms.h"


//2OPT
/**
 * Applies a swap policy till the solution cannot be improved further
 * 
 * @param path Path considered (will be changed if it finds a swap)
 * @param cost Cost of the current path (will be changed if it finds a swap)
 * @param swap_function The swap function to use
*/
void tsp_2opt(int* path, double* cost, int (*swap_function)(int*, double*)) {

    while ((*swap_function)(path, cost) > 0); //repeat until I can't find new swaps that improve the cost of my solution

}

int tsp_find_2opt_swap(int* path, double* cost) {

    for (int i = 0; i < inst.nnodes - 2; i++) {
        for (int j = i + 2; j < inst.nnodes; j++) {
            if (i == 0 && j+1 == inst.nnodes) continue;
            int k = (j+1 == inst.nnodes) ? 0 : j+1;  //allow for the loop over the edge

            double improvement = (inst.costs[path[i] * inst.nnodes + path[i+1]] + inst.costs[path[j] * inst.nnodes + path[k]]) - (inst.costs[path[i] * inst.nnodes + path[j]] + inst.costs[path[i+1] * inst.nnodes + path[k]]);

            if (improvement > TSP_EPSILON) {
                *cost = *cost - improvement;    //update the cost
                tsp_reverse(path, i+1, j);  //reverse the part in the middle of the swap
                return 1;   //found swap
            }

        }
    }

    return -1; //no swap found

}

int tsp_find_2opt_best_swap(int* path, double* cost) {

    double best_improvement = -INFINITY;
    int best_start = -1, best_end = -1;

    for (int i = 0; i < inst.nnodes - 2; i++) {        
        for (int j = i + 2; j < inst.nnodes; j++) {
            if (i == 0 && j+1 == inst.nnodes) continue;
            int k = (j+1 == inst.nnodes) ? 0 : j+1;  //allow for the loop over the edge

            double improvement = (inst.costs[path[i] * inst.nnodes + path[i+1]] + inst.costs[path[j] * inst.nnodes + path[k]]) - (inst.costs[path[i] * inst.nnodes + path[j]] + inst.costs[path[i+1] * inst.nnodes + path[k]]);

            if (improvement > TSP_EPSILON && improvement > best_improvement + TSP_EPSILON) {
                best_improvement = improvement;
                best_start = i; best_end = j;
            }

        }
    }

    if (best_improvement < -TSP_EPSILON) return -1; //no swap found

    *cost = *cost - best_improvement;    //update the cost

    tsp_reverse(path, best_start+1, best_end);  //reverse the part in the middle of the swap

    return 1;    //found swap

}


//GREEDY
/**
 * Finds the greedy solution starting from a specified node and saves it in path
 * 
 * @param path Path considered (the greedy sol will be saved here)
 * @param start_node Starting node for this solution
 * 
 * @returns The cost of the created path
*/
double tsp_greedy_path_from_node(int* path, const int start_node) {

    int frontier = start_node, next = -1;
    double cost = 0;
    int* visited = (int*)calloc(inst.nnodes, sizeof(int));

    path[0] = start_node;
    visited[start_node] = 1;

    for (int i = 1; i < inst.nnodes; i++) {  //repeat nnodes times

        for (int j = 0; j < inst.nnodes - 1; j++) {

            next = inst.sort_edges[frontier * (inst.nnodes - 1) + j]; //check the min_edges in the tsp_instance struct for more info
            if (!visited[next]) break;    //if I didn't explore that node yet then it's the closest (new) node

        }
        
        path[i] = next;
        visited[next] = 1;
        cost += inst.costs[frontier * inst.nnodes + next];
        frontier = next;

    }

    if (visited != NULL) { free(visited); visited = NULL; }

    cost += inst.costs[frontier * inst.nnodes + start_node];    //add the cost of the last edge

    return cost;

}

/**
 * (THREAD SPECIFIC) Finds the greedy solution starting from a specified node and applies (if specified) a 2opt optimization
 * 
 * @param params tsp_mt_parameters type structture using: ->t_index, ->inst, ->s_node, ->swap_function
*/
void* tsp_greedy_from_node(void* params) {

    int* path = (int*)calloc(inst.nnodes, sizeof(int));
    double cost = tsp_greedy_path_from_node(path, ((tsp_mt_parameters*)params)->s_node); //greedy starting from specified starting node

    #if TSP_VERBOSE >= 100
    tsp_check_integrity(path, cost);
    #endif

    if (((tsp_mt_parameters*)params)->swap_function != NULL) {   //if I want to use the 2opt

        tsp_2opt(path, &cost, ((tsp_mt_parameters*)params)->swap_function);  //fix solution using 2opt

        #if TSP_VERBOSE >= 100
        tsp_check_integrity(path, cost);
        #endif

    }

    tsp_check_best_sol(path, cost, tsp_time_elapsed()); //if this solution is better than the best seen so far update it

    if(path != NULL) { free(path); path = NULL; }

    tsp_free_thread(((tsp_mt_parameters*)params)->t_index);

}

int tsp_solve_greedy(int (*swap_function)(int*, double*)) {

    tsp_mt_parameters params[N_THREADS];

    for (int i = 0; i < inst.nnodes; i++) {

        if (tsp_time_elapsed() > tsp_time_limit) {
            tsp_wait_all_threads();
            return -1;
        }

        int thread = tsp_wait_for_thread();

        params[thread].t_index = thread;
        params[thread].s_node = i;
        params[thread].swap_function = swap_function;

        pthread_create(&tsp_threads[thread], NULL, tsp_greedy_from_node, (void*)&params[thread]);

    }

    tsp_wait_all_threads();

    return 0;

}

int tsp_solve_greedy_st(int (*swap_function)(int*, double*)) {

    int* path = (int*)calloc(inst.nnodes, sizeof(int));

    for (int i = 0; i < inst.nnodes; i++) {  //for each starting node
        
        double cost = tsp_greedy_path_from_node(path, i);

        #if TSP_VERBOSE >= 100
        tsp_check_integrity(path, cost);
        #endif

        #if TSP_VERBOSE >= 10
        printf("Greedy solution from %4d:\t%15.4f\t\t", i, cost);
        #endif

        if (swap_function != NULL) {   //if I want to use the 2opt

            tsp_2opt(path, &cost, swap_function);  //fix solution using 2opt

            #if TSP_VERBOSE >= 10
            printf("(%s) ->\t%15.4f", tsp_alg_type+1, cost);
            #endif

            #if TSP_VERBOSE >= 100
            tsp_check_integrity(path, cost);
            #endif

        }

        double elapsed_time = tsp_time_elapsed();

        #if TSP_VERBOSE >= 10
        printf("\t\t/\tTime elapsed: %15.4f\n", elapsed_time);
        #endif

        #if TSP_VERBOSE >= 100
        printf("Integrity check passed.\n");    //if I get here I've passed the integrity check
        #endif
        
        tsp_check_best_sol(path, cost, elapsed_time); //if this solution is better than the best seen so far update it

        if (elapsed_time > tsp_time_limit) { if(path != NULL) { free(path); path = NULL; } return -1; } //if I exceeded the time limit

    }

    if (path != NULL) { free(path); path = NULL; }

    return 0;

}


// TABU
/**
 * Finds the best 2opt swap checking and updating the tabu list after performing it
 * 
 * @param path Path considered (will be changed if it finds a swap)
 * @param cost Cost of the current path (will be changed if it finds a swap)
 * @param t_index Index of the thread working on this task
*/
void tsp_find_tabu_swap(int* path, double* cost, const int t_index) {

    double best_improvement = -INFINITY;
    int best_start = -1, best_end = -1;

    for (int i = 0; i < inst.nnodes - 2; i++) {

        if (tsp_check_tabu(t_index, path[i], path[i+1])) continue;

        for (int j = i + 2; j < inst.nnodes; j++) {
            if (i == 0 && j+1 == inst.nnodes) continue;
            int k = (j+1 == inst.nnodes) ? 0 : j+1;  //allow for the loop over the edge

            if (tsp_check_tabu(t_index, path[j], path[k])) continue;

            double improvement = (inst.costs[path[i] * inst.nnodes + path[i+1]] + inst.costs[path[j] * inst.nnodes + path[k]]) - (inst.costs[path[i] * inst.nnodes + path[j]] + inst.costs[path[i+1] * inst.nnodes + path[k]]);
            
            if (improvement > best_improvement + TSP_EPSILON) {
                best_improvement = improvement;
                best_start = i; best_end = j;
            }

        }
    }

    if (best_improvement < -TSP_EPSILON)
        tsp_add_tabu(t_index, path[best_start], path[best_end]);

    *cost = *cost - best_improvement;

    tsp_reverse(path, best_start+1, best_end);

}

/**
 * (THREAD SPECIFIC) Looks for 2opt swaps till there's time left
 * 
 * @param params tsp_mt_parameters type structture using: ->t_index, ->inst, ->path, ->cost
*/
void* tsp_tabu_from_node(void* params) {

    int* path = ((tsp_mt_parameters*)params)->path;
    double* cost = ((tsp_mt_parameters*)params)->cost;
    int t_index = ((tsp_mt_parameters*)params)->t_index;
    double time = 0;

    while (time < tsp_time_limit) {

        tsp_find_tabu_swap(path, cost, t_index);

        time = tsp_time_elapsed();

        #if TSP_VERBOSE >= 100
        tsp_check_integrity(path, *cost);
        #endif

        tsp_check_best_sol(path, *cost, time);

    }

    tsp_free_thread(t_index);

}

int tsp_solve_tabu() {
    
    tsp_mt_parameters params[N_THREADS];
    double cost[N_THREADS];

    for (int i = 0; i < N_THREADS; i++) {

        int thread = tsp_wait_for_thread();

        int start_node = rand() % inst.nnodes;
        int* path = calloc(inst.nnodes, sizeof(int));

        cost[thread] = tsp_greedy_path_from_node(path, start_node);

        params[thread].t_index = thread;
        params[thread].path = path;
        params[thread].cost = &cost[thread];

        pthread_create(&tsp_threads[thread], NULL, tsp_tabu_from_node, (void*)&params[thread]);

    }

    tsp_wait_all_threads();

    for (int thread = 0; thread < N_THREADS; thread++)
        if (params[thread].path != NULL) { free(params[thread].path); params[thread].path = NULL; }

    return -1;

}