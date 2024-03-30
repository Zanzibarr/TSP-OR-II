#include "../include/algorithms.h"


//2OPT
/**
 * @brief Applies a swap policy till the solution cannot be improved further
 * 
 * @param path Path considered (will be changed if it finds a swap)
 * @param cost Cost of the current path (will be changed if it finds a swap)
 * @param swap_function The swap function to use
*/
void tsp_2opt(int* path, double* cost, int (*swap_function)(int*, double*)) {

    while ((*swap_function)(path, cost) > 0 && tsp_time_elapsed() < tsp_time_limit) { //repeat until I can't find new swaps that improve the cost of my solution

        #if TSP_VERBOSE >= 50
        tsp_save_intermediate_cost(0, *cost);
        tsp_save_intermediate_cost(1, tsp_inst.best_cost);
        #endif
    
    }

}

int tsp_find_2opt_swap(int* path, double* cost) {

    for (int i = 0; i < tsp_inst.nnodes - 2; i++) {
        for (int j = i + 2; j < tsp_inst.nnodes; j++) {
            if (i == 0 && j+1 == tsp_inst.nnodes) continue;
            int k = (j+1 == tsp_inst.nnodes) ? 0 : j+1;  //allow for the loop over the edge

            double improvement =    (tsp_inst.costs[path[i] * tsp_inst.nnodes + path[i+1]] + tsp_inst.costs[path[j]   * tsp_inst.nnodes + path[k]]) -
                                    (tsp_inst.costs[path[i] * tsp_inst.nnodes + path[j]]   + tsp_inst.costs[path[i+1] * tsp_inst.nnodes + path[k]]);

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

    for (int i = 0; i < tsp_inst.nnodes - 2; i++) {        
        for (int j = i + 2; j < tsp_inst.nnodes; j++) {
            if (i == 0 && j+1 == tsp_inst.nnodes) continue;
            int k = (j+1 == tsp_inst.nnodes) ? 0 : j+1;  //allow for the loop over the edge

            double improvement =    (tsp_inst.costs[path[i] * tsp_inst.nnodes + path[i+1]] + tsp_inst.costs[path[j]   * tsp_inst.nnodes + path[k]]) -
                                    (tsp_inst.costs[path[i] * tsp_inst.nnodes + path[j]]   + tsp_inst.costs[path[i+1] * tsp_inst.nnodes + path[k]]);

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
 * @brief Finds the greedy solution starting from a specified node and saves it in path
 * 
 * @param path Path considered (the greedy sol will be saved here)
 * @param start_node Starting node for this solution
 * 
 * @return The cost of the created path
*/
double tsp_greedy_path_from_node(int* path, const int start_node) {

    int frontier = start_node, next = -1;
    double cost = 0;
    int* visited = (int*)calloc(tsp_inst.nnodes, sizeof(int));

    path[0] = start_node;
    visited[start_node] = 1;

    for (int i = 1; i < tsp_inst.nnodes; i++) {  //repeat nnodes times

        for (int j = 0; j < tsp_inst.nnodes - 1; j++) {

            next = tsp_inst.sort_edges[frontier * (tsp_inst.nnodes - 1) + j]; //check the min_edges in the tsp_instance struct for more info
            if (!visited[next]) break;    //if I didn't explore that node yet then it's the closest (new) node

        }
        
        path[i] = next;
        visited[next] = 1;
        cost += tsp_inst.costs[frontier * tsp_inst.nnodes + next];
        frontier = next;

    }

    if (visited != NULL) { free(visited); visited = NULL; }

    cost += tsp_inst.costs[frontier * tsp_inst.nnodes + start_node];    //add the cost of the last edge

    return cost;

}

/**
 * @brief (THREAD SPECIFIC) Finds the greedy solution starting from a specified node and applies (if specified) a 2opt optimization
 * 
 * @param params tsp_mt_parameters type structure using: ->t_index, ->s_node, ->swap_function
*/
void* tsp_greedy_from_node(void* params) {

    int* path = (int*)calloc(tsp_inst.nnodes, sizeof(int));
    double cost = tsp_greedy_path_from_node(path, ((tsp_mt_parameters*)params)->s_node); //greedy starting from specified starting node

    #if TSP_VERBOSE >= 100
    tsp_check_integrity(path, cost, "algorithms.c: tsp_greedy_from_node - 1");
    #endif

    if (((tsp_mt_parameters*)params)->swap_function != NULL) {   //if I want to use the 2opt

        tsp_2opt(path, &cost, ((tsp_mt_parameters*)params)->swap_function);  //fix solution using 2opt

        #if TSP_VERBOSE >= 100
        tsp_check_integrity(path, cost, "algorithms.c: tsp_greedy_from_node - 2");
        #endif

    }

    tsp_check_best_sol(path, cost, tsp_time_elapsed()); //if this solution is better than the best seen so far update it

    #if TSP_VERBOSE >= 50
    tsp_save_intermediate_cost(0, cost);
    #endif

    if(path != NULL) { free(path); path = NULL; }

    tsp_free_thread(((tsp_mt_parameters*)params)->t_index);

}

int tsp_solve_greedy(int (*swap_function)(int*, double*)) {

    tsp_mt_parameters params[N_THREADS];

    for (int i = 0; i < tsp_inst.nnodes; i++) {

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

    int* path = (int*)calloc(tsp_inst.nnodes, sizeof(int));

    for (int i = 0; i < tsp_inst.nnodes; i++) {  //for each starting node
        
        double cost = tsp_greedy_path_from_node(path, i);

        #if TSP_VERBOSE >= 100
        tsp_check_integrity(path, cost, "algorithms.c: tsp_solve_greedy_st - 1");
        #endif

        if (swap_function != NULL) {   //if I want to use the 2opt

            tsp_2opt(path, &cost, swap_function);  //fix solution using 2opt

            #if TSP_VERBOSE >= 100
            tsp_check_integrity(path, cost, "algorithms.c: tsp_solve_greedy_st - 2");
            #endif

        }

        double elapsed_time = tsp_time_elapsed();
        
        tsp_check_best_sol(path, cost, elapsed_time); //if this solution is better than the best seen so far update it

        #if TSP_VERBOSE >= 50
        tsp_save_intermediate_cost(0, cost);
        #endif

        if (elapsed_time > tsp_time_limit) { if(path != NULL) { free(path); path = NULL; } return -1; } //if I exceeded the time limit

    }

    if (path != NULL) { free(path); path = NULL; }

    return 0;

}


// TABU
/**
 * @brief Finds the best 2opt swap checking and updating the tabu list after performing it
 * 
 * @param path Path considered (will be changed if it finds a swap)
 * @param cost Cost of the current path (will be changed if it finds a swap)
 * @param t_index Index of the thread working on this task
*/
void tsp_find_tabu_swap(int* path, double* cost, const int t_index) {

    double best_improvement = -INFINITY;
    int best_start = -1, best_end = -1;

    for (int i = 0; i < tsp_inst.nnodes - 2; i++) {

        if (tsp_check_tabu(t_index, path[i], path[i+1])) continue;

        for (int j = i + 2; j < tsp_inst.nnodes; j++) {
            if (i == 0 && j+1 == tsp_inst.nnodes) continue;
            int k = (j+1 == tsp_inst.nnodes) ? 0 : j+1;  //allow for the loop over the edge

            if (tsp_check_tabu(t_index, path[j], path[k])) continue;

            double improvement =    (tsp_inst.costs[path[i] * tsp_inst.nnodes + path[i+1]] + tsp_inst.costs[path[j]   * tsp_inst.nnodes + path[k]]) -
                                    (tsp_inst.costs[path[i] * tsp_inst.nnodes + path[j]]   + tsp_inst.costs[path[i+1] * tsp_inst.nnodes + path[k]]);
            
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
 * @brief (THREAD SPECIFIC) Looks for 2opt swaps till there's time left
 * 
 * @param params tsp_mt_parameters type structure using: ->t_index, ->path, ->cost
*/
void* tsp_tabu_from_node(void* params) {

    int* path = ((tsp_mt_parameters*)params)->path;
    double* cost = ((tsp_mt_parameters*)params)->cost;
    int t_index = ((tsp_mt_parameters*)params)->t_index;
    double time = 0;

    while (time < tsp_time_limit) {

        tsp_find_tabu_swap(path, cost, t_index);

        time = tsp_time_elapsed();

        #if TSP_VERBOSE >= 50
        tsp_save_intermediate_cost(t_index, *cost);
        #endif

        #if TSP_VERBOSE >= 100
        tsp_check_integrity(path, *cost, "algorithms.c: tsp_tabu_from_node");
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

        int start_node = rand() % tsp_inst.nnodes;
        int* path = calloc(tsp_inst.nnodes, sizeof(int));

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


// VNS
/**
 * @brief Performs a 3opt kick on the path updating the cost
 * 
 * @param path The path on which to perform the kick
 * @param cost The cost to be updated after the kick
 * @param start The starting index
 * @param end The ending index
 */
void tsp_vns_3kick(int* path, double* cost, int start, int end) {

    if (end - start < 10) return;

    int i = -1, j = -1, k = -1;

    // pick random values that are at least 2 positions apart from each other
    i = (rand() % (end - start - 1)) + start;
    j = i;
    while (j >= i - 1 && j <= i + 1) j = (rand() % (end - start - 1)) + start;
    k = i;
    while (k >= i - 1 && k <= i + 1 || k >= j - 1 && k <= j + 1) k = (rand() % (end - start - 1)) + start;

    // sort them
    if (i > k) { int c = i; i = k; k = c; }
    if (i > j) { int c = i; i = j; j = c; }
    if (j > k) { int c = j; j = k; k = c; }

    // update cost
    *cost = *cost -
        (   //removing the costs of the old edges
            (tsp_inst.costs[path[i] * tsp_inst.nnodes + path[i+1]]) +
            (tsp_inst.costs[path[j] * tsp_inst.nnodes + path[j+1]]) +
            (tsp_inst.costs[path[k] * tsp_inst.nnodes + path[(k+1==tsp_inst.nnodes)?0:k+1]])
        ) +
        (   //adding the costs of the new edges
            (tsp_inst.costs[path[i] * tsp_inst.nnodes + path[k]]) +
            (tsp_inst.costs[path[j+1] * tsp_inst.nnodes + path[i+1]]) +
            (tsp_inst.costs[path[j] * tsp_inst.nnodes + path[(k+1==tsp_inst.nnodes)?0:k+1]])
        );

    // kick
    int* temp_path = (int*)calloc(end-start, sizeof(int));
    for (int c = start;     c <= i;        c++) temp_path[c-start]         = path[c];         //copy start--i
    for (int c = 0;        c < k-j;        c++) temp_path[i+1+c-start]     = path[k-c];       //connect i--k and reverse j+1--k
    for (int c = 0;        c < j-i;        c++) temp_path[i+k-j+1+c-start] = path[i+1+c];     //connect j+1--i+1 and copy i+1--j
    for (int c = k+1;      c < end;        c++) temp_path[c-start]         = path[c];         //connect j--k+1 and copy k+1--end
    for (int c = start; c < end; c++) path[c] = temp_path[c-start];
    if (temp_path != NULL) { free(temp_path); temp_path = NULL; }

}

/**
 * @brief (THREAD SPECIFIC) Performs a variable number of kicks and then fixes the solution using the 2opt algorithm
 * 
 * @param params tsp_mt_parameters type structure using: ->t_index, -> path, -> cost
*/
void* tsp_vns_kicknsolve(void* params) {

    int* path = ((tsp_mt_parameters*)params)->path;
    double* cost = ((tsp_mt_parameters*)params)->cost;
    int t_index = ((tsp_mt_parameters*)params)->t_index;

    for (int i = 0; i < t_index/4 + 1; i++) tsp_vns_3kick(path, cost, 0, tsp_inst.nnodes);
    tsp_2opt(path, cost, tsp_find_2opt_best_swap);

    tsp_free_thread(t_index);

}

/**
 * @brief (MULTITHREAD) Performs multiple kicks in parallel when we reach local optima
 * 
 * @param path The path to improve
 * @param cost The cost associated to the path
*/
void tsp_multi_kick(int* path, double* cost) {

    tsp_mt_parameters params[N_THREADS];
    int* multi_kick_paths[N_THREADS];
    double multi_kick_cost[N_THREADS];

    for (int i = 0; i < N_THREADS; i++) multi_kick_paths[i] = (int*)calloc(tsp_inst.nnodes, sizeof(int));

    while (tsp_time_elapsed() < tsp_time_limit) {

        for (int i = 0; i < N_THREADS; i++) {

            int thread = tsp_wait_for_thread();

            for (int j = 0; j < tsp_inst.nnodes; j++) multi_kick_paths[thread][j] = path[j];
            multi_kick_cost[thread] = *cost;

            params[thread].t_index = thread;
            params[thread].path = multi_kick_paths[thread];
            params[thread].cost = &multi_kick_cost[thread];

            pthread_create(&tsp_threads[thread], NULL, tsp_vns_kicknsolve, (void*)&params[thread]);

        }

        tsp_wait_all_threads();

        int min_thread = 0;
        double min_cost = multi_kick_cost[min_thread];

        for (int i = 1; i < N_THREADS; i++)
            if (multi_kick_cost[i] < min_cost) {
                min_thread = i;
                min_cost = multi_kick_cost[min_thread];
            }

        for (int i = 0; i < tsp_inst.nnodes; i++) path[i] = multi_kick_paths[min_thread][i];
        *cost = min_cost;

        #if TSP_VERBOSE >= 100
        tsp_check_integrity(path, *cost, "algorithms.c: tsp_solve_f2opt - 2");
        #endif

        tsp_check_best_sol(path, *cost, tsp_time_elapsed());

        #if TSP_VERBOSE >= 50
        tsp_save_intermediate_cost(0, *cost);
        tsp_save_intermediate_cost(1, tsp_inst.best_cost);
        #endif

    }

    for (int i = 0; i < N_THREADS; i++) if (multi_kick_paths[i] != NULL) { free(multi_kick_paths[i]); multi_kick_paths[i] = NULL; }

}

int tsp_solve_vns() {

    int* path = (int*)calloc(tsp_inst.nnodes, sizeof(int));
    double cost = tsp_greedy_path_from_node(path, rand() % tsp_inst.nnodes);

    tsp_multi_kick(path, &cost);

    if (path != NULL) { free(path); path = NULL; }

    return -1;

}


// FVNS
/**
 * @brief Calculate the block a node belongs to based on it's depth
 * 
 * @param depth The depth used to calculate the block index
 * @param node The node whose block index is needed
 * 
 * @return The block index of the node, based on the depth
*/

int tsp_calculate_block(int depth, int node) {

    return floor((depth % 2 == 1 ? tsp_inst.coords[node].y : tsp_inst.coords[node].x) * pow(2, (int)ceil(((double)depth + 1)/ 2)) / TSP_GRID_SIZE);

}

/**
 * @brief Orders the node in the path based on a grid scheme
 * 
 * @param path The path to order
*/
void tsp_partition_path(int* path) {

    tsp_entry* list = (tsp_entry*)calloc(tsp_inst.nnodes, sizeof(tsp_entry));

    int mapping[] = {   // built for max depth 8, works for any lower number
          0,   2,   8,  10,    32,  34,  40,  42,    128, 130, 136, 138,   160, 162, 168, 170,
          1,   3,   9,  11,    33,  35,  41,  43,    129, 131, 137, 139,   161, 163, 169, 171,
          4,   6,  12,  14,    36,  38,  44,  46,    132, 134, 140, 142,   164, 166, 172, 174,
          5,   7,  13,  15,    37,  39,  45,  47,    133, 135, 141, 143,   165, 167, 173, 175,
           
         16,  18,  24,  26,    48,  50,  56,  58,    144, 146, 152, 154,   176, 178, 184, 186,
         17,  19,  25,  27,    49,  51,  57,  59,    145, 147, 153, 155,   177, 179, 185, 187,
         20,  22,  28,  30,    52,  54,  60,  62,    148, 150, 156, 158,   180, 182, 188, 190,
         21,  23,  29,  31,    53,  55,  61,  63,    149, 151, 157, 159,   181, 183, 189, 191,

         64,  66,  72,  74,    96,  98, 104, 106,    192, 194, 200, 202,    224, 226, 232, 234,
         65,  67,  73,  75,    97,  99, 105, 107,    193, 195, 201, 203,    225, 227, 233, 235,
         68,  70,  76,  78,   100, 102, 108, 110,    196, 198, 204, 206,    228, 230, 236, 238,
         69,  71,  77,  79,   101, 103, 109, 111,    197, 199, 205, 207,    229, 231, 237, 239,

         80,  82,  88,  90,   112, 114, 120, 122,    208, 210, 216, 218,    240, 242, 248, 250,
         81,  83,  89,  91,   113, 115, 121, 123,    209, 211, 217, 219,    241, 243, 249, 251,
         84,  86,  92,  94,   116, 118, 124, 126,    212, 214, 220, 222,    244, 246, 252, 254,
         85,  87,  93,  95,   117, 119, 125, 127,    213, 215, 221, 223,    245, 247, 253, 255
        };

    for (int i = 0; i < tsp_inst.nnodes; i++)
        list[i] = (tsp_entry){path[i], mapping[tsp_calculate_block(7, path[i]) * 16 + tsp_calculate_block(6, path[i])] + ((double)rand())/RAND_MAX};

    qsort((void*)list, (size_t)tsp_inst.nnodes, sizeof(tsp_entry), compare_tsp_entries);

    for (int i = 0; i < tsp_inst.nnodes; i++) path[i] = list[i].key;

    free(list);

}

/**
 * @brief Finds the best swap that can be applied in the path between a start node and an end node
 * 
 * @param path The path where to operate the swap
 * @param start The starting index
 * @param end The ending index (excluded)
 * 
 * @return -1 if no improving swap has been found, 1 otherwise
*/
int tsp_f2opt_swap(int* path, int start, int end) {

    double best_improvement = -INFINITY;
    int best_start = -1, best_end = -1;

    for (int i = start; i < end - 2; i++) {        
        for (int j = i + 2; j < end; j++) {
            if (i == start && j+1 == end) continue;
            int k = (j+1 == end) ? start : j+1;  //allow for the loop over the edge

            double improvement =    (tsp_inst.costs[path[i] * tsp_inst.nnodes + path[i+1]] + tsp_inst.costs[path[j]   * tsp_inst.nnodes + path[k]]) -
                                    (tsp_inst.costs[path[i] * tsp_inst.nnodes + path[j]]   + tsp_inst.costs[path[i+1] * tsp_inst.nnodes + path[k]]);

            if (improvement > TSP_EPSILON && improvement > best_improvement + TSP_EPSILON) {
                best_improvement = improvement;
                best_start = i; best_end = j;
            }

        }
    }

    if (best_improvement < -TSP_EPSILON) return -1; //no swap found

    tsp_reverse(path, best_start+1, best_end);  //reverse the part in the middle of the swap

    return 1;    //found swap

}

/**
 * @brief Finds the splitting point based on a depth, between a start and an end index
 * 
 * @param path The path where to find the splitting index
 * @param depth The depth at which to search the splitting index
 * @param start The starting index
 * @param end The ending index
 * 
 * @return The splitting index referred to the specified depth
*/
int tsp_f2opt_split(const int* path, int depth, int start, int end) {

    int start_block = tsp_calculate_block(depth, path[start]);

    for (int i = start + 1; i < end; i++)
        if (tsp_calculate_block(depth, path[i]) != start_block) return i;

    return start;

}

/**
 * @brief (THREAD SPECIFIC) Applies swaps till no better swap can be found
 * 
 * @param params tsp_mt_parameters type structure using: ->t_index, ->s_node, ->e_node
*/
void* tsp_f2opt_block(void* params) {

    int* path = ((tsp_mt_parameters*)params) -> path;
    int start = ((tsp_mt_parameters*)params) -> s_node;
    int end = ((tsp_mt_parameters*)params) -> e_node;

    while (tsp_f2opt_swap(path, start, end) > 0 && tsp_time_elapsed() < tsp_time_limit);

    if (end - start > 10) {
        double _ = 0;
        for (int j = 0; j < 5; j++)
            tsp_vns_3kick(path, &_, start, end);

        while (tsp_f2opt_swap(path, start, end) > 0 && tsp_time_elapsed() < tsp_time_limit);
    }

    tsp_free_thread(((tsp_mt_parameters*)params)->t_index);

}

int tsp_solve_fvns() {

    int* path = (int*)calloc(tsp_inst.nnodes, sizeof(int));
    for (int i = 0; i < tsp_inst.nnodes; i++) path[i] = i;

    tsp_partition_path(path);
    
    int* partitions[TSP_F2OPT_MAX_DEPTH];
    for (int i = 0; i < TSP_F2OPT_MAX_DEPTH; i++)
        partitions[i] = (int*)calloc(pow(2, i+1)+1, sizeof(int));

    partitions[0][0] = 0;
    partitions[0][1] = tsp_f2opt_split(path, 0, 0, tsp_inst.nnodes);
    partitions[0][2] = tsp_inst.nnodes;

    for (int i = 1; i < TSP_F2OPT_MAX_DEPTH; i++)
        for (int j = 0; j < pow(2, i+1)+1; j++)
            partitions[i][j] = (j%2 == 0) ? partitions[i-1][j/2] : tsp_f2opt_split(path, i, partitions[i-1][(j-1)/2], partitions[i-1][(j+1)/2]);
    

    tsp_mt_parameters params[N_THREADS];

    for (int i = TSP_F2OPT_MAX_DEPTH-1; i >= 0; i--) {
        for (int j = 0; j < pow(2, i+1); j++) {

            int thread = tsp_wait_for_thread();

            params[thread].t_index = thread;
            params[thread].path = path;
            params[thread].s_node = partitions[i][j];
            params[thread].e_node = partitions[i][j+1];
            
            pthread_create(&tsp_threads[thread], NULL, tsp_f2opt_block, (void*)&params[thread]);

        }
        tsp_wait_all_threads();
    }

    for (int i = 0; i < TSP_F2OPT_MAX_DEPTH; i++) if (partitions[i] != NULL) { free(partitions[i]); partitions[i] = NULL; }

    double cost = tsp_compute_path_cost(path);
    tsp_2opt(path, &cost, tsp_find_2opt_best_swap);

    #if TSP_VERBOSE >= 100
    tsp_check_integrity(path, cost, "algorithms.c: tsp_solve_f2opt - 1");
    #endif

    tsp_check_best_sol(path, cost, tsp_time_elapsed());

    #if TSP_VERBOSE >= 10
    printf("---------------------- Finished f2opt, starting vns.\n");
    #endif

    tsp_multi_kick(path, &cost);

    if (path != NULL) { free(path); path = NULL; }

    return -1;

}


// CPLEX

/**
 * @brief Converts coordinates to cplex x_pos
 * 
 * @param i Row coordinate
 * @param j Col coordinate
 * 
 * @return The index used by cplex to locate the edge (i, j)
*/
int tsp_cplex_coords_to_xpos(const int i, const int j) {

	if ( i == j ) { printf("ERROR: i == j in xpos"); exit(1); }
	if ( i > j ) return tsp_cplex_coords_to_xpos(j,i);

	return i * tsp_inst.nnodes + j - (( i + 1 ) * ( i + 2 )) / 2;

}

/**
 * @brief Builds the cplex model from the tsp_inst initialized
 * 
 * @param env cplex pointer to the cplex environment
 * @param lp cplex pointer to the cplex linear problem
*/
void tsp_cplex_build_model(CPXENVptr env, CPXLPptr lp) {

	int izero = 0;
	char binary = 'B'; 
	
	char** cname = (char**)calloc(1, sizeof(char *));	// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));

	// add binary var.s x(i,j) for i < j  
	for ( int i = 0; i < tsp_inst.nnodes; i++ ) {

		for ( int j = i+1; j < tsp_inst.nnodes; j++ ) {

			sprintf(cname[0], "x(%d,%d)", i+1,j+1);  // x(1,2), x(1,3) ....
			double obj = tsp_inst.costs[i * tsp_inst.nnodes + j]; // cost = distance   
			double lb = 0.0;
			double ub = 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) { printf("ERROR: wrong CPXnewcols on x var.s"); exit(1); }
    		if ( CPXgetnumcols(env,lp)-1 != tsp_cplex_coords_to_xpos(i,j) ) { printf("ERROR: wrong position for x var.s"); exit(1); }

		}

	} 

	// add degree constr.s 
	int* index = (int*)malloc(tsp_inst.nnodes * sizeof(int));
	double* value = (double*)malloc(tsp_inst.nnodes * sizeof(double));  
	
	// add the degree constraints
	for ( int h = 0; h < tsp_inst.nnodes; h++ ) { // degree constraints

		double rhs = 2.0;
		char sense = 'E';                     // 'E' for equality constraint 
		sprintf(cname[0], "degree(%d)", h+1); 
		int nnz = 0;
		for ( int i = 0; i < tsp_inst.nnodes; i++ ) {
			if ( i == h ) continue;
			index[nnz] = tsp_cplex_coords_to_xpos(i,h);
			value[nnz] = 1.0;
			nnz++;
		}
		
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) ) { printf("ERROR: wrong CPXaddrows [degree]"); exit(1); }

	} 

    free(value);
    free(index);	

	#if TSP_VERBOSE >= 100
    CPXwriteprob(env, lp, "model.lp", NULL);
    #endif

	free(cname[0]);
	free(cname);

}

void tsp_cplex_set_params(CPXENVptr env, CPXLPptr lp) {

    //TODO: do stuff

}

void tsp_cplex_save_solution(CPXENVptr env, CPXLPptr lp) {

    //TODO: Now it just prints it in the cplex format, must convert it to the tsp_inst.best_solution list...

    int ncols = CPXgetnumcols(env, lp);
	double *xstar = (double *) calloc(ncols, sizeof(double));
	
    //CPXgetx get the solution x* of our instance
	if ( CPXgetx(env, lp, xstar, 0, ncols-1) ) { printf("ERROR: CPXgetx() error"); exit(1); }	

	for ( int i = 0; i < tsp_inst.nnodes; i++ ) {

		for ( int j = i+1; j < tsp_inst.nnodes; j++ ) {

			if ( xstar[tsp_cplex_coords_to_xpos(i,j)] > 0.5 ) printf("  ... x(%3d,%3d) = 1\n", i+1,j+1);

		}

	}

	free(xstar);

}

int tsp_solve_cplex() {

    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

    tsp_cplex_build_model(env, lp);

    tsp_cplex_set_params(env, lp);

    if (CPXmipopt(env,lp)) { printf("ERROR: CPXmipopt() error"); exit(1); }

    tsp_cplex_save_solution(env, lp);

	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env); 

    //TODO: check the time limit with cplex error code
    return 1;

}