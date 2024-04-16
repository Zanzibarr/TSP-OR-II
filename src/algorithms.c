#include "../include/algorithms.h"


#pragma region HEURISTIC_2OPT
/**
 * @brief Applies a swap policy till the solution cannot be improved further
 * 
 * @param path Path considered (will be changed if it finds a swap)
 * @param cost Cost of the current path (will be changed if it finds a swap)
 * @param swap_function The swap function to use
*/
void tsp_2opt(int* path, double* cost, int (*swap_function)(int*, double*)) {

    while ((*swap_function)(path, cost) > 0 && tsp_time_elapsed() < tsp_time_limit) { //repeat until I can't find new swaps that improve the cost of my solution

        if (tsp_verbose >= 50) {
            tsp_save_intermediate_cost(0, *cost);
            tsp_save_intermediate_cost(1, tsp_inst.best_cost);
        }
    
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
#pragma endregion


#pragma region HEURISTIC_GREEDY
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

    if (tsp_verbose >= 100) tsp_check_integrity(path, cost, "algorithms.c: tsp_greedy_from_node - 1");

    if (((tsp_mt_parameters*)params)->swap_function != NULL) {   //if I want to use the 2opt

        tsp_2opt(path, &cost, ((tsp_mt_parameters*)params)->swap_function);  //fix solution using 2opt

        if (tsp_verbose >= 100) tsp_check_integrity(path, cost, "algorithms.c: tsp_greedy_from_node - 2");

    }

    tsp_check_best_sol(path, cost, tsp_time_elapsed()); //if this solution is better than the best seen so far update it

    if (tsp_verbose >= 50) tsp_save_intermediate_cost(0, cost);

    if(path != NULL) { free(path); path = NULL; }

    tsp_free_thread(((tsp_mt_parameters*)params)->t_index);

    return NULL;

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

        if (tsp_verbose >= 100) tsp_check_integrity(path, cost, "algorithms.c: tsp_solve_greedy_st - 1");

        if (swap_function != NULL) {   //if I want to use the 2opt

            tsp_2opt(path, &cost, swap_function);  //fix solution using 2opt

            if (tsp_verbose >= 100) tsp_check_integrity(path, cost, "algorithms.c: tsp_solve_greedy_st - 2");

        }

        double elapsed_time = tsp_time_elapsed();
        
        tsp_check_best_sol(path, cost, elapsed_time); //if this solution is better than the best seen so far update it

        if (tsp_verbose >= 50) tsp_save_intermediate_cost(0, cost);

        if (elapsed_time > tsp_time_limit) { if(path != NULL) { free(path); path = NULL; } return -1; } //if I exceeded the time limit

    }

    if (path != NULL) { free(path); path = NULL; }

    return 0;

}
#pragma endregion


#pragma region METAHEURISTIC_TABU
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

        if (tsp_verbose >= 50) tsp_save_intermediate_cost(t_index, *cost);

        if (tsp_verbose >= 100) tsp_check_integrity(path, *cost, "algorithms.c: tsp_tabu_from_node");

        tsp_check_best_sol(path, *cost, time);

    }

    tsp_free_thread(t_index);

    return NULL;

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
#pragma endregion


#pragma region METAHEURISTIC_VNS
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

    return NULL;

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

        if (tsp_verbose >= 100) tsp_check_integrity(path, *cost, "algorithms.c: tsp_solve_f2opt - 2");

        tsp_check_best_sol(path, *cost, tsp_time_elapsed());

        if (tsp_verbose >= 50) {
            tsp_save_intermediate_cost(0, *cost);
            tsp_save_intermediate_cost(1, tsp_inst.best_cost);
        }

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
#pragma endregion


#pragma region METAHEURISTIC_FVNS
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

    qsort((void*)list, (size_t)tsp_inst.nnodes, sizeof(tsp_entry), tsp_compare_entries);

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

    return NULL;

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

    if (tsp_verbose >= 100) tsp_check_integrity(path, cost, "algorithms.c: tsp_solve_f2opt - 1");

    tsp_check_best_sol(path, cost, tsp_time_elapsed());

    if (tsp_verbose >= 10) tsp_print_info("Finished f2opt, starting vns.\n");

    tsp_multi_kick(path, &cost);

    if (path != NULL) { free(path); path = NULL; }

    return -1;

}
#pragma endregion


#pragma region CPLEX

int tsp_cplex_solve_model(CPXENVptr env, CPXLPptr lp, double* xstar, int* ncomp, int* comp, int* succ, double* cost) {

    CPXsetdblparam(env, CPXPARAM_TimeLimit, tsp_time_limit - tsp_time_elapsed());
    /*int mipopt_output = CPXmipopt(env, lp);
    if ( mipopt_output==CPXERR_DETTILIM_STRONGBRANCH || mipopt_output==CPXERR_PRESLV_DETTIME_LIM ||
        mipopt_output==CPXERR_PRESLV_TIME_LIM || mipopt_output==CPXERR_TILIM_CONDITION_NO ||
        mipopt_output==CPXERR_TILIM_STRONGBRANCH) return -1;
    else if ( mipopt_output ) { printf("CPXmipopt error.\n"); exit(1); }*/
    if ( CPXmipopt(env, lp) ) { tsp_print_error("CPXmipopt error.\n"); }
    tsp_cplex_save_solution(env, lp, xstar, cost);
    tsp_cplex_build_solution(xstar, ncomp, comp, succ); 
    //CPXgetbestobjval(env, lp, cost);
    return 0;

}

int tsp_cplex_benders_loop(CPXENVptr env, CPXLPptr lp, double* xstar, int* ncomp, int* comp, int* succ, double* cost, char patching) {

    int iter = 1;
    char infeasible_solution = 0;
    while (tsp_time_elapsed() < tsp_time_limit) {

        if (tsp_cplex_solve_model(env, lp, xstar, ncomp, comp, succ, cost)) {  // method exceeded time limit while solving model
            if (infeasible_solution) return 0;  // we already have an infeasible solution
            else return -1; // no solution found yet
        }
        else infeasible_solution = 1;

        if (tsp_verbose >= 10) {
            tsp_print_info("Iteration number %d; %d connected components; %f elapsed time; %f current incumbent\n", iter++, *ncomp, tsp_time_elapsed(), *cost);
        }

        if (*ncomp==1) break;   // feasible solution found

        tsp_cplex_add_sec(env, lp, ncomp, comp, succ);

        if (patching) tsp_cplex_patch_comp(xstar, ncomp, comp, succ, cost);

    }

    tsp_cplex_convert_solution(ncomp, succ, cost);

    if (patching) tsp_2opt(tsp_inst.best_solution, &tsp_inst.best_cost, tsp_find_2opt_best_swap);

    //if (tsp_verbose>=100 && tsp_cplex_sol.ncomp==1) tsp_check_integrity();

    if (*ncomp==1) return 1;
    else return 0;   // infeasible solution found

}

void tsp_cplex_patch_comp(double* xstar, int* ncomp, int* comp, int* succ, double* cost) {

    int starts[*ncomp];
    for (int i=0; i<*ncomp; i++) starts[i]=-1;

    // glue components together until exceeding time limit or we have only one component left
    while (/*tsp_time_elapsed() < tsp_time_limit &&*/ (*ncomp)!=1) {

        // edge to be removed expressed by first node in succ order
        int best_k1 = 0, best_k2 = 0, best_edge_k1 = 0, best_edge_k2 = 0;
        double best_delta = -INFINITY, delta_N, delta_R;
        // move stores how to glue components: given edges (i,succ[i]) and (j, succ[j]), replace them with:
        // - move='R' (Reverse): (i,j) and (succ[j],succ[i]) -> reverse k2 afterwards
        // - move='N' (Not reverse): (i,succ[j]) and (j,succ[i]) -> do not reverse k2 afterwards
        char best_move;

        // look for best patching move
        for (int k1=1; k1<=(*ncomp); k1++) {

            int start_k1 = starts[k1-1];
            if (start_k1==-1) {
                for (start_k1=0; comp[start_k1]!=k1; start_k1++);
                starts[k1-1] = start_k1;
            }
            for (int k2=k1+1; k2<=(*ncomp); k2++) {

                int start_k2 = starts[k2-1];
                if (start_k2==-1) {
                    for (start_k2=0; comp[start_k2]!=k2; start_k2++);
                    starts[k2-1] = start_k2;
                }
                int current_k1 = start_k1, current_k2 = start_k2,
                    succ_k1 = succ[current_k1], succ_k2 = succ[current_k2];
                do {
                    
                    do {

                        delta_N =   (tsp_get_edge_cost(current_k1,succ_k1) + tsp_get_edge_cost(current_k2, succ_k2)) -
                                    (tsp_get_edge_cost(current_k1,succ_k2) + tsp_get_edge_cost(current_k2, succ_k1));
                        delta_R =   (tsp_get_edge_cost(current_k1,succ_k1) + tsp_get_edge_cost(current_k2, succ_k2)) -
                                    (tsp_get_edge_cost(current_k1,current_k2) + tsp_get_edge_cost(succ_k2, succ_k1));
                                        
                        if (delta_N>=delta_R && delta_N>best_delta) {
                            best_k1 = k1; best_k2 = k2; best_edge_k1 = current_k1; best_edge_k2 = current_k2;
                            best_delta = delta_N;
                            best_move = 'N';
                        }
                        if (delta_R>delta_N && delta_R>best_delta) {
                            best_k1 = k1; best_k2 = k2; best_edge_k1 = current_k1; best_edge_k2 = succ_k2;
                            best_delta = delta_R;
                            best_move = 'R';
                        }

                        current_k2 = succ_k2; succ_k2 = succ[current_k2];

                    } while (current_k2!=start_k2);
                    current_k1 = succ_k1; succ_k1 = succ[current_k1];

                } while (current_k1!=start_k1);

            }

        }

        /*if (tsp_verbose >= 10)
            printf("Best swap: %d (comp %d), %d (comp %d), improv. %f, move %c\n",
                best_edge_k1, best_k1, best_edge_k2, best_k2, best_delta, best_move);*/
        // glue components
        // if needed, reverse orientation of k2
        if (best_move=='R') {
            int start;
            for (start=0; comp[start]!=best_k2; start++);
            int next = start, carry = succ[start], curr;
            do {
                curr = next;
                next = carry;
                if (next!=start) carry = succ[next];
                succ[next] = curr;
            } while (next!=start);
        }
        // change edges
        int tmp = succ[best_edge_k1];
        succ[best_edge_k1] = succ[best_edge_k2];
        succ[best_edge_k2] = tmp;
        // update comp
        for (int i=0; i<tsp_inst.nnodes; i++) {
            if (comp[i]==best_k2) comp[i]=best_k1;
            else if (comp[i]>best_k2) comp[i]--;
        }
        // update ncomp
        (*ncomp)--;
        // update starts
        for (int i=best_k2-1; i<(*ncomp)-1; i++) starts[i] = starts[i+1];
        starts[(*ncomp)-1]=-1;
        // update cost
        (*cost) -= best_delta;

    }

}

int tsp_solve_cplex() {

    int error;
    CPXENVptr tsp_cplex_env = CPXopenCPLEX(&error);
	CPXLPptr tsp_cplex_lp = CPXcreateprob(tsp_cplex_env, &error, "TSP");

    // set cplex log file
    CPXsetdblparam(tsp_cplex_env, CPX_PARAM_SCRIND, CPX_OFF);
    char cplex_log_file[100];
    sprintf(cplex_log_file, "%s/%d-%d-%s.log", TSP_CPLEX_LOG_FOLDER, tsp_seed, tsp_inst.nnodes, tsp_alg_type);
    remove(cplex_log_file);
    if ( CPXsetlogfilename(tsp_cplex_env, cplex_log_file, "w") )
        { tsp_print_error("CPXsetlogfilename error\n"); }

    // build cplex model
    tsp_cplex_build_model(tsp_cplex_env, tsp_cplex_lp);

    // create lp file from cplex model
    char cplex_lp_file[100];
    sprintf(cplex_lp_file, "%s/%d-%d-%s.lp", TSP_CPLEX_LP_FOLDER, tsp_seed, tsp_inst.nnodes, tsp_alg_type);
    if ( CPXwriteprob(tsp_cplex_env, tsp_cplex_lp, cplex_lp_file, NULL) )
        { tsp_print_error("CPXwriteprob error\n"); }

    // pick algorithm
    double* xstar = (double*) calloc(CPXgetnumcols(tsp_cplex_env, tsp_cplex_lp), sizeof(double));
    int ncomp = 0;
    int* comp = (int*) calloc(tsp_inst.nnodes, sizeof(int));
    int* succ = (int*) calloc(tsp_inst.nnodes, sizeof(int));
    double cplex_cost;
    int output;
    switch (tsp_find_alg(tsp_alg_type)) {
        case 6:
            output = tsp_cplex_solve_model(tsp_cplex_env, tsp_cplex_lp, xstar, &ncomp, comp, succ, &cplex_cost);
            tsp_cplex_convert_solution(&ncomp, comp, &cplex_cost);
            break;
        case 7:
            output = tsp_cplex_benders_loop(tsp_cplex_env, tsp_cplex_lp, xstar, &ncomp, comp, succ, &cplex_cost, 0);
            break;
        case 8:
            output = tsp_cplex_benders_loop(tsp_cplex_env, tsp_cplex_lp, xstar, &ncomp, comp, succ, &cplex_cost, 1);
            break;
    }

    // store solution in special data structure if it has multiple tours
    tsp_multi_sol.ncomp = ncomp;
    if (tsp_multi_sol.ncomp>1) {

        tsp_multi_sol.comp = (int*) calloc(tsp_inst.nnodes, sizeof(int));
        for (int i=0; i<tsp_inst.nnodes; i++) tsp_multi_sol.comp[i]=comp[i];
        tsp_multi_sol.succ = (int*) calloc(tsp_inst.nnodes, sizeof(int));
        for (int i=0; i<tsp_inst.nnodes; i++) tsp_multi_sol.succ[i]=succ[i];

    }

    // free memory
	CPXfreeprob(tsp_cplex_env, &tsp_cplex_lp);
	CPXcloseCPLEX(&tsp_cplex_env);
    free(xstar); free(comp); free(succ);

    // remove "clone<x>.log" files generated by cplex
    int file_number = 1;
    char clone_file[50];
    sprintf(clone_file, "clone0.log");
    remove(clone_file);
    do { sprintf(clone_file, "clone%d.log", file_number++); } while (!remove(clone_file)); 

    return output;

}
#pragma endregion