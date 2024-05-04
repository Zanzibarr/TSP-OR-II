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

            double improvement =    (tsp_get_edge_cost(path[i], path[i+1]) + tsp_get_edge_cost(path[j], path[k])) -
                                    (tsp_get_edge_cost(path[i], path[j])   + tsp_get_edge_cost(path[i+1], path[k]));

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

            double improvement =    (tsp_get_edge_cost(path[i], path[i+1]) + tsp_get_edge_cost(path[j], path[k])) -
                                    (tsp_get_edge_cost(path[i], path[j])   + tsp_get_edge_cost(path[i+1], path[k]));

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
        cost += tsp_get_edge_cost(frontier, next);
        frontier = next;

    }

    if (visited != NULL) { free(visited); visited = NULL; }

    cost += tsp_get_edge_cost(frontier, start_node);    //add the cost of the last edge

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

            double improvement =    (tsp_get_edge_cost(path[i], path[i+1]) + tsp_get_edge_cost(path[j], path[k])) -
                                    (tsp_get_edge_cost(path[i], path[j])   + tsp_get_edge_cost(path[i+1], path[k]));
            
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
    i = (rand() % (end - start - 1)) + start;   //FIXME: This is not thread safe
    j = i;
    while (j >= i - 1 && j <= i + 1) j = (rand() % (end - start - 1)) + start;   //FIXME: This is not thread safe
    k = i;
    while (k >= i - 1 && k <= i + 1 || k >= j - 1 && k <= j + 1) k = (rand() % (end - start - 1)) + start;   //FIXME: This is not thread safe

    // sort them
    if (i > k) { int c = i; i = k; k = c; }
    if (i > j) { int c = i; i = j; j = c; }
    if (j > k) { int c = j; j = k; k = c; }

    // update cost
    *cost = *cost -
        (   //removing the costs of the old edges
            (tsp_get_edge_cost(path[i], path[i+1])) +
            (tsp_get_edge_cost(path[j], path[j+1])) +
            (tsp_get_edge_cost(path[k], path[(k+1==tsp_inst.nnodes)?0:k+1]))
        ) +
        (   //adding the costs of the new edges
            (tsp_get_edge_cost(path[i], path[k])) +
            (tsp_get_edge_cost(path[j+1], path[i+1])) +
            (tsp_get_edge_cost(path[j], path[(k+1==tsp_inst.nnodes)?0:k+1]))
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


#pragma region HEURISTIC_F2OPT

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

            double improvement =    (tsp_get_edge_cost(path[i], path[i+1]) + tsp_get_edge_cost(path[j], path[k])) -
                                    (tsp_get_edge_cost(path[i], path[j])   + tsp_get_edge_cost(path[i+1], path[k]));

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

/**
 * @brief Apply f2opt
 * 
 * @param path An array to store the solution into
 * @param cost A variable to store the cost of the solution found
*/
double tsp_f2opt(int* path) {

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

    return cost;

}

#pragma endregion


#pragma region METAHEURISTIC_FVNS

int tsp_solve_fvns() {

    int* path = (int*)calloc(tsp_inst.nnodes, sizeof(int));

    double cost = tsp_f2opt(path);

    tsp_check_best_sol(path, cost, tsp_time_elapsed());

    if (tsp_verbose >= 10) tsp_print_info("Finished f2opt, starting vns.\n");

    tsp_multi_kick(path, &cost);

    if (path != NULL) { free(path); path = NULL; }

    return -1;

}

#pragma endregion


#pragma region CPLEX


/**
 * @brief solve the current model through cplex and computes the connected components of solution found
 * 
 * @return int 0 if model was solved before timelimit, -1 if time limit exceeded and no solution found, -2 if time limit exceeded and feasible solution to model found, -3 if time limit exceeded and infeasible solution to model found
 */
int tsp_cplex_solve_model(CPXENVptr env, CPXLPptr lp, double* xstar, int* ncomp, int* comp, int* succ, double* cost) {

    CPXsetdblparam(env, CPXPARAM_TimeLimit, tsp_time_limit - tsp_time_elapsed());
    if ( CPXmipopt(env, lp) ) tsp_raise_error("CPXmipopt error.\n");

    double tmp;
    if ( CPXgetbestobjval(env, lp, &tmp) ) tsp_raise_error("No MIP solution available in cplex model.");
    int output;     //FIXME: Guarda il metodo tsp_cplex_solve() e sistema i return code
    switch ( CPXgetstat(env, lp) ) {
        case CPXMIP_TIME_LIM_FEAS:
            output = -2;
            break;
        case CPXMIP_TIME_LIM_INFEAS:
            output = -1;
            break;
        default:
            output = 0;
            break;
    }

    tsp_cplex_save_solution(env, lp, xstar, cost);          // save the solution into 'xstar' and saves the cost into 'cost'
    tsp_cplex_build_solution(xstar, ncomp, comp, succ);     // split 'xstar' into 'comp' and 'succ', storing number of connected components into 'ncomp'
    if (output==-2 && *ncomp>1) output = -3;

    return output;

}

/**
 * @brief Applies the bender loop to add SECs
 * 
 * @return int 1 if feasible solution was found before timelimit, 0 if timelimit exceeded but infeasible solution found, -1 otherwise
*/
int tsp_cplex_benders_loop(CPXENVptr env, CPXLPptr lp, double* xstar, int* ncomp, int* comp, int* succ, double* cost, char patching) {

    int iter = 1;
    int output;

    while (tsp_time_elapsed() < tsp_time_limit) {

        output = tsp_cplex_solve_model(env, lp, xstar, ncomp, comp, succ, cost);

        iter++;
        if (tsp_verbose >= 10) tsp_print_info("Iteration number: %4d - connected components: %4d - current incumbent: %10.4f\n", iter, *ncomp, *cost);

        if (*ncomp==1) break;   // feasible solution found
        else if (output==-1) return output;

        tsp_cplex_add_sec(env, lp, ncomp, comp, succ);

        if (patching) {
            tsp_cplex_patch_comp(xstar, ncomp, comp, succ, cost);
            if (output==-3) output=-2;
        }

        if (output==-2 || output==-3) break;

    }

    tsp_cplex_convert_solution(ncomp, succ, cost);

    if (patching) tsp_2opt(tsp_inst.best_solution, &tsp_inst.best_cost, tsp_find_2opt_best_swap);

    //if (tsp_verbose>=100 && tsp_cplex_sol.ncomp==1) tsp_check_integrity();

    if (!output && *ncomp==1) return 0;
    else return output;

}

int tsp_solve_cplex_base(CPXENVptr env, CPXLPptr lp, double* xstar, int* ncomp, int* comp, int* succ, double* cost) {

    int output = tsp_cplex_solve_model(env, lp, xstar, ncomp, comp, succ, cost);
    tsp_cplex_convert_solution(ncomp, comp, cost);
    if ( output==-1 ) return output;
    if ( *ncomp>1 ) return -4;
    return 0;

}

int tsp_solve_cplex() {

    int error;
    CPXENVptr tsp_cplex_env = CPXopenCPLEX(&error);
	CPXLPptr tsp_cplex_lp = CPXcreateprob(tsp_cplex_env, &error, "TSP");

    // set cplex log file
    CPXsetdblparam(tsp_cplex_env, CPX_PARAM_SCRIND, CPX_OFF);
    char cplex_log_file[100];
    sprintf(cplex_log_file, "%s/%lu-%d-%s.log", TSP_CPLEX_LOG_FOLDER, tsp_seed, tsp_inst.nnodes, tsp_alg_type);
    remove(cplex_log_file);
    if ( CPXsetlogfilename(tsp_cplex_env, cplex_log_file, "w") ) tsp_raise_error("CPXsetlogfilename error.\n");

    // build cplex model
    tsp_cplex_build_model(tsp_cplex_env, tsp_cplex_lp);

    // create lp file from cplex model
    char cplex_lp_file[100];
    sprintf(cplex_lp_file, "%s/%lu-%d-%s.lp", TSP_CPLEX_LP_FOLDER, tsp_seed, tsp_inst.nnodes, tsp_alg_type);
    if ( CPXwriteprob(tsp_cplex_env, tsp_cplex_lp, cplex_lp_file, NULL) ) tsp_raise_error("CPXwriteprob error\n");

    // pick algorithm
    double* xstar = (double*) calloc(CPXgetnumcols(tsp_cplex_env, tsp_cplex_lp), sizeof(double));
    int ncomp = 0;
    int* comp = (int*) calloc(tsp_inst.nnodes, sizeof(int));
    int* succ = (int*) calloc(tsp_inst.nnodes, sizeof(int));
    double cplex_cost;
    // values: 0 feasible within timelimit, -1 no solution found because of time limit,
    // -2 feasible outside time limit, -3 infeasible outside time limit
    // -4 infeasible within timelimit (used only by cplex-base)
    int output;
    switch (tsp_find_alg(tsp_alg_type)) {
        case 6:
            output = tsp_solve_cplex_base(tsp_cplex_env, tsp_cplex_lp, xstar, &ncomp, comp, succ, &cplex_cost);
            break;
        case 7:
            output = tsp_cplex_benders_loop(tsp_cplex_env, tsp_cplex_lp, xstar, &ncomp, comp, succ, &cplex_cost, 0);    //FIXME: Wrong costs
            break;
        case 8:
            output = tsp_cplex_benders_loop(tsp_cplex_env, tsp_cplex_lp, xstar, &ncomp, comp, succ, &cplex_cost, 1);    //FIXME: Wrong costs
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
    if (xstar != NULL) { free(xstar); xstar = NULL; }
    if (comp != NULL) { free(comp); comp = NULL; }
    if (succ != NULL) { free(succ); succ = NULL; }

    // remove "clone<x>.log" files generated by cplex
    int file_number = 1;
    char clone_file[50];
    sprintf(clone_file, "clone0.log");
    remove(clone_file);
    do { sprintf(clone_file, "clone%d.log", file_number++); } while (!remove(clone_file)); 

    return output;

}

#pragma endregion


#pragma region BRANCH&CUT

/**
 * @brief Solve using cplex
 * 
 * @param env cplex environment
 * @param lp cplex lp
 * @param succ successors type list containing the solution found by cplex
 * 
 * @return int 0 if the model was solved before the time limit, -2 if an intermediate solution has been found but cplex couldn't end, -3 if no solution has been found, -4 if the problem has been proven infeasible
*/
int tsp_cplex_solve(CPXENVptr env, CPXLPptr lp, int* succ) {

    // check time limit (precomputing + first heuristic might have required all the time)
    if (tsp_time_limit - tsp_time_elapsed() < TSP_EPSILON) {
        tsp_print_warn("The time limit is too short for this algorithm.\n");
        return -3;
    }

    // set the time limit
    if ( CPXsetdblparam(env, CPXPARAM_TimeLimit, tsp_time_limit - tsp_time_elapsed()) ) tsp_raise_error("CPXsetdblparam error.\n");

    // solve the model using cplex
    if ( CPXmipopt(env, lp) ) tsp_raise_error("CPXmipopt error.\n");
    
    // get the output status (time limit (intermediate solution found / not found), infeasible)
    int output;
    switch ( CPXgetstat(env, lp) ) {
        case CPXMIP_TIME_LIM_FEAS:      // exceeded time limit, found intermediate solution
            output = -2;
            break;
        case CPXMIP_TIME_LIM_INFEAS:    // exceeded time limit, no intermediate solution found
            return -3;
        case CPXMIP_INFEASIBLE:         // proven to be unfeasible
            return -4;
        default:                        // found optimal solution
            output = 0;
            break;
    }

    // space for temporary data structures 
    double* xstar = (double*) calloc(CPXgetnumcols(env, lp), sizeof(double));
    int* comp = (int*) calloc(tsp_inst.nnodes, sizeof(int));    // I don't need this, I allocate it but then don't use it...
    int ncomp = 0;

    // convert xstar to successors type list
	if ( CPXgetx(env, lp, xstar, 0, CPXgetnumcols(env, lp)-1) ) tsp_raise_error("CPXgetx() error\n");
    tsp_cplex_decompose_xstar(xstar, comp, succ, &ncomp);

    // free memory
    if (xstar != NULL) { free(xstar); xstar = NULL; }
    if (comp != NULL) { free(comp); comp = NULL; }
    
    // return status code from cplex
    return output;

}

/**
 * @brief Generic callback function for cplex
 * 
 * @param context cplex context
 * @param context_id cplex context id
 * @param user_handle user data to pass to this function
 * 
 * @return cplex error code: 0 if everything went according to plan
*/
static int CPXPUBLIC tsp_cplex_callback(CPXCALLBACKCONTEXTptr context, CPXLONG context_id, void* user_handle) {

    switch (context_id) {

        case CPX_CALLBACKCONTEXT_CANDIDATE:

            // check for connected components in the candidate solution
            return tsp_cplex_callback_candidate(context, *((int*)user_handle));

        case CPX_CALLBACKCONTEXT_RELAXATION:

            // do stuff with fractionary solution
            return tsp_cplex_callback_relaxation(context, *((int*)user_handle));

        default:
            tsp_raise_error("Callback called for wrong reason: %d.\n", context_id);

    }

    // cplex error code
	return 1;

}

int tsp_solve_cplex_bnc() {

    int error;
    CPXENVptr env = NULL;
	CPXLPptr lp = NULL;

    tsp_cplex_init(&env, &lp, &error);

    // fast heuristic to get a solution in case cplex can't find one
    int* path = (int*) malloc(tsp_inst.nnodes * sizeof(int));
    double cost = tsp_f2opt(path);

    if (tsp_verbose >= 0) tsp_print_info("--------- Finished f2opt (cost: %10.4f), starting cplex-bnc. ---------\n", cost);

    int ncols = CPXgetnumcols(env, lp);
    
    // set callback function
    CPXLONG context_id = CPX_CALLBACKCONTEXT_CANDIDATE;
    if (tsp_cplex_rel_cb) {
        if (tsp_verbose >= 0) tsp_print_info("Adding callback function for relaxation to cplex.\n");
        context_id = context_id | CPX_CALLBACKCONTEXT_RELAXATION;
    }
    if ( CPXcallbacksetfunc(env, lp, context_id, tsp_cplex_callback, (void*)&tsp_inst.nnodes) ) tsp_raise_error("CPXcallbacksetfunc() error.");

    // add a starting heuristic to cplex
    //TODO: Since this is way better than cplex's heuristics, shall I keep improving this heuristic with vns or tabu and give the new solutions to cplex?
    //  : Do I allocate ncols or nnodes here?
    if (tsp_cplex_mipstart) {

        if (tsp_verbose >= 0) tsp_print_info("Using an heuristic as mipstart for cplex.\n");

        int *indexes = (int *) malloc(ncols * sizeof(int));
        double *values = (double *) malloc(ncols * sizeof(double));
        int effortlevel = CPX_MIPSTART_NOCHECK;
        int izero = 0;
        tsp_cplex_path_to_xstar(ncols, path, indexes, values);
        if ( CPXaddmipstarts(env, lp, 1, tsp_inst.nnodes, &izero, indexes, values, &effortlevel, NULL) ) tsp_raise_error("CPXaddmipstarts() error");
        if (indexes != NULL) { free(indexes); indexes = NULL; }
        if (values != NULL) { free(values); values = NULL; }
    }

    // space for data structures
    int* succ = (int*) calloc(tsp_inst.nnodes, sizeof(int));

    // solve cplex
    int ret_code = tsp_cplex_solve(env, lp, succ);

    // handle cplex response
    switch (ret_code) {
        case 0:
            tsp_cplex_store_solution(succ);
            break;
        case -2:
            tsp_print_warn("cplex exceeded the time limit, using an intermediate integer solution improved with 2opt.\n");
            int* path = (int*)malloc(tsp_inst.nnodes * sizeof(int));
            double cost = tsp_succ_to_path(succ, path);

            while (tsp_find_2opt_best_swap(path, &cost) > 0);   //2opt to fix solution //FIXME: I'm going out of the time limit...

            tsp_check_best_sol(path, cost, tsp_time_elapsed());
            if (path != NULL) { free(path); path = NULL; }
            break;
        case -3:
            tsp_print_warn("cplex couldn't find any feasible solution within the time limit. Returning output of f2opt.\n");
            tsp_check_best_sol(path, cost, tsp_time_elapsed()); //TODO: I could simply store the cplex solution since I gave him my initial heuristic
            break;
        case -4:
            return -4;
        default:
            tsp_raise_error("Unexpected return code from cplex_solve.\n");
    }
    
    // free memory
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
    if (succ != NULL) { free(succ); succ = NULL; }

    // remove "clone<x>.log" files generated by cplex
    int file_number = 1;
    char clone_file[50];
    sprintf(clone_file, "clone0.log");
    remove(clone_file);
    do { sprintf(clone_file, "clone%d.log", file_number++); } while (!remove(clone_file)); 

    return ret_code;

}

#pragma endregion