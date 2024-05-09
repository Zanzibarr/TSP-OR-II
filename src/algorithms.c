#include "../include/algorithms.h"


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

    safe_free(visited);

    cost += tsp_get_edge_cost(frontier, start_node);    //add the cost of the last edge

    if (tsp_verbose >= 100) tsp_check_integrity(path, cost, "algorithms.c - tsp_greedy_path_from_node - 1");

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

    if (((tsp_mt_parameters*)params)->swap_function != NULL)   //if I want to use the 2opt
        tsp_2opt(path, &cost, ((tsp_mt_parameters*)params)->swap_function);  //fix solution using 2opt

    tsp_check_best_sol(path, cost, time_elapsed()); //if this solution is better than the best seen so far update it

    safe_free(path);

    tsp_free_thread(((tsp_mt_parameters*)params)->t_index);

    return NULL;

}

void tsp_solve_greedy() {

    tsp_mt_parameters params[N_THREADS];

    for (int i = 0; i < tsp_inst.nnodes; i++) {

        if (time_elapsed() > tsp_env.time_limit) break;

        int thread = tsp_wait_for_thread();

        params[thread].t_index = thread;
        params[thread].s_node = i;
        params[thread].swap_function = NULL;

        pthread_create(&tsp_threads[thread], NULL, tsp_greedy_from_node, (void*)&params[thread]);

    }

    tsp_wait_all_threads();

    if (time_elapsed() > tsp_env.time_limit) tsp_env.status = 1;
    else tsp_env.status = 0;

}

#pragma endregion


#pragma region HEURISTIC_G2OPT

/**
 * @brief Performs a 3opt kick on the path updating the cost
 * 
 * @param path The path on which to perform the kick
 * @param cost The cost to be updated after the kick
 * @param start The starting index
 * @param end The ending index
 * @param seed The seed to use for thread-safe random
 */
void tsp_vns_3kick(int* path, double* cost, const int start, const int end, unsigned int* seed) {

    if (end - start < 10) return;

    int i = -1, j = -1, k = -1;

    // pick random values that are at least 2 positions apart from each other
    i = (rand_r(seed) % (end - start - 1)) + start;
    j = i;
    while (j >= i - 1 && j <= i + 1) j = (rand_r(seed) % (end - start - 1)) + start;
    k = i;
    while (k >= i - 1 && k <= i + 1 || k >= j - 1 && k <= j + 1) k = (rand_r(seed) % (end - start - 1)) + start;

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
    safe_free(temp_path);

}

/**
 * @brief Calculate the block a node belongs to based on it's depth
 * 
 * @param depth The depth used to calculate the block index
 * @param node The node whose block index is needed
 * 
 * @return The block index of the node, based on the depth
*/
int tsp_f2opt_calculate_block(const int depth, const int node) {

    return floor((depth % 2 == 1 ? tsp_inst.coords[node].y : tsp_inst.coords[node].x) * pow(2, (int)ceil(((double)depth + 1)/ 2)) / TSP_GRID_SIZE);

}

/**
 * @brief Orders the node in the path based on a grid scheme
 * 
 * @param path The path to order
*/
void tsp_f2opt_partition_path(int* path) {

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
        list[i] = (tsp_entry){path[i], mapping[tsp_f2opt_calculate_block(7, path[i]) * 16 + tsp_f2opt_calculate_block(6, path[i])] + ((double)rand())/RAND_MAX};

    qsort((void*)list, (size_t)tsp_inst.nnodes, sizeof(tsp_entry), tsp_compare_entries);

    for (int i = 0; i < tsp_inst.nnodes; i++) path[i] = list[i].key;

    if (tsp_verbose >= 100) tsp_check_integrity(path, tsp_compute_path_cost(path), "algorithms.c - tsp_f2opt_partition_path - 1");

    safe_free(list);

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
int tsp_f2opt_swap(int* path, const int start, const int end) {

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
int tsp_f2opt_split(const int* path, int depth, const int start, const int end) {

    int start_block = tsp_f2opt_calculate_block(depth, path[start]);

    for (int i = start + 1; i < end; i++)
        if (tsp_f2opt_calculate_block(depth, path[i]) != start_block) return i;

    return start;

}

/**
 * @brief (THREAD SPECIFIC) Applies swaps till no better swap can be found
 * 
 * @param params tsp_mt_parameters type structure using: ->t_index, ->s_node, ->e_node, ->seed
*/
void* tsp_f2opt_block(void* params) {

    int* path = ((tsp_mt_parameters*)params) -> path;
    int start = ((tsp_mt_parameters*)params) -> s_node;
    int end = ((tsp_mt_parameters*)params) -> e_node;
    unsigned int* seed = ((tsp_mt_parameters*)params)->seed;

    while (tsp_f2opt_swap(path, start, end) > 0 && time_elapsed() < tsp_env.time_limit);

    if (end - start > 10) {
        double _ = 0;
        for (int j = 0; j < 5; j++)
            tsp_vns_3kick(path, &_, start, end, seed);

        while (tsp_f2opt_swap(path, start, end) > 0 && time_elapsed() < tsp_env.time_limit);
    }

    safe_free(seed);
    tsp_free_thread(((tsp_mt_parameters*)params)->t_index);

    return NULL;

}

/**
 * @brief Apply f2opt
 * 
 * @param path An array to store the solution into
 * 
 * @return The cost of the solution found
*/
double tsp_f2opt(int* path) {

    for (int i = 0; i < tsp_inst.nnodes; i++) path[i] = i;

    tsp_f2opt_partition_path(path);
    
    int* partitions[TSP_F2OPT_MAX_DEPTH];
    for (int i = 0; i < TSP_F2OPT_MAX_DEPTH; i++) partitions[i] = (int*)calloc(pow(2, i+1)+1, sizeof(int));

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
            unsigned int* seed = (unsigned int*)malloc(1 * sizeof(unsigned int));
            *seed = rand();

            params[thread].t_index = thread;
            params[thread].path = path;
            params[thread].s_node = partitions[i][j];
            params[thread].e_node = partitions[i][j+1];
            params[thread].seed = seed;
            
            pthread_create(&tsp_threads[thread], NULL, tsp_f2opt_block, (void*)&params[thread]);

        }
        tsp_wait_all_threads();
    }

    for (int i = 0; i < TSP_F2OPT_MAX_DEPTH; i++) safe_free(partitions[i]);

    double cost = tsp_compute_path_cost(path);
    tsp_2opt(path, &cost, tsp_find_2opt_best_swap);

    return cost;

}

void tsp_solve_g2opt() {

    if (tsp_env.g2opt_f2opt) {

        int* path = (int*) malloc(tsp_inst.nnodes * sizeof(int));
        double cost = tsp_f2opt(path);
        double time = time_elapsed();
        tsp_check_best_sol(path, cost, time);
        if (time > tsp_env.time_limit) tsp_env.status = 1;
        safe_free(path);
        return;

    }

    int (*swap_function)(int*, double*);
    switch (tsp_env.g2opt_swap_pol) {
        case 1: swap_function = tsp_find_2opt_swap; break;
        case 2: swap_function = tsp_find_2opt_best_swap; break;
        default: raise_error("Something went wrong choosing the swap policy.\n");
    }

    tsp_mt_parameters params[N_THREADS];

    for (int i = 0; i < tsp_inst.nnodes; i++) {

        if (time_elapsed() > tsp_env.time_limit) break;

        int thread = tsp_wait_for_thread();

        params[thread].t_index = thread;
        params[thread].s_node = i;
        params[thread].swap_function = swap_function;

        pthread_create(&tsp_threads[thread], NULL, tsp_greedy_from_node, (void*)&params[thread]);

    }

    tsp_wait_all_threads();

    if (time_elapsed() > tsp_env.time_limit) tsp_env.status = 1;

}

#pragma endregion


#pragma region METAHEURISTIC_TABU

/**
 * @brief (THREAD SPECIFIC) Finds the best 2opt swap checking and updating the tabu list after performing it
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

    tsp_check_best_sol(path, *cost, time_elapsed());

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

    while (time < tsp_env.time_limit) 
        tsp_find_tabu_swap(path, cost, t_index);

    tsp_free_thread(t_index);

    return NULL;

}

void tsp_solve_tabu() {
    
    tsp_mt_parameters params[N_THREADS];
    double cost[N_THREADS];

    tsp_allocate_tabu_space();

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

    for (int thread = 0; thread < N_THREADS; thread++) safe_free(params[thread].path);

    tsp_env.status = 1;

}

#pragma endregion


#pragma region METAHEURISTIC_VNS

/**
 * @brief (THREAD SPECIFIC) Performs a variable number of kicks and then fixes the solution using the 2opt algorithm
 * 
 * @param params tsp_mt_parameters type structure using: ->t_index, -> path, -> cost, -> seed
*/
void* tsp_vns_kicknsolve(void* params) {

    int* path = ((tsp_mt_parameters*)params)->path;
    double* cost = ((tsp_mt_parameters*)params)->cost;
    int t_index = ((tsp_mt_parameters*)params)->t_index;
    unsigned int* seed = ((tsp_mt_parameters*)params)->seed;

    for (int i = 0; i < t_index/4 + 1; i++) tsp_vns_3kick(path, cost, 0, tsp_inst.nnodes, seed);
    tsp_2opt(path, cost, tsp_find_2opt_best_swap);

    safe_free(seed);
    tsp_free_thread(t_index);

    return NULL;

}

/**
 * @brief (MULTITHREAD) Performs multiple kicks in parallel when we reach local optima
 * 
 * @param path The path to improve
 * @param cost The cost associated to the path
*/
void tsp_vns_multi_kicknsolve(int* path, double* cost) {

    tsp_mt_parameters params[N_THREADS];
    int* multi_kick_paths[N_THREADS];
    double multi_kick_cost[N_THREADS];

    for (int i = 0; i < N_THREADS; i++) multi_kick_paths[i] = (int*)calloc(tsp_inst.nnodes, sizeof(int));

    while (time_elapsed() < tsp_env.time_limit) {

        for (int i = 0; i < N_THREADS; i++) {

            int thread = tsp_wait_for_thread();
            unsigned int* seed = (unsigned int*)malloc(1 * sizeof(unsigned int));
            *seed = rand();

            for (int j = 0; j < tsp_inst.nnodes; j++) multi_kick_paths[thread][j] = path[j];
            multi_kick_cost[thread] = *cost;

            params[thread].t_index = thread;
            params[thread].path = multi_kick_paths[thread];
            params[thread].cost = &multi_kick_cost[thread];
            params[thread].seed = seed;

            pthread_create(&tsp_threads[thread], NULL, tsp_vns_kicknsolve, (void*)&params[thread]);

        }

        tsp_wait_all_threads();

        // choose best solution found -> improve that one
        int min_thread = 0;
        double min_cost = multi_kick_cost[min_thread];

        for (int i = 1; i < N_THREADS; i++)
            if (multi_kick_cost[i] < min_cost) {
                min_thread = i;
                min_cost = multi_kick_cost[min_thread];
            }

        for (int i = 0; i < tsp_inst.nnodes; i++) path[i] = multi_kick_paths[min_thread][i];
        *cost = min_cost;

        tsp_check_best_sol(path, *cost, time_elapsed());

    }

    for (int i = 0; i < N_THREADS; i++) safe_free(multi_kick_paths[i]);

}

void tsp_solve_vns() {

    int* path = (int*)calloc(tsp_inst.nnodes, sizeof(int));
    double cost = 0;
    switch (tsp_env.vns_fvns) {
        case 0:
            cost = tsp_greedy_path_from_node(path, rand() % tsp_inst.nnodes);
            tsp_2opt(path, &cost, tsp_find_2opt_best_swap);
            break;
        case 1:
            cost = tsp_f2opt(path);
            break;
        default: raise_error("Error choosing vns options.\n");
    }

    tsp_check_best_sol(path, cost, time_elapsed());

    tsp_vns_multi_kicknsolve(path, &cost);

    safe_free(path);

    tsp_env.status = 1;

}

#pragma endregion


#pragma region EXACT_CPLEX

/**
 * @brief Solve using cplex
 * 
 * @param env cplex environment
 * @param lp cplex lp
 * @param xstar xstar cplex's solution
 * @param ncomp number of components
 * @param comp list containing the component index of each node
 * @param succ successors type list containing the solution found by cplex
 * @param cost cost of the solution found
 * 
 * @return int 0 if the model was solved before the time limit, 1 if an intermediate solution has been found but cplex couldn't end, 2 if no solution has been found, 3 if the problem has been proven infeasible, 4 if the execution has been terminated by the user and a solution has been found, 5 if the execution has been terminated by the user and no solution has been found, 6 if cplex didn't even have the time to start
*/
int tsp_cplex_solve(CPXENVptr env, CPXLPptr lp, double* xstar, int* ncomp, int* comp, int* succ, double* cost) {

    int cpxerror = 0;

    // check time limit
    if (tsp_env.time_limit - time_elapsed() < TSP_EPSILON) {
        print_warn("The time limit is too short for this algorithm.\n");
        return 6;
    }

    // set the time limit
    cpxerror = CPXsetdblparam(env, CPXPARAM_TimeLimit, tsp_env.time_limit - time_elapsed());
    if (cpxerror) raise_error("CPXsetdblparam error (%d).\n", cpxerror);

    // solve the model using cplex
    cpxerror = CPXmipopt(env, lp);
    if (cpxerror) raise_error("CPXmipopt error (%d).\n", cpxerror);
    
    // get the output status (time limit (intermediate solution found / not found), infeasible)
    int output, status = CPXgetstat(env, lp);
    switch ( status ) {
        case CPXMIP_TIME_LIM_FEAS:      // exceeded time limit, found intermediate solution
            output = 1;
            break;
        case CPXMIP_TIME_LIM_INFEAS:    // exceeded time limit, no intermediate solution found
            return 2;
        case CPXMIP_INFEASIBLE:         // proven to be unfeasible
            return 3;
        case CPXMIP_ABORT_FEAS:         // terminated by user, found solution
            output = 4;
            break;
        case CPXMIP_ABORT_INFEAS:       // terminated by user, not found solution
            return 5;
        case CPXMIP_OPTIMAL_TOL:        // found optimal within the tollerance
            output = 0;
            break;
        case CPXMIP_OPTIMAL:            // found optimal
            output = 0;
            break;
        default:                        // unhandled status
            raise_error("Unhandled cplex status: %d.\n", status);
    }

    print_warn("--------------------- This is the output: %d.\n", output);

    // convert xstar to successors type list (save into succ)
	cpxerror = CPXgetx(env, lp, xstar, 0, CPXgetnumcols(env, lp)-1);
    if (cpxerror) raise_error("CPXgetx() error (%d).\n", cpxerror);

    *cost = tsp_cplex_compute_xstar_cost(xstar);
    tsp_cplex_decompose_xstar(xstar, comp, succ, ncomp);
    
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
            raise_error("Callback called for wrong reason: %d.\n", context_id);

    }

    // cplex error code
	return 1;

}

void tsp_solve_cplex() {

    if (!tsp_env.cplex_benders && !tsp_env.cplex_can_cb && !tsp_env.cplex_patching) print_warn("Neither benders, candidate callback or set. Solution might contain cycles.\n");

    int cpxerror;
    CPXENVptr env = NULL;
	CPXLPptr lp = NULL;

    tsp_cplex_init(&env, &lp, &cpxerror);

    int ncols = CPXgetnumcols(env, lp);
    
    // set callback function
    if (tsp_env.cplex_can_cb || tsp_env.cplex_rel_cb) {

        CPXLONG context_id = 0;

        if (tsp_env.cplex_can_cb) {
            if (tsp_verbose >= 10) print_info("Adding callback function for candidate solutions to cplex.\n");
            context_id = context_id | CPX_CALLBACKCONTEXT_CANDIDATE;
        }
        if (tsp_env.cplex_rel_cb) {
            if (tsp_verbose >= 10) print_info("Adding callback function for relaxation to cplex.\n");
            context_id = context_id | CPX_CALLBACKCONTEXT_RELAXATION;
        }

        cpxerror = CPXcallbacksetfunc(env, lp, context_id, tsp_cplex_callback, (void*)&tsp_inst.nnodes);
        if (cpxerror) raise_error("CPXcallbacksetfunc() error (%d).\n", cpxerror);

    }

    // add a starting heuristic to cplex
    if (tsp_env.cplex_mipstart) {

        int* path = (int*) malloc(tsp_inst.nnodes * sizeof(int));
        double cost = tsp_f2opt(path);

        tsp_check_best_sol(path, cost, time_elapsed());

        if (tsp_verbose >= 100) print_info("Finished f2opt (cost: %10.4f), starting cplex.\n", cost);

        if (tsp_verbose >= 10) print_info("Using an heuristic as mipstart for cplex.\n");

        int* index = (int *) malloc(ncols * sizeof(int));
        double* value = (double *) malloc(ncols * sizeof(double));
        int effortlevel = CPX_MIPSTART_NOCHECK;
        int izero = 0;
        tsp_cplex_path_to_ind_val(ncols, path, index, value);
        cpxerror = CPXaddmipstarts(env, lp, 1, tsp_inst.nnodes, &izero, index, value, &effortlevel, NULL);
        if (cpxerror) raise_error("CPXaddmipstarts() error (%d).\n", cpxerror);

        safe_free(value);
        safe_free(index);
        safe_free(path);

    }

    // space for data structures
    double* xstar = (double*) calloc(ncols, sizeof(double));
    int ncomp = 1;
    int* comp = (int*) calloc(tsp_inst.nnodes, sizeof(int));
    int* succ = (int*) calloc(tsp_inst.nnodes, sizeof(int));
    double cost = 0;

    // solve with cplex
    int ret = -1;
    if (tsp_env.cplex_benders) {

        if (tsp_verbose >= 10) print_info("Starting bender's loop.\n");

        int iter = 0;
        while (time_elapsed() < tsp_env.time_limit) {

            ret = tsp_cplex_solve(env, lp, xstar, &ncomp, comp, succ, &cost);

            if (tsp_verbose >= 100) print_info("Iteration number: %4d - connected components: %4d - current incumbent: %15.4f\n", iter++, ncomp, cost);

            if (ncomp == 1 || ret > 0) break;

            tsp_cplex_add_sec(env, lp, &ncomp, comp, succ);

            if (ncomp != 1 && ret != 4 && tsp_env.cplex_patching) {

                tsp_cplex_patching(xstar, &ncomp, comp, succ, &cost);
                //TODO: Give to cplex the patched solution
                
            }

        }

    } else
        ret = tsp_cplex_solve(env, lp, xstar, &ncomp, comp, succ, &cost);

    // apply patching if I have time left
    if (ret != 3 && ncomp != 1 && tsp_env.cplex_patching)  {

        if (ret == 0) tsp_cplex_patching(xstar, &ncomp, comp, succ, &cost);
        else print_warn("Couldn't apply patching due to time limit.\n");

    }

    // status from tsp_cplex_solve:
    /*
    0 : OK
    1 : time limit, found solution
    2 : time limit, didn't find solution
    3 : infeasible
    */

    // handle tsp_cplex_solve return status
    switch (ret) {
        case 0:
            tsp_cplex_check_best_sol(ncomp, comp, succ, cost);
            break;

        case 1:
            print_warn("cplex exceeded the time limit.\n");
            tsp_cplex_check_best_sol(ncomp, comp, succ, cost);            //TODO: Find a way to improve the solution within the time limit
            tsp_env.status = 1;
            break;

        case 2:
            print_warn("cplex couldn't find any feasible solution within the time limit.\n");
            if (tsp_env.cplex_mipstart) {
                print_info("Using the solution for the mipstart.\n");
                tsp_cplex_check_best_sol(ncomp, comp, succ, cost);
                tsp_env.status = 1;
            } else
                tsp_env.status = 3;
            break;

        case 3:
            tsp_env.status = 4;
            break;

        case 4:
            print_warn("cplex stopped by the user.\n");
            tsp_cplex_check_best_sol(ncomp, comp, succ, cost);            //TODO: Find a way to improve the solution within the time limit
            break;

        case 5:
            print_warn("cplex stopped by the user and couldn't find any solution.\n");
            if (tsp_env.cplex_mipstart) {
                print_info("Using the solution for the mipstart.\n");
                tsp_cplex_check_best_sol(ncomp, comp, succ, cost);
            }
            tsp_env.status = 2;
            break;

        case 6:
            if (tsp_env.cplex_mipstart)
                print_info("Using the solution for the mipstart.\n");
            
            break;

        default:
            raise_error("Unexpected return code from cplex_solve (%d).\n", ret);
    }
    
    tsp_cplex_close(env, lp, xstar, comp, succ);

}

#pragma endregion