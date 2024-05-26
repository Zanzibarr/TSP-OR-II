#include "../include/inst_gen.h"
#include "../include/tsp.h"
#include "../include/mincut.h"

#pragma region GLOBALS DEFINITIONS

tsp_environment tsp_env;
tsp_instance    tsp_inst;
tsp_statistics  tsp_stat;

#pragma endregion


#pragma region MEMORY MANAGEMENT

/**
 * @brief Dinamically allocate the space for the costs list
*/
void tsp_allocate_costs_space() {

    if (!tsp_inst.nnodes) raise_error("Error in tsp_allocate_coords_space: The nnodes variable hasn't been assigned yet.\n");

    tsp_inst.costs = (double*)malloc(tsp_inst.nnodes * tsp_inst.nnodes * sizeof(double));

}

/**
 * @brief Dinamically allocate the space for the sort_edges list
*/
void tsp_allocate_sort_edges_space() {

    if (!tsp_inst.nnodes) raise_error("Error in tsp_allocate_sort_edges_space: The nnodes variable hasn't been assigned yet.\n");

    tsp_inst.sort_edges = (int*)malloc(tsp_inst.nnodes * (tsp_inst.nnodes - 1) * sizeof(int));

}

/**
 * @brief Dinamically allocate the space for the best solution
*/
void tsp_allocate_best_sol_space() {

    if (!tsp_inst.nnodes) raise_error("Error in tsp_allocate_best_sol_space: The nnodes variable hasn't been assigned yet.\n");

    tsp_inst.solution_succ = (int*)malloc(tsp_inst.nnodes * sizeof(int));
    
}

void tsp_allocate_tabu_space() {

    for (int thread = 0; thread < N_THREADS; thread++) {
        tsp_env.tabu_tables[thread].list = (tsp_tabu_entry*)malloc(tsp_inst.nnodes * sizeof(tsp_tabu_entry));
        for (int i = 0; i < tsp_inst.nnodes; i++) { tsp_env.tabu_tables[thread].list[i] = (tsp_tabu_entry){-1, -1, -1, -1}; }
    }

}

void tsp_allocate_coords_space() {

    if (!tsp_inst.nnodes) raise_error("Error in tsp_allocate_coords_space: The nnodes variable hasn't been assigned yet.\n");

    tsp_inst.coords = (tsp_pair*)malloc(tsp_inst.nnodes * sizeof(tsp_pair));

}

void tsp_free_all() {

    safe_free(tsp_inst.coords);
    safe_free(tsp_inst.costs);
    safe_free(tsp_inst.sort_edges);
    safe_free(tsp_inst.solution_succ);

    for (int thread = 0; thread < N_THREADS; thread++) safe_free(tsp_env.tabu_tables[thread].list);

    pthread_mutex_destroy(&tsp_mutex_update_sol);

}

#pragma endregion


#pragma region PRECOMPUTING

/**
 * @brief Debugging tool used to check whether the sorting of the sort_edges list is correct
*/
void tsp_check_sort_edges_integrity() {
    
    for (int i = 0; i < tsp_inst.nnodes; i++) {
        double min_cost = 0;
        for (int j = 0; j < tsp_inst.nnodes - 1; j++) {
            double checked_cost = tsp_get_edge_cost(i, tsp_inst.sort_edges[i * (tsp_inst.nnodes - 1) + j]);

            if (checked_cost < min_cost - TSP_EPSILON) raise_error("INTEGRITY CHECK: Error in tsp_check_sort_edges_integrity: integrity compromised.\n");
            
            min_cost = checked_cost;
        }
    }

}

//TODO: Parallelize
void tsp_precompute_sort_edges() {

    tsp_allocate_sort_edges_space();
    tsp_entry* list = (tsp_entry*)malloc(tsp_inst.nnodes * sizeof(tsp_entry));

    for (int i = 0; i < tsp_inst.nnodes; i++) {

        for (int j = 0; j < tsp_inst.nnodes; j++) {  //saving the entries to be sorted
            list[j].key = j;
            list[j].value = tsp_get_edge_cost(i, j);    //considering only the costs of the edges leaving node i
        }

        qsort((void*)list, (size_t)tsp_inst.nnodes, sizeof(tsp_entry), tsp_compare_entries); //sort by cost of the edge
        
        for (int j = 1; j < tsp_inst.nnodes; j++)
            tsp_inst.sort_edges[i * (tsp_inst.nnodes - 1) + j-1] = list[j].key;    //populate the ith row with the nodes ordered by increasing distance

    }

    safe_free(list);

    if (tsp_env.effort_level >= 100) tsp_check_sort_edges_integrity();

}

/**
 * @brief Compute the euclidian distance between the points i and j
 * 
 * @param i First node
 * @param j Second node
 * 
 * @return The euclidian distance between the two points
*/
double tsp_compute_distance(const int i, const int j) {
    
    return (int)(sqrt(pow(tsp_inst.coords[i].x - tsp_inst.coords[j].x, 2) + pow(tsp_inst.coords[i].y - tsp_inst.coords[j].y, 2)) + .5);
    
}

//TODO: Parallelize
void tsp_precompute_costs() {

    tsp_allocate_costs_space();

    for (int i = 0; i < tsp_inst.nnodes; i++) for (int j = 0; j < tsp_inst.nnodes; j++)
        tsp_inst.costs[i * tsp_inst.nnodes + j] = tsp_compute_distance(i, j);

}

#pragma endregion


#pragma region CONVERSIONS

int tsp_convert_coord_to_xpos(const int i, const int j) {

	if ( i == j ) raise_error("INTEGRITY CHECK: Error in tsp_convert_coord_to_xpos: i == j.\n");
	if ( i > j ) return tsp_convert_coord_to_xpos(j,i);

	return i * tsp_inst.nnodes + j - (( i + 1 ) * ( i + 2 )) / 2;

}

void tsp_convert_path_to_indval(const int ncols, const int* path, int* ind, double* val) {

    if (ncols <= 0 || path == NULL || ind == NULL || val == NULL) raise_error("Error in tsp_convert_path_to_indval.\n");

    if (tsp_env.effort_level >= 150) print_info("Converting path to indval.\n");

    double t_start = time_elapsed();

    int k = 0;

    for (int i = 0; i < tsp_inst.nnodes-1; i++) {

        int xpos = tsp_convert_coord_to_xpos(path[i], path[i+1]);
        ind[k] = xpos;
        val[k++] = 1.0;

    }

    ind[k] = tsp_convert_coord_to_xpos(path[tsp_inst.nnodes - 1], path[0]);
    val[k++] = 1.0;

    // Integrity check
    if (tsp_env.effort_level >= 100) {
        if (k != tsp_inst.nnodes) raise_error("INTEGRITY CHECK: Error in tsp_convert_path_to_indval: k != nnodes (%d != %d).\n", k, tsp_inst.nnodes);
        for (int e = 0; e < k; e++) if (ind[e] < 0 || ind[e] >= ncols || val[e] != 1.0) raise_error("INTEGRITY CHECK: Error in tsp_convert_path_to_indval: filling ind or val (%d - %f).\n", ind[e], val[e]);
    }

    pthread_mutex_lock(&tsp_mutex_update_stat);
    tsp_stat.time_for_conversions += time_elapsed() - t_start;
    pthread_mutex_unlock(&tsp_mutex_update_stat);

}

void tsp_convert_comp_to_cutindval(const int kcomp, const int ncomp, const int ncols, const int* comp, int* ind, double* val, int* nnz, double* rhs) {

    if (ncols <= 0 || kcomp <= 0 || kcomp > ncomp || comp == NULL || ind == NULL || val == NULL || nnz == NULL || rhs == NULL) raise_error("Error in tsp_convert_comp_to_cutindval.\n");

    if (tsp_env.effort_level >= 150) print_info("Converting comp to indval.\n");

    double t_start = time_elapsed();

    *nnz = 0;
    *rhs = -1.0;

    for(int i = 0; i < tsp_inst.nnodes; i++) {

        if(comp[i] != kcomp) continue;
        (*rhs)++;

        for(int j = i+1; j < tsp_inst.nnodes; j++){
            if(comp[j] != kcomp) continue;
            ind[*nnz]=tsp_convert_coord_to_xpos(i,j);
            val[*nnz]=1.0;
            (*nnz)++;
        }

    }

    // Integrity check
    if (tsp_env.effort_level >= 100) {
        if (*nnz < 0 || *nnz > ncols) raise_error("INTEGRITY CHECK: Error in tsp_convert_comp_to_cutindval: calculating nnz (%d).\n", *nnz);
        for (int e = 0; e < *nnz; e++) if (ind[e] < 0 || ind[e] >= ncols || val[e] != 1.0) raise_error("INTEGRITY CHECK: Error in tsp_convert_comp_to_cutindval: filling ind or val (%d - %f).\n", ind[e], val[e]);
    }

    pthread_mutex_lock(&tsp_mutex_update_stat);
    tsp_stat.time_for_conversions += time_elapsed() - t_start;
    pthread_mutex_unlock(&tsp_mutex_update_stat);

}

void tsp_convert_succ_to_solindval(const int* succ, const int ncols, int* ind, double* val) {

    for (int i = 0; i < ncols; i++) { ind[i] = i; val[i] = 0; }

    for (int i = 0; i < tsp_inst.nnodes; i++) {
        int xpos = tsp_convert_coord_to_xpos(i, succ[i]);
        ind[xpos] = xpos;
        val[xpos] = 1.0;
    }
    
}

void tsp_convert_xstar_to_compsucc(const double* xstar, int* comp, int* ncomp, int* succ) {

    if (xstar == NULL || comp == NULL || ncomp == NULL || succ == NULL) raise_error("Error in tsp_convert_xstar_to_compsucc.\n");

    if (tsp_env.effort_level >= 150) print_info("Converting xstar to compsucc.\n");

    double t_start = time_elapsed();

    //initialize comp, succ
    *ncomp = 0;
	for (int i = 0; i < tsp_inst.nnodes; i++) {
		succ[i] = -1;
		comp[i] = -1;
	}
	
    //convert xstar to comp, succ
	for ( int start = 0; start < tsp_inst.nnodes; start++ ) {
		
        if ( comp[start] >= 0 ) continue;

		(*ncomp)++;
		int i = start;
		int done = 0;
		while ( !done ) {
			comp[i] = (*ncomp);
			done = 1;
			for ( int j = 0; j < tsp_inst.nnodes; j++ ) {
				if ( i!=j && xstar[tsp_convert_coord_to_xpos(i,j)] > TSP_CPLEX_ZERO_THRESHOLD && comp[j] == -1 ) {
                    
					succ[i] = j;
					i = j;
					done = 0;
					break;

				}
			}
		}	
		succ[i] = start;

	}
    
    // Integrity check
    if (tsp_env.effort_level >= 100) {
        int* check = (int*)calloc(tsp_inst.nnodes, sizeof(int));
        for (int i = 0; i < tsp_inst.nnodes; i++) {
            if (check[succ[i]]) raise_error("INTEGRITY CHECK: Error in tsp_convert_xstar_to_compsucc: double node in succ.\n");
            check[succ[i]] = 1;
        }
        safe_free(check);
    }

    pthread_mutex_lock(&tsp_mutex_update_stat);
    tsp_stat.time_for_conversions += time_elapsed() - t_start;
    pthread_mutex_unlock(&tsp_mutex_update_stat);

}

void tsp_convert_succ_to_path(const int* succ, const double ncomp, int* path) {

    if (succ == NULL || ncomp != 1 || path == NULL) raise_error("Error in tsp_convert_succ_to_path.\n");

    if (tsp_env.effort_level >= 150) print_info("Converting succ to path.\n");

    double t_start = time_elapsed();

    for (int i=0, current_node=0; i<tsp_inst.nnodes; i++, current_node=succ[current_node]) path[i] = current_node;

    // Integrity check
    if (tsp_env.effort_level >= 100) tsp_check_integrity(path, tsp_compute_path_cost(path), "tsp.c - tsp_succ_to_path.\n");

    pthread_mutex_lock(&tsp_mutex_update_stat);
    tsp_stat.time_for_conversions += time_elapsed() - t_start;
    pthread_mutex_unlock(&tsp_mutex_update_stat);
    
}

void tsp_convert_path_to_succ(const int* path, int* succ) {

    if (succ == NULL || path == NULL) raise_error("Error in tsp_convert_path_to_succ.\n");

    if (tsp_env.effort_level >= 150) print_info("Converting path to succ.\n");

    double t_start = time_elapsed();

    int current_node = path[tsp_inst.nnodes - 1];
    for (int i = 0; i < tsp_inst.nnodes; i++) {
        succ[current_node] = path[i];
        current_node = path[i];
    }

    //Integrity check
    if (tsp_env.effort_level >= 100) {
        int* tmp_path = (int*) malloc(tsp_inst.nnodes * sizeof(int));
        int ncomp = 1;
        tsp_convert_succ_to_path(succ, ncomp, tmp_path);    //this has an integrity check inside
        safe_free(tmp_path);
    }

    pthread_mutex_lock(&tsp_mutex_update_stat);
    tsp_stat.time_for_conversions += time_elapsed() - t_start;
    pthread_mutex_unlock(&tsp_mutex_update_stat);
    
}

void tsp_convert_xstar_to_elistnxstar(const double* xstar, const int nnodes, int* elist, double* nxstar, int* nedges) {

    if (xstar == NULL || nnodes == 0|| elist == NULL || nxstar == NULL || nedges == NULL) raise_error("Error in tsp_convert_xstar_to_elistnxstar.\n");

    if (tsp_env.effort_level >= 150) print_info("Converting xstar to elistnxstar.\n");

    double t_start = time_elapsed();

    int k = 0;

    for (int i = 0; i < nnodes; i++) for (int j = i+1; j < nnodes; j++) {
        int xpos = tsp_convert_coord_to_xpos(i, j);
        if (xstar[xpos] > TSP_EPSILON) {
            elist[k++] = i;
            elist[k++] = j;
            nxstar[(*nedges)++] = xstar[xpos];
        }
    }

    // Integrity check
    if (tsp_env.effort_level >= 100) {
        if (*nedges <= 0) raise_error("INTEGRITY CHECK: Error in tsp_convert_xstar_to_elistnxstar: calculating the number of edges.\n");
        for (int e = 0; e < *nedges; e++) if (nxstar[e] <= TSP_EPSILON || nxstar[e] > 1 + TSP_EPSILON) raise_error("INTEGRITY CHECK: Error in tsp_convert_xstar_to_elistnxstar: passing xstar to nxstar (%f).\n", nxstar[e]);
        for (int i = 0; i < k; i++) if (elist[i] < 0 || elist[i] >= tsp_inst.nnodes) raise_error("INTEGRITY CHECK: Error in tsp_convert_xstar_to_elistnxstar: computing elist.\n");
    }

    pthread_mutex_lock(&tsp_mutex_update_stat);
    tsp_stat.time_for_conversions += time_elapsed() - t_start;
    pthread_mutex_unlock(&tsp_mutex_update_stat);
    
}

void tsp_convert_cutindex_to_indval(const int* cut_index, const int cut_nnodes, int* ind, double* val, int* nnz) {

    if (cut_index == NULL || cut_nnodes <= 0 || cut_nnodes > tsp_inst.nnodes || ind == NULL || val == NULL || nnz == NULL) raise_error("Error in tsp_convert_cutindex_to_indval.\n");

    if (tsp_env.effort_level >= 150) print_info("Converting cutindex to indval.\n");

    double t_start = time_elapsed();

	for(int i=0; i<cut_nnodes; i++){
		for(int j=i+1; j<cut_nnodes; j++){
            ind[*nnz] = tsp_convert_coord_to_xpos(cut_index[i], cut_index[j]);
            val[*nnz] = 1.0;
            (*nnz)++;
		}
	}

    // Integrity check
    if (tsp_env.effort_level >= 100) {
        if (*nnz < 0 || *nnz > cut_nnodes * (cut_nnodes - 1) / 2) raise_error("INTEGRITY CHECK: Error in tsp_convert_cutindex_to_indval: calculating nnz (%d).\n", *nnz);
        for (int e = 0; e < *nnz; e++) if (ind[e] < 0 || ind[e] >= tsp_inst.nnodes * (tsp_inst.nnodes - 1) / 2 || val[e] != 1.0) raise_error("INTEGRITY CHECK: Error in tsp_convert_cutindex_to_indval: filling ind or val (%d - %f).\n", ind[e], val[e]);
    }

    pthread_mutex_lock(&tsp_mutex_update_stat);
    tsp_stat.time_for_conversions += time_elapsed() - t_start;
    pthread_mutex_unlock(&tsp_mutex_update_stat);

}

#pragma endregion


#pragma region ALGORITHMS TOOLS

int tsp_find_alg() {

    if (!strcmp(tsp_env.alg_type, TSP_PARSING_GREEDY))  return 0;
    if (!strcmp(tsp_env.alg_type, TSP_PARSING_G2OPT))   return 1;
    if (!strcmp(tsp_env.alg_type, TSP_PARSING_TABU))    return 2;
    if (!strcmp(tsp_env.alg_type, TSP_PARSING_VNS))     return 3;
    if (!strcmp(tsp_env.alg_type, TSP_PARSING_CPLEX))   return 4;

    raise_error("Error in tsp_find_alg: choosing the algorithm to use.\n");
    return -1;

}

int tsp_compare_entries(const void* arg1, const void* arg2) {

    double diff = ((tsp_entry*)arg1)->value - ((tsp_entry*)arg2)->value;
    return (fabs(diff) < TSP_EPSILON) ? 0 : ((diff < 0) ? -1 : 1);

}

void tsp_check_best_sol(const int* path, const int* succ, const int* ncomp, const double* cost, const double time) {

    int* sol = (int*) malloc(tsp_inst.nnodes * sizeof(int));
    int sol_ncomp = 0;
    double sol_cost = -1;
    
    if (path != NULL) {

        if (succ != NULL || succ != NULL || ncomp != NULL) raise_error("Error in tsp_check_best_sol: path.\n");

        sol_ncomp = 1;

        tsp_convert_path_to_succ(path, sol);

    }

    if (succ != NULL && ncomp != NULL) {

        if (path != NULL) raise_error("Error in tsp_check_best_sol: succ.\n");

        for (int i = 0; i < tsp_inst.nnodes; i++) sol[i] = succ[i];
        sol_ncomp = *ncomp;

    }

    if (cost == NULL) sol_cost = tsp_compute_succ_cost(sol);
    else sol_cost = *cost;

    // Integrity check
    if (tsp_env.effort_level >= 100) {
        if (sol == NULL || sol_ncomp == 0 || sol_cost < 0) raise_error("INTEGRITY CHECK: Error in tsp_check_best_sol: integrity.\n");
        if (sol_ncomp == 1) {
            int* tmp_path = (int*) malloc(tsp_inst.nnodes * sizeof(int));
            int ncomp = 1;
            tsp_convert_succ_to_path(sol, ncomp, tmp_path); // this has an integrity check inside
            safe_free(tmp_path);
        } else {
            int* check = (int*) calloc(tsp_inst.nnodes, sizeof(int));
            for (int i = 0; i < tsp_inst.nnodes; i++) {
                if (check[succ[i]]) raise_error("INTEGRITY CHECK: Error in tsp_check_best_sol: double node in succ.\n");
                check[succ[i]] = 1;
            }
            safe_free(check);
        }
    }

    // statistics
    if (sol_ncomp == 1) {
        pthread_mutex_lock(&tsp_mutex_update_stat);
        tsp_stat.n_solutions_found++;
        pthread_mutex_unlock(&tsp_mutex_update_stat);
    }

    // If I have ONE component and you have a better ONE component solution, you win
    // If I have MORE components and you have ONE component, you win
    // If I have MORE components and you have a better MORE components solution, you win
    if (
        (
            sol_ncomp == 1 && tsp_inst.ncomp == 1 && sol_cost >= tsp_inst.best_cost - TSP_EPSILON
        ) || (
            sol_ncomp != 1 && (
                tsp_inst.ncomp == 1
                ||
                (tsp_inst.ncomp != 1 && sol_cost >= tsp_inst.best_cost - TSP_EPSILON)
            )
        )
    ) return;

    pthread_mutex_lock(&tsp_mutex_update_sol);

    if (
        (sol_ncomp == 1 && 
            (
                tsp_inst.ncomp == 1 && sol_cost < tsp_inst.best_cost - TSP_EPSILON  // I have a better ONE component solution than you, I win
                ||
                tsp_inst.ncomp != 1 // I have a ONE component solution and you don't, I win
            )
        ) || (
            sol_ncomp > 1 && tsp_inst.ncomp != 1 && sol_cost < tsp_inst.best_cost - TSP_EPSILON // We both have a MORE components solutions, but mine is better, I win
        )
    ) {

        if (tsp_env.effort_level >= 1) print_info("New best solution found (cost: %15.4f).\n", sol_cost);

        for (int i = 0; i < tsp_inst.nnodes; i++) tsp_inst.solution_succ[i] = sol[i];
        tsp_inst.best_cost = sol_cost;
        tsp_inst.best_time = time;
        tsp_inst.ncomp = sol_ncomp;

        tsp_stat.n_solutions_improv++;

    }

    pthread_mutex_unlock(&tsp_mutex_update_sol);

    if (path != NULL) safe_free(sol);

}

void tsp_reverse(int* path, int start, int end) {

    while (start < end) {
        int c = path[start];
        path[start] = path[end];
        path[end] = c;
        start++; end--;
    }

}

/**
 * @brief Computes the dinamic tenure based on the counter
 * 
 * @param counter Counter used to calculate the tenure
 * 
 * @return The new tenure
 */
double tsp_dinamic_tenue(const int counter) {

    return tsp_env.tabu_tenure_a * sin((double)counter * tsp_env.tabu_tenure_f) + tsp_env.tabu_tenure;

}

int tsp_check_tabu(const int t_index, const int from, const int to) {

    if (from > to) return tsp_check_tabu(t_index, to, from);

    return 
        tsp_env.tabu_tables[t_index].list[from].counter_1 != -1 &&  //if I have a tabu saved from the "from" node
        (
            tsp_env.tabu_tables[t_index].list[from].node_1 == to && //if the node "to" creates a tabu "from"-"to" (first option)...
            tsp_env.tabu_tables[t_index].counter - tsp_env.tabu_tables[t_index].list[from].counter_1 < tsp_dinamic_tenue(tsp_env.tabu_tables[t_index].counter)   //... and I still remember that as a tabu 
            ||
            tsp_env.tabu_tables[t_index].list[from].node_2 == to && //if the node "to" creates a tabu "from"-"to" (second option)...
            tsp_env.tabu_tables[t_index].counter - tsp_env.tabu_tables[t_index].list[from].counter_2 < tsp_dinamic_tenue(tsp_env.tabu_tables[t_index].counter)   //... and I still remember that as a tabu 
        );

}

void tsp_add_tabu(const int t_index, const int from, const int to) {

    if (from > to) return tsp_add_tabu(t_index, to, from);

    if (tsp_env.tabu_tables[t_index].list[from].counter_1 != -1 && tsp_env.tabu_tables[t_index].counter - tsp_env.tabu_tables[t_index].list[from].counter_1 < tsp_dinamic_tenue(tsp_env.tabu_tables[t_index].counter)) {

        tsp_env.tabu_tables[t_index].list[from].node_2 = to;
        tsp_env.tabu_tables[t_index].list[from].counter_2 = tsp_env.tabu_tables[t_index].counter++;

    } else {

        tsp_env.tabu_tables[t_index].list[from].node_1 = to;
        tsp_env.tabu_tables[t_index].list[from].counter_1 = tsp_env.tabu_tables[t_index].counter++;

    }

}

double tsp_compute_path_cost(const int* path) {

    double cost = 0;

    for (int i = 0; i < tsp_inst.nnodes - 1; i++) cost += tsp_get_edge_cost(path[i], path[i+1]);
    cost += tsp_get_edge_cost(path[tsp_inst.nnodes - 1], path[0]);

    return cost;

}

double tsp_compute_xstar_cost(const double* xstar) {

    double cost = 0;

    for ( int i = 0; i < tsp_inst.nnodes; i++ ) for ( int j = i+1; j < tsp_inst.nnodes; j++ )
        if ( xstar[tsp_convert_coord_to_xpos(i,j)] > TSP_CPLEX_ZERO_THRESHOLD ) cost += tsp_get_edge_cost(i, j);

    return cost;

}

double tsp_compute_succ_cost(const int* succ) {

    double cost = 0;

    for (int i = 0; i < tsp_inst.nnodes; i++) cost += tsp_get_edge_cost(i, succ[i]);

    return cost;

}

double tsp_get_edge_cost(const int i, const int j) {

    return tsp_inst.costs[i * tsp_inst.nnodes + j];

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

void tsp_2opt(int* path, double* cost, int (*swap_function)(int*, double*)) {

    while ((*swap_function)(path, cost) > 0 && time_elapsed() < tsp_env.time_limit) {//repeat until I can't find new swaps that improve the cost of my solution
        // Integrity check
        if (tsp_env.effort_level >= 100) tsp_check_integrity(path, *cost, "tsp.c - tsp_2opt - 1");
    }

}

#pragma endregion


#pragma region INITIALIZATIONS

/**
 * @brief Initializes the problem environment
*/
void tsp_init_env() {

    tsp_env.status = 0;
    tsp_env.effort_level = TSP_DEF_VERBOSE;

    strcpy(tsp_env.file_name, "NONE");
    tsp_env.seed = 0;
    strcpy(tsp_env.alg_type, "greedy");

    tsp_env.time_limit = (double)TSP_DEF_TL/1000;
    tsp_env.cplex_terminate = 0;

    struct timeval tv; gettimeofday(&tv, NULL);
    tsp_env.time_start = ((double)tv.tv_sec)+((double)tv.tv_usec/1e+6);
    tsp_env.time_total = .0;

    tsp_env.noplot = 0;
    tsp_env.tmp_choice = 0;
    
    tsp_env.g2opt_swap_pol = 0;
    tsp_env.g2opt_f2opt = 0;
    tsp_env.tabu_tenure = TSP_DEF_TABU_TENURE;
    tsp_env.tabu_tenure_a = 0;
    tsp_env.tabu_tenure_f = .02;
    tsp_env.vns_fvns = 0;
    tsp_env.cplex_mipstart = 0;
    tsp_env.cplex_benders = 0;
    tsp_env.cplex_patching = 0;
    tsp_env.cplex_can_cb = 0;
    tsp_env.cplex_rel_cb = 0;
    tsp_env.cplex_cb_patching = 0;
    tsp_env.cplex_hard_fixing = 0;
    tsp_env.cplex_hard_fixing_pfix = TSP_DEFAULT_PFIX;
    tsp_env.cplex_local_branching = 0;
    tsp_env.cplex_local_branching_k = 100;

}

/**
 * @brief Initializes the problem instance
*/
void tsp_init_inst() {

    tsp_inst.nnodes = TSP_DEF_NNODES;
    tsp_inst.best_cost = INFINITY;
    tsp_inst.best_time = 0;
    tsp_inst.ncomp = 0;

}

void tsp_init_stat() {

    tsp_stat.n_solutions_found = 0;
    tsp_stat.n_solutions_improv = 0;

    tsp_stat.time_for_conversions = .0;
    tsp_stat.time_for_candidate_callback = .0;
    tsp_stat.time_for_relaxation_callback = .0;
    
}

void tsp_init_defs() {

    tsp_init_env();
    tsp_init_inst();
    tsp_init_stat();
    tsp_init_threads();

}

void tsp_init_solution() {

    // build solution
    if (tsp_env.seed > 0)   //if random
        tsp_gen_random_instance();
    else    //if from instance
        tsp_gen_instance_from_file();

    //precomputing
    tsp_precompute_costs();
    tsp_precompute_sort_edges();
    tsp_allocate_best_sol_space();

}

#pragma endregion


#pragma region SAVING FILES

void tsp_save_solution() {
    
    FILE *solution_file;

    // Stuff for the file name
    char prefix[200];

    char timestamp[50];
    time_t now = time (0);
    strftime (timestamp, 100, "%Y-%m-%d--%H:%M:%S", localtime (&now));

    if (tsp_env.seed > 0) 
        sprintf(prefix, "%lu_%d_%s_%s", tsp_env.seed, tsp_inst.nnodes, tsp_env.alg_type, timestamp);
    else
        sprintf(prefix, "%s_%s_%s", tsp_env.file_name, tsp_env.alg_type, timestamp);
    sprintf(tsp_env.solution_file, "%s/%s_%s", TSP_SOL_FOLDER, prefix, TSP_SOLUTION_FILE);  //where to save the file

    solution_file = fopen(tsp_env.solution_file, "w");
    if (solution_file == NULL) raise_error("Error in tsp_save_solution: writing the file for the solution.");

    // Printing to file the solution and other info
    fprintf(solution_file, "Algorithm: %s\n", tsp_env.alg_type);

    if (tsp_env.g2opt_swap_pol) fprintf(solution_file, "Swap policy: %s.\n", ((tsp_env.g2opt_swap_pol == 1) ? "first swap" : "best swap"));
    if (tsp_env.g2opt_f2opt) fprintf(solution_file, "Using f2opt.\n");
    if (!strcmp(tsp_env.alg_type, TSP_PARSING_TABU)) fprintf(solution_file, "Tabu tenure: %4d.\nTabu variability: %4d.\nTabu variability frequency: %10.4f.\n", tsp_env.tabu_tenure, tsp_env.tabu_tenure_a, tsp_env.tabu_tenure_f);
    if (tsp_env.vns_fvns) fprintf(solution_file, "Fast vns enabled.\n");
    if (tsp_env.cplex_mipstart) fprintf(solution_file, "Using a mipstart.\n");
    if (tsp_env.cplex_benders) fprintf(solution_file, "Using benders loop.\n");
    switch (tsp_env.cplex_patching) {
        case 1:
            fprintf(solution_file, "Using normal patching on solutions.\n");
            break;
        case 2:
            fprintf(solution_file, "Using greedy patching on solutions.\n");
            break;
        case 0: break;
        default:
            raise_error("Error in tsp_instance_info: choosing the patching function.\n");
    }
    switch (tsp_env.cplex_cb_patching) {
        case 1:
            fprintf(solution_file, "Using patching on candidate callback.\n");
            break;
        case 2:
            fprintf(solution_file, "Using patching on relaxation callback.\n");
            break;
        case 3:
            fprintf(solution_file, "Using patching on both candidate and relaxation callback.\n");
            break;
        case 0: break;
        default:
            raise_error("Error in tsp_instance_info: choosing the callback patching function.\n");
    }
    if (tsp_env.cplex_can_cb) fprintf(solution_file, "Using candidate callback.\n");
    if (tsp_env.cplex_rel_cb) fprintf(solution_file, "Using relaxation callback.\n");
    if (tsp_env.tmp_choice) fprintf(solution_file, "Added temporary option (%d).\n", tsp_env.tmp_choice);

    fprintf(solution_file, "Number of components: %4d\n", tsp_inst.ncomp);
    fprintf(solution_file, "Cost: %15.4f\n", tsp_inst.best_cost);
    fprintf(solution_file, "Time: %15.4fs\n", tsp_inst.best_time);
    fprintf(solution_file, "Total execution time: %15.4fs\n", tsp_env.time_total);
    fprintf(solution_file, "Number of feasible integer solutions found: %4d.\n", tsp_stat.n_solutions_found);
    fprintf(solution_file, "Number of solution improvements: %4d.\n", tsp_stat.n_solutions_improv);
    fprintf(solution_file, "Time lost for conversions: %10.8fs\n", tsp_stat.time_for_conversions);
    fprintf(solution_file, "Time spent in candidate callback: %10.8fs\n", tsp_stat.time_for_candidate_callback);
    fprintf(solution_file, "Time spent in relaxation callback: %10.8fs\n", tsp_stat.time_for_relaxation_callback);

    switch (tsp_env.status) {
        case 1:
            fprintf(solution_file, "The algorithm exceeded the time limit and has been stopped.\n");
            break;
        case 2:
            fprintf(solution_file, "The algorithm has been terminated by the user.\n");
            break;
        case 3:
            fprintf(solution_file, "cplex couldn't find any solution within the time limit.\n");
            break;
        case 4:
            fprintf(solution_file, "The problem has been proven to be infeasible.\n");
            break;
        case 0: break;
        default: raise_error("Error in tsp_save_solution: Unexpected status.\n");
    }
    
    fprintf(solution_file, "--------------------\n");

    for (int i = 0; i < tsp_inst.nnodes; i++) {
        int to = tsp_inst.solution_succ[i];
        fprintf(solution_file, "(%15.4f, %15.4f) -> (%15.4f, %15.4f)\n", tsp_inst.coords[i].x, tsp_inst.coords[i].y, tsp_inst.coords[to].x, tsp_inst.coords[to].y);
    }

    fclose(solution_file);

}

void tsp_plot_solution() {

    if (tsp_env.noplot) return;

    print_info("Plotting solution.\n");

    char command[600];

    sprintf(command, "python3 plotting/plot_solution.py %s", tsp_env.solution_file);
    system(command);

}

#pragma endregion


#pragma region DEBUGGING TOOLS

void tsp_instance_info() {

    printf("--------------------\n");
    printf("Type of Instance: %s\n", ((tsp_env.seed == 0) ? "from file" : "random"));
    if (tsp_env.seed == 0) printf("File name: %s\n", tsp_env.file_name);
    else printf("Seed: %lu\n", tsp_env.seed);
    printf("Time limit: %10.4fs\n", tsp_env.time_limit);
    printf("Number of nodes: %4d\n", tsp_inst.nnodes);
    printf("Edge weight type: %s\n", TSP_EDGE_W_TYPE);

    printf("--------------------\n");

    printf("Algorithm: %s.\n", tsp_env.alg_type);

    if (tsp_env.g2opt_swap_pol) printf("Swap policy: %s.\n", ((tsp_env.g2opt_swap_pol == 1) ? "first swap" : "best swap"));
    if (tsp_env.g2opt_f2opt) printf("Fast 2opt enabled.\n");
    if (!strcmp(tsp_env.alg_type, TSP_PARSING_TABU)) printf("Tabu tenure: %4d.\nTabu variability: %4d.\nTabu variability frequency: %10.4f.\n", tsp_env.tabu_tenure, tsp_env.tabu_tenure_a, tsp_env.tabu_tenure_f);
    if (tsp_env.vns_fvns) printf("Fast vns enabled.\n");
    if (tsp_env.cplex_mipstart) printf("Using a mipstart.\n");
    if (tsp_env.cplex_benders) printf("Using benders loop.\n");
    switch (tsp_env.cplex_patching) {
        case 1:
            printf("Using normal patching on solutions.\n");
            break;
        case 2:
            printf("Using greedy patching on solutions.\n");
            break;
        case 0: break;
        default:
            raise_error("Error in tsp_instance_info: choosing the patching function.\n");
    }
    switch (tsp_env.cplex_cb_patching) {
        case 1:
            printf("Using patching on candidate callback.\n");
            break;
        case 2:
            printf("Using patching on relaxation callback.\n");
            break;
        case 3:
            printf("Using patching on both candidate and relaxation callback.\n");
            break;
        case 0: break;
        default:
            raise_error("Error in tsp_instance_info: choosing the callback patching function.\n");
    }
    if (tsp_env.cplex_can_cb) printf("Using candidate callback.\n");
    if (tsp_env.cplex_rel_cb) printf("Using relaxation callback.\n");
    if (tsp_env.tmp_choice) printf("Added temporary option: %d.\n", tsp_env.tmp_choice);
    
    printf("--------------------\n");

    if (tsp_env.effort_level >= 100) {
        printf("Integrity checks enabled (all kind).\n");
        printf("\033[93m\033[1m[ WARN  ]:\033[0m Execution might be slowed down.\n");
        printf("--------------------\n");
    }

    if (tsp_env.effort_level < 1000) return;

    printf("NODES:\n");
    for (int i = 0; i < tsp_inst.nnodes; i++) printf("node[%4d]: (%15.4f, %15.4f)\n", i, tsp_inst.coords[i].x, tsp_inst.coords[i].y);
    printf("--------------------\n");
    printf("COSTS\n");
    for (int i = 0; i < tsp_inst.nnodes; i++) for (int j = 0; j < tsp_inst.nnodes; j++) printf("v%4d - v%4d : %15.4f\n", i, j, tsp_get_edge_cost(i, j));
    printf("--------------------\n");
    printf("MIN EDGES\n");
    for (int i = 0; i < tsp_inst.nnodes; i++) {
        for (int j = 0; j < tsp_inst.nnodes - 1; j++)
            printf("%4d(%15.4f) ", tsp_inst.sort_edges[i * (tsp_inst.nnodes - 1) + j], tsp_get_edge_cost(i, tsp_inst.sort_edges[i * (tsp_inst.nnodes - 1) + j]));
        printf("\n");
    }
    printf("--------------------\n");

}

void tsp_print_solution() {

    printf("--------------------\nBEST SOLUTION:\n");
    printf("Number of components: %4d\n", tsp_inst.ncomp);
    printf("Cost: %15.4f\n", tsp_inst.best_cost);
    printf("Time:\t%15.4fs\n", tsp_inst.best_time);
    printf("Execution time: %8.4fs\n", tsp_env.time_total);
    
    switch (tsp_env.status) {
        case 1:
            printf("The algorithm exceeded the time limit and has been stopped.\n");
            break;
        case 2:
            printf("The algorithm has been terminated by the user.\n");
            break;
        case 3:
            printf("cplex couldn't find any solution within the time limit.\n");
            break;
        case 4:
            printf("The problem has been proven to be infeasible.\n");
            break;
        case 0: break;
        default: raise_error("Error in tsp_print_solution: unexpected status.\n");
    }

    printf("--------------------\nSTATISTICS:\n");
    printf("Number of feasible integer solutions found: %4d.\n", tsp_stat.n_solutions_found);
    printf("Number of solution improvements: %4d.\n", tsp_stat.n_solutions_improv);
    printf("Time lost for conversions: %10.8fs\n", tsp_stat.time_for_conversions);
    printf("Time spent in candidate callback: %10.8fs\n", tsp_stat.time_for_candidate_callback);
    printf("Time spent in relaxation callback: %10.8fs\n", tsp_stat.time_for_relaxation_callback);
    if (tsp_env.effort_level >= 500) {
        for (int i = 0; i < tsp_inst.nnodes; i++) {
            int to = tsp_inst.solution_succ[i];
            printf("(%15.4f, %15.4f) -> (%15.4f, %15.4f)\n", tsp_inst.coords[i].x, tsp_inst.coords[i].y, tsp_inst.coords[to].x, tsp_inst.coords[to].y);
        }
    }
    printf("--------------------\n");

}

void tsp_check_integrity(const int* path, const double cost, const char* message) {

    int error = 0;
    int* visited = (int*)calloc(tsp_inst.nnodes, sizeof(int));

    for (int i = 0; i < tsp_inst.nnodes; i++) {
        if (path[i] < 0 || path[i] >= tsp_inst.nnodes) { error = 1; break; }
        if (visited[path[i]]) { error = 2; break; }
        visited[path[i]] += 1;
    }

    safe_free(visited);

    double c_cost = tsp_compute_path_cost(path);
    if (error == 0 && fabs(c_cost - cost) > TSP_EPSILON) error = 3;

    if (error >= 1) {
        print_warn("INTEGRITY COMPROMISED - error_code: %d ----- %s\n", error, message);
        if (error == 1) raise_error("INTEGRITY CHECK: Non-existent node in path.\n");
        else if (error == 2) raise_error("INTEGRITY CHECK: Double node in path.\n");
        else if (error == 3) raise_error("INTEGRITY CHECK: Cost: %.10f, Checked cost: %.10f, Difference: %.10f, Threshold: %.10f\n", cost, c_cost, fabs(c_cost - cost), TSP_EPSILON);
        else raise_error("INTEGRITY CHECK: Unknown error.\n");
    }

}

#pragma endregion


#pragma region USEFUL METHODS

void safe_free(void* ptr) {

    if (ptr != NULL) { free(ptr); ptr = NULL; }

}

void init_rand() { for (int i = 0; i < 100; i++) rand(); }

double rnd_coord() { return (double)rand()/RAND_MAX*TSP_GRID_SIZE; }

double time_elapsed() {

    struct timeval tv;

    gettimeofday(&tv, NULL);
    return ((double)tv.tv_sec)+((double)tv.tv_usec/1e+6) - tsp_env.time_start;

}

void print_info(const char* str, ...) {

    // initializing list pointer 
    va_list ptr; 
    va_start(ptr, str);

    fprintf(stdout, "\033[92m\033[1m[ INFO  ]:\033[0m Time: %10.4fs   -   ", time_elapsed());
  
    // char array to store token 
    char token[1000]; 
    // index of where to store the characters of str in 
    // token 
    int k = 0; 
  
    // parsing the formatted string 
    for (int i = 0; str[i] != '\0'; i++) { 
        token[k++] = str[i]; 
  
        if (str[i + 1] == '%' || str[i + 1] == '\0') { 
            token[k] = '\0'; 
            k = 0; 
            if (token[0] != '%') { 
                fprintf( 
                    stdout, "%s", 
                    token); // printing the whole token if 
                            // it is not a format specifier 
            } 
            else { 
                int j = 1; 
                char ch1 = 0; 
  
                // this loop is required when printing 
                // formatted value like 0.2f, when ch1='f' 
                // loop ends 
                while ((ch1 = token[j++]) < 58) { 
                } 
                // for integers 
                if (ch1 == 'i' || ch1 == 'd' || ch1 == 'u'
                    || ch1 == 'h') { 
                    fprintf(stdout, token, 
                            va_arg(ptr, int)); 
                } 
                // for characters 
                else if (ch1 == 'c') { 
                    fprintf(stdout, token, 
                            va_arg(ptr, int)); 
                } 
                // for float values 
                else if (ch1 == 'f') { 
                    fprintf(stdout, token, 
                            va_arg(ptr, double)); 
                } 
                else if (ch1 == 'l') { 
                    char ch2 = token[2]; 
  
                    // for long int 
                    if (ch2 == 'u' || ch2 == 'd'
                        || ch2 == 'i') { 
                        fprintf(stdout, token, 
                                va_arg(ptr, long)); 
                    } 
  
                    // for double 
                    else if (ch2 == 'f') { 
                        fprintf(stdout, token, 
                                va_arg(ptr, double)); 
                    } 
                } 
                else if (ch1 == 'L') { 
                    char ch2 = token[2]; 
  
                    // for long long int 
                    if (ch2 == 'u' || ch2 == 'd'
                        || ch2 == 'i') { 
                        fprintf(stdout, token, 
                                va_arg(ptr, long long)); 
                    } 
  
                    // for long double 
                    else if (ch2 == 'f') { 
                        fprintf(stdout, token, 
                                va_arg(ptr, long double)); 
                    } 
                } 
  
                // for strings 
                else if (ch1 == 's') { 
                    fprintf(stdout, token, 
                            va_arg(ptr, char*)); 
                } 
  
                // print the whole token 
                // if no case is matched 
                else { 
                    fprintf(stdout, "%s", token); 
                } 
            } 
        } 
    }
  
    // ending traversal 
    va_end(ptr);

}

void print_warn(const char* str, ...) {

    // initializing list pointer 
    va_list ptr; 
    va_start(ptr, str);

    fprintf(stdout, "\033[93m\033[1m[ WARN  ]:\033[0m Time: %10.4fs   -   ", time_elapsed());
  
    // char array to store token 
    char token[1000]; 
    // index of where to store the characters of str in 
    // token 
    int k = 0; 
  
    // parsing the formatted string 
    for (int i = 0; str[i] != '\0'; i++) { 
        token[k++] = str[i]; 
  
        if (str[i + 1] == '%' || str[i + 1] == '\0') { 
            token[k] = '\0'; 
            k = 0; 
            if (token[0] != '%') { 
                fprintf( 
                    stdout, "%s", 
                    token); // printing the whole token if 
                            // it is not a format specifier 
            } 
            else { 
                int j = 1; 
                char ch1 = 0; 
  
                // this loop is required when printing 
                // formatted value like 0.2f, when ch1='f' 
                // loop ends 
                while ((ch1 = token[j++]) < 58) { 
                } 
                // for integers 
                if (ch1 == 'i' || ch1 == 'd' || ch1 == 'u'
                    || ch1 == 'h') { 
                    fprintf(stdout, token, 
                            va_arg(ptr, int)); 
                } 
                // for characters 
                else if (ch1 == 'c') { 
                    fprintf(stdout, token, 
                            va_arg(ptr, int)); 
                } 
                // for float values 
                else if (ch1 == 'f') { 
                    fprintf(stdout, token, 
                            va_arg(ptr, double)); 
                } 
                else if (ch1 == 'l') { 
                    char ch2 = token[2]; 
  
                    // for long int 
                    if (ch2 == 'u' || ch2 == 'd'
                        || ch2 == 'i') { 
                        fprintf(stdout, token, 
                                va_arg(ptr, long)); 
                    } 
  
                    // for double 
                    else if (ch2 == 'f') { 
                        fprintf(stdout, token, 
                                va_arg(ptr, double)); 
                    } 
                } 
                else if (ch1 == 'L') { 
                    char ch2 = token[2]; 
  
                    // for long long int 
                    if (ch2 == 'u' || ch2 == 'd'
                        || ch2 == 'i') { 
                        fprintf(stdout, token, 
                                va_arg(ptr, long long)); 
                    } 
  
                    // for long double 
                    else if (ch2 == 'f') { 
                        fprintf(stdout, token, 
                                va_arg(ptr, long double)); 
                    } 
                } 
  
                // for strings 
                else if (ch1 == 's') { 
                    fprintf(stdout, token, 
                            va_arg(ptr, char*)); 
                } 
  
                // print the whole token 
                // if no case is matched 
                else { 
                    fprintf(stdout, "%s", token); 
                } 
            } 
        } 
    }
  
    // ending traversal 
    va_end(ptr);

}

void raise_error(const char* str, ...) {

    // initializing list pointer 
    va_list ptr; 
    va_start(ptr, str);

    fprintf(stdout, "\033[91m\033[1m[ ERROR ]:\033[0m Time: %10.4fs - ", time_elapsed());
  
    // char array to store token 
    char token[1000]; 
    // index of where to store the characters of str in 
    // token 
    int k = 0; 
  
    // parsing the formatted string 
    for (int i = 0; str[i] != '\0'; i++) { 
        token[k++] = str[i]; 
  
        if (str[i + 1] == '%' || str[i + 1] == '\0') { 
            token[k] = '\0'; 
            k = 0; 
            if (token[0] != '%') { 
                fprintf( 
                    stdout, "%s", 
                    token); // printing the whole token if 
                            // it is not a format specifier 
            } 
            else { 
                int j = 1; 
                char ch1 = 0; 
  
                // this loop is required when printing 
                // formatted value like 0.2f, when ch1='f' 
                // loop ends 
                while ((ch1 = token[j++]) < 58) { 
                } 
                // for integers 
                if (ch1 == 'i' || ch1 == 'd' || ch1 == 'u'
                    || ch1 == 'h') { 
                    fprintf(stdout, token, 
                            va_arg(ptr, int)); 
                } 
                // for characters 
                else if (ch1 == 'c') { 
                    fprintf(stdout, token, 
                            va_arg(ptr, int)); 
                } 
                // for float values 
                else if (ch1 == 'f') { 
                    fprintf(stdout, token, 
                            va_arg(ptr, double)); 
                } 
                else if (ch1 == 'l') { 
                    char ch2 = token[2]; 
  
                    // for long int 
                    if (ch2 == 'u' || ch2 == 'd'
                        || ch2 == 'i') { 
                        fprintf(stdout, token, 
                                va_arg(ptr, long)); 
                    } 
  
                    // for double 
                    else if (ch2 == 'f') { 
                        fprintf(stdout, token, 
                                va_arg(ptr, double)); 
                    } 
                } 
                else if (ch1 == 'L') { 
                    char ch2 = token[2]; 
  
                    // for long long int 
                    if (ch2 == 'u' || ch2 == 'd'
                        || ch2 == 'i') { 
                        fprintf(stdout, token, 
                                va_arg(ptr, long long)); 
                    } 
  
                    // for long double 
                    else if (ch2 == 'f') { 
                        fprintf(stdout, token, 
                                va_arg(ptr, long double)); 
                    } 
                } 
  
                // for strings 
                else if (ch1 == 's') { 
                    fprintf(stdout, token, 
                            va_arg(ptr, char*)); 
                } 
  
                // print the whole token 
                // if no case is matched 
                else { 
                    fprintf(stdout, "%s", token); 
                } 
            } 
        } 
    }
  
    // ending traversal 
    va_end(ptr);

    exit(1);

}

#pragma endregion