#include "../include/inst_gen.h"
#include "../include/tsp.h"
#include "../include/mincut.h"

#pragma region GLOBALS DEFINITIONS

int tsp_verbose = TSP_DEFAULT_VERBOSE;

tsp_environment tsp_env;
tsp_instance    tsp_inst;

#pragma endregion


#pragma region MEMORY MANAGEMENT

/**
 * @brief Dinamically allocate the space for the costs list
*/
void tsp_allocate_costs_space() {

    if (!tsp_inst.nnodes) raise_error("The nnodes variable hasn't been assigned yet.");

    tsp_inst.costs = (double*)calloc(tsp_inst.nnodes * tsp_inst.nnodes, sizeof(double));

}

/**
 * @brief Dinamically allocate the space for the sort_edges list
*/
void tsp_allocate_sort_edges_space() {

    if (!tsp_inst.nnodes) raise_error("The nnodes variable hasn't been assigned yet.");

    tsp_inst.sort_edges = (int*)calloc(tsp_inst.nnodes * (tsp_inst.nnodes - 1), sizeof(int));

}

/**
 * @brief Dinamically allocate the space for the best solution
*/
void tsp_allocate_best_sol_space() {

    if (!tsp_inst.nnodes) raise_error("The nnodes variable hasn't been assigned yet.");

    tsp_inst.best_solution = (int*)calloc(tsp_inst.nnodes, sizeof(int));
    
}

void tsp_allocate_tabu_space() {

    for (int thread = 0; thread < N_THREADS; thread++) {
        tsp_env.tabu_tables[thread].list = (tsp_tabu_entry*)calloc(tsp_inst.nnodes, sizeof(tsp_tabu_entry));
        for (int i = 0; i < tsp_inst.nnodes; i++) { tsp_env.tabu_tables[thread].list[i] = (tsp_tabu_entry){-1, -1, -1, -1}; }
    }

}

void tsp_allocate_coords_space() {

    if (!tsp_inst.nnodes) raise_error("The nnodes variable hasn't been assigned yet.");

    tsp_inst.coords = (tsp_pair*)calloc(tsp_inst.nnodes, sizeof(tsp_pair));

}

void tsp_free_all() {

    if (tsp_inst.coords != NULL) { free(tsp_inst.coords); tsp_inst.coords = NULL; }
    if (tsp_inst.costs != NULL) { free(tsp_inst.costs); tsp_inst.costs = NULL; }
    if (tsp_inst.sort_edges != NULL) { free(tsp_inst.sort_edges); tsp_inst.sort_edges = NULL; }
    if (tsp_inst.best_solution != NULL) { free(tsp_inst.best_solution); tsp_inst.best_solution = NULL; }

    for (int thread = 0; thread < N_THREADS; thread++)
        if (tsp_env.tabu_tables[thread].list != NULL) { free(tsp_env.tabu_tables[thread].list); tsp_env.tabu_tables[thread].list = NULL; }

    pthread_mutex_destroy(&tsp_mutex_update_sol);

    if (tsp_inst.comp != NULL) { free(tsp_inst.comp); tsp_inst.comp = NULL; }
    if (tsp_inst.succ != NULL) { free(tsp_inst.succ); tsp_inst.succ = NULL; }

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

            if (checked_cost < min_cost - TSP_EPSILON) raise_error("SORT_EDGES INTEGRITY COMPROMISED\n");
            
            min_cost = checked_cost;
        }
    }

    if (tsp_verbose >= 1000) print_info("sort_edges integrity check passed.\n");

}

//TODO: Parallelize
void tsp_precompute_sort_edges() {

    tsp_allocate_sort_edges_space();
    tsp_entry* list = (tsp_entry*)calloc(tsp_inst.nnodes, sizeof(tsp_entry));

    for (int i = 0; i < tsp_inst.nnodes; i++) {

        for (int j = 0; j < tsp_inst.nnodes; j++) {  //saving the entries to be sorted
            list[j].key = j;
            list[j].value = tsp_get_edge_cost(i, j);    //considering only the costs of the edges leaving node i
        }

        qsort((void*)list, (size_t)tsp_inst.nnodes, sizeof(tsp_entry), tsp_compare_entries); //sort by cost of the edge
        
        for (int j = 1; j < tsp_inst.nnodes; j++)
            tsp_inst.sort_edges[i * (tsp_inst.nnodes - 1) + j-1] = list[j].key;    //populate the ith row with the nodes ordered by increasing distance

    }

    if (list != NULL) { free(list); list = NULL; }

    if (tsp_verbose >= 100) tsp_check_sort_edges_integrity();

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


#pragma region ALGORITHMS TOOLS

int tsp_find_alg() {

    if (!strcmp(tsp_env.alg_type, TSP_PARSING_GREEDY))  return 0;
    if (!strcmp(tsp_env.alg_type, TSP_PARSING_G2OPT))   return 1;
    if (!strcmp(tsp_env.alg_type, TSP_PARSING_TABU))    return 2;
    if (!strcmp(tsp_env.alg_type, TSP_PARSING_VNS))     return 3;
    if (!strcmp(tsp_env.alg_type, TSP_PARSING_CPLEX))   return 4;

    raise_error("Error choosing the algorithm to use.\n");

}

int tsp_compare_entries(const void* arg1, const void* arg2) {

    double diff = ((tsp_entry*)arg1)->value - ((tsp_entry*)arg2)->value;
    return (fabs(diff) < TSP_EPSILON) ? 0 : ((diff < 0) ? -1 : 1);

}

//[DONE]
void tsp_check_best_sol(const int* path, const double cost, const double time) {

    if (cost > tsp_inst.best_cost + TSP_EPSILON) return;

    pthread_mutex_lock(&tsp_mutex_update_sol);

    if (cost < tsp_inst.best_cost - TSP_EPSILON) {

        for (int i = 0; i < tsp_inst.nnodes; i++) tsp_inst.best_solution[i] = path[i];
        tsp_inst.best_cost = cost;
        tsp_inst.best_time = time;

        if (tsp_verbose >= 10) print_info("New best solution : %15.4f\n", cost);

    }

    pthread_mutex_unlock(&tsp_mutex_update_sol);

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
 * @return The new tenure
 */
double tsp_dinamic_tenue(int counter) {

    return tsp_env.tabu_tenure_a * sin((double)counter * tsp_env.tabu_tenure_f) + tsp_env.tabu_tenure;

}

int tsp_check_tabu(int t_index, int from, int to) {

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

void tsp_add_tabu(int t_index, int from, int to) {

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

double tsp_get_edge_cost(int i, int j) {

    return tsp_inst.costs[i * tsp_inst.nnodes + j];

}

double tsp_succ_to_path(const int* succ, int* path) {

    int* visited = (int*)calloc(tsp_inst.nnodes, sizeof(int));

    if (visited != NULL) { free(visited); visited = NULL; }

    for (int i=0, current_node=0; i<tsp_inst.nnodes; i++, current_node=succ[current_node]) path[i] = current_node;
    double cost = 0.0;
    for (int i = 0; i < tsp_inst.nnodes - 1; i++) cost += tsp_get_edge_cost(path[i], path[i+1]);
    cost += tsp_get_edge_cost(path[tsp_inst.nnodes - 1], path[0]);

    if (tsp_verbose >= 100) tsp_check_integrity(path, cost, "tsp.c - tsp_succ_to_path.\n");

    return cost;

}

#pragma endregion


#pragma region CPLEX

/**
 * @brief Converts coordinates to cplex x_pos
 * 
 * @param i Row coordinate
 * @param j Col coordinate
 * 
 * @return The index used by cplex to locate the edge (i, j)
*/
int tsp_cplex_coords_to_xpos(const int i, const int j) {

	if ( i == j ) raise_error("ERROR: i == j in xpos");
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
			double obj = tsp_get_edge_cost(i, j); // cost = distance   
			double lb = 0.0;
			double ub = 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) raise_error("ERROR: wrong CPXnewcols on x var.s");
    		if ( CPXgetnumcols(env, lp)-1 != tsp_cplex_coords_to_xpos(i,j) ) raise_error("ERROR: wrong position for x var.s");

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
		
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) ) raise_error("ERROR: wrong CPXaddrows [degree]");

	} 

    if (value != NULL) { free(value); value = NULL; }
    if (index != NULL) { free(index); index = NULL; }	

	if (tsp_verbose >= 100) CPXwriteprob(env, lp, "cplex_outputs/lp/model.lp", NULL);

	if (cname[0] != NULL) { free(cname[0]); cname[0] = NULL; }
	if (cname != NULL) { free(cname); cname = NULL; }

}

void tsp_cplex_init(CPXENVptr* env, CPXLPptr* lp, int* error) {

    *env = CPXopenCPLEX(error);
	*lp = CPXcreateprob(*env, error, "TSP");

    // set cplex log file
    CPXsetdblparam(*env, CPX_PARAM_SCRIND, CPX_OFF);
    char cplex_log_file[100];
    sprintf(cplex_log_file, "%s/%lu-%d-%s.log", TSP_CPLEX_LOG_FOLDER, tsp_env.seed, tsp_inst.nnodes, tsp_env.alg_type);
    remove(cplex_log_file);
    if ( CPXsetlogfilename(*env, cplex_log_file, "w") ) raise_error("CPXsetlogfilename error.\n");

    // build cplex model
    tsp_cplex_build_model(*env, *lp);

    // create lp file from cplex model
    char cplex_lp_file[100];
    sprintf(cplex_lp_file, "%s/%lu-%d-%s.lp", TSP_CPLEX_LP_FOLDER, tsp_env.seed, tsp_inst.nnodes, tsp_env.alg_type);
    if ( CPXwriteprob(*env, *lp, cplex_lp_file, NULL) ) raise_error("CPXwriteprob error\n");

}

void tsp_cplex_compute_xstar_cost(double* xstar, double* cost) {

    *cost = 0;

    for ( int i = 0; i < tsp_inst.nnodes; i++ ) for ( int j = i+1; j < tsp_inst.nnodes; j++ )
        if ( xstar[tsp_cplex_coords_to_xpos(i,j)] > TSP_CPLEX_ZERO_THRESHOLD ) *cost += tsp_get_edge_cost(i, j);

}

void tsp_cplex_add_sec(CPXENVptr env, CPXLPptr lp, int* ncomp, int* comp, int* succ) {

    if ((*ncomp)==1) raise_error("ERROR: add_sec() error");

    int izero = 0;
    char** cname = (char**)calloc(1, sizeof(char *));	// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));

    for (int k=0; k<(*ncomp); k++) {
        
        int* comp_nodes = (int*) malloc(0);
        int comp_size = 0;

        int start_node;
        for (start_node=0; start_node<tsp_inst.nnodes && comp[start_node]!=k+1; start_node++);
        int current_node = start_node, j = 0;
        do {
            comp_nodes = (int*) realloc(comp_nodes, ++comp_size*sizeof(int));
            comp_nodes[comp_size-1] = current_node;
            current_node = succ[current_node];
        } while (start_node!=current_node);

        int* index = (int*) calloc(CPXgetnumcols(env, lp), sizeof(double));
        double* value = (double*) calloc(CPXgetnumcols(env, lp), sizeof(double));
        int nnz = 0;
        char sense = 'L';
        double rhs = comp_size-1;
        // per nome vincolo, ottenere numeri righe, 
        for (int i=0; i<comp_size; i++) {
            for (int j=i+1; j<comp_size; j++) {
                index[nnz] = tsp_cplex_coords_to_xpos(comp_nodes[i], comp_nodes[j]);
                value[nnz] = 1.0;
                nnz++;
            }
        }
        if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) ) raise_error("ERROR: wrong CPXaddrows [degree]");

        if (comp_nodes != NULL) { free(comp_nodes); comp_nodes = NULL; }
        if (index != NULL) { free(index); index = NULL; }
        if (value != NULL) { free(value); value = NULL; }

    }

    if (cname[0] != NULL) { free(cname[0]); }
    if (cname != NULL) { free(cname); cname = NULL; }

}

void tsp_cplex_patch_comp(int* ncomp, int* comp, int* succ, double* cost) {

    int starts[*ncomp];
    for (int i=0; i<*ncomp; i++) starts[i]=-1;

    // glue components together until exceeding time limit or we have only one component left
    while (/*time_elapsed() < tsp_time_limit &&*/ (*ncomp)!=1) {

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

void tsp_cplex_check_best_sol(const int ncomp, const int* comp, const int* succ) {

    double time = time_elapsed();
    
    if (ncomp == 1) {

        int* solution = (int*) malloc(tsp_inst.nnodes * sizeof(int));
        double cost = tsp_succ_to_path(succ, solution);

        if (tsp_verbose >= 100) tsp_check_integrity(solution, cost, "tsp.c: tsp_cplex_check_best_sol - 1");

        tsp_check_best_sol(solution, cost, time);

        if (solution != NULL) { free(solution); solution = NULL; }

    } else {

        tsp_inst.best_time = time;

        tsp_inst.ncomp = ncomp;
        if (tsp_inst.comp == NULL) tsp_inst.comp = (int*)calloc(tsp_inst.nnodes, sizeof(int));
        if (tsp_inst.succ == NULL) tsp_inst.succ = (int*)calloc(tsp_inst.nnodes, sizeof(int));
        
        for (int i = 0; i < tsp_inst.nnodes; i++) { tsp_inst.comp[i] = comp[i]; tsp_inst.succ[i] = succ[i]; }

    }

}

void tsp_cplex_path_to_xstar(const int ncols, const int* path, int* indexes, double* values) {

    int k = 0;

    for (int i = 0; i < tsp_inst.nnodes-1; i++) {

        int from = path[i], to = path[i+1];
        int xpos = tsp_cplex_coords_to_xpos(from, to);
        indexes[k] = xpos;
        values[k++] = 1.0;

    }

    indexes[k] = tsp_cplex_coords_to_xpos(path[tsp_inst.nnodes - 1], path[0]);
    indexes[k++] = 1;

    if (k != tsp_inst.nnodes) raise_error("Something went wrong inside tsp_cplex_path_to_xstar.\n");

}

void tsp_cplex_decompose_xstar(const double* xstar, int* comp, int* succ, int* ncomp) {
    
    //initialize comp, succ
    *ncomp = 0;
	for ( int i = 0; i < tsp_inst.nnodes; i++ ) {
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
				if ( i!=j && xstar[tsp_cplex_coords_to_xpos(i,j)] > TSP_CPLEX_ZERO_THRESHOLD && comp[j] == -1 ) {
                    
					succ[i] = j;
					i = j;
					done = 0;
					break;

				}
			}
		}	
		succ[i] = start;

	}

}

int tsp_cplex_callback_candidate(CPXCALLBACKCONTEXTptr context, const int nnodes) {

    // get candidate point
    int ncols = nnodes * (nnodes - 1) / 2;
    double* xstar = (double*) malloc(ncols * sizeof(double));
    double objval = CPX_INFBOUND; if ( CPXcallbackgetcandidatepoint(context, xstar, 0, ncols-1, &objval) ) raise_error("CPXcallbackgetcandidatepoint() error.\n");
    
    if (objval == CPX_INFBOUND) raise_error("CPXcallbackgetcandidatepoint() error, no candidate objval returned.\n");

    // space for data structures
    int* comp = (int*) calloc(tsp_inst.nnodes, sizeof(int));
    int* succ = (int*) calloc(tsp_inst.nnodes, sizeof(int));
    int ncomp = 0;

    // convert xstar to comp and succ
    tsp_cplex_decompose_xstar(xstar, comp, succ, &ncomp);
    if(xstar != NULL) { free(xstar); xstar = NULL; }

    if (ncomp == 1) {   // feasible solution (only one tour)
    
        if(comp != NULL) { free(comp); comp = NULL; } 
        if(succ != NULL) { free(succ); succ = NULL; }

        // user info
        double lower_bound = CPX_INFBOUND; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_BND, &lower_bound);
        double incumbent = CPX_INFBOUND; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent);
        double gap = (1 - lower_bound/incumbent) * 100;

        //FIXME: This prints local info... not global incumbent
        //  : The first incumbent is CPX_INFBOUND if I don't give it the initial heuristic
        if (tsp_verbose >= 100) print_info("found feasible solution   -   lower_bound: %15.4f   -   incumbent: %15.4f   -   gap: %6.2f%c.\n", lower_bound, incumbent, gap, '%');

        return 0;

    }

    if (tsp_verbose >= 100) print_info("adding SEC    -   number of SEC: %d.\n", ncomp);

    // add as many SEC as connected components
    const char sense = 'L';
    const int izero = 0;
    int* index = (int*) calloc(ncols, sizeof(int));
    double* value = (double*) calloc(ncols, sizeof(double));

    for(int k = 1; k <= ncomp; k++) {

        int nnz = 0;
        double rhs = -1.0;

        for(int i = 0; i < tsp_inst.nnodes; i++) {

            if(comp[i]!=k) continue;
            rhs++;
            for(int j = i+1; j < tsp_inst.nnodes; j++){
                if(comp[j]!=k) continue;
                index[nnz]=tsp_cplex_coords_to_xpos(i,j);
                value[nnz]=1.0;
                nnz++;
            }

        }

        // reject candidate and add SEC
        if ( CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, value) ) raise_error("CPXcallbackrejectcandidate() error.\n");

    }

    // free the memory
    if(index != NULL) { free(index); index = NULL; }
    if(value != NULL) { free(value); value = NULL; }
    if(comp != NULL) { free(comp); comp = NULL; } 
    if(succ != NULL) { free(succ); succ = NULL; }

    return 0;

}

int tsp_concorde_callback_add_cplex_sec(double cut_value, int cut_nnodes, int* cut_index, void* userhandle) {

    CPXCALLBACKCONTEXTptr context = *(CPXCALLBACKCONTEXTptr*)userhandle;
    int cut_nedges = cut_nnodes * (cut_nnodes - 1) / 2;

	int* index = (int*) calloc(cut_nedges, sizeof(int));
	double* value = (double*) calloc(cut_nedges, sizeof(double));
	int nnz=0;

	for(int i=0; i<cut_nnodes; i++){
		for(int j=i+1; j<cut_nnodes; j++){
            index[nnz] = tsp_cplex_coords_to_xpos(cut_index[i], cut_index[j]);
            value[nnz] = 1.0;
            nnz++;
		}
	}

	const char sense = 'L';
	const double rhs = cut_nnodes - 1;
	const int izero = 0;
	const int purgeable = CPX_USECUT_PURGE;
	const int local = 0;

    if (CPXcallbackaddusercuts(context, 1, cut_nedges, &rhs, &sense, &izero, index, value, &purgeable, &local)) raise_error("CPXcallbackaddusercuts() error.\n");

	if (index != NULL) { free(index); index = NULL; }
	if (value != NULL) { free(value); value = NULL; }

	return 0;

}

int tsp_cplex_callback_relaxation(CPXCALLBACKCONTEXTptr context, const int nnodes) {

    int nodeuid = -1; if (CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODEUID, &nodeuid)) raise_error("CPXcallbackgetinfoint() error.\n");
    if (nodeuid % 10) return 0;

    if (tsp_verbose >= 1000) print_info("nodeuid: %d.\n", nodeuid);

    int ncols = (nnodes * (nnodes - 1) / 2);
	double* xstar = (double*) malloc(ncols * sizeof(double));  
	double objval = CPX_INFBOUND; if (CPXcallbackgetrelaxationpoint(context, xstar, 0, ncols-1, &objval)) raise_error("CPXcallbackgetcandidatepoint() error.\n");

    int ncomp = -1;
    int* elist = (int*) calloc(2 * ncols, sizeof(int));
    int k = 0;
    
    //TODO(ASK): elist should contain all edges in the graph?
    for (int i = 0; i < nnodes; i++) for (int j = i+1; j < nnodes; j++) {
        elist[k++] = i;
        elist[k++] = j;
    }

    int* comps = NULL;
    int* compscount = NULL;
    if(CCcut_connect_components(nnodes, ncols, elist, xstar, &ncomp, &compscount, &comps)) raise_error("CCcut_connect_components() error.\n");

    if (ncomp == 1) {
        
        if (tsp_verbose >= 100) print_info("Adding SEC for relaxation    -   number of SEC: %d.\n", 1);
        if(CCcut_violated_cuts(nnodes, ncols, elist, xstar, 1.9, tsp_concorde_callback_add_cplex_sec, (void*)&context)) raise_error("CCcut_violated_cuts() error.\n");
    
    } else {
        
        if (tsp_verbose >= 100) print_info("Adding SEC for relaxation    -   number of SEC: %d.\n", ncomp);
        
        int start = 0;
        for(int c=0; c<ncomp; ++c) {
            
            int* subtour = (int*) malloc(compscount[c] * sizeof(int));
            
            for(int i=0; i<compscount[c]; ++i) {
                subtour[i] = comps[i+start];
            }

            tsp_concorde_callback_add_cplex_sec(0, compscount[c], subtour, (void*)&context);

            start += compscount[c];

            if (subtour != NULL) { free(subtour); subtour = NULL; }

        }

    }

    if (comps != NULL) { free(comps); comps = NULL; }
    if (compscount != NULL) { free(compscount); compscount = NULL; }

    if (elist != NULL) { free(elist); elist = NULL; }
    if (xstar != NULL) { free(xstar); xstar = NULL; }

    return 0;
    
}

void tsp_cplex_close(CPXENVptr env, CPXLPptr lp, int* comp, int* succ) {

    // free memory
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
    if (comp != NULL) { free(comp); comp = NULL; }
    if (succ != NULL) { free(succ); succ = NULL; }

    // remove "clone<x>.log" files generated by cplex
    int file_number = 1;
    char clone_file[50];
    sprintf(clone_file, "clone0.log");
    remove(clone_file);
    do { sprintf(clone_file, "clone%d.log", file_number++); } while (!remove(clone_file));

}

#pragma endregion


#pragma region INITIALIZATIONS

/**
 * @brief Initializes the problem environment
*/
void tsp_init_env() {

    tsp_env.status = 0;

    strcpy(tsp_env.file_name, "NONE");
    tsp_env.seed = 0;
    strcpy(tsp_env.alg_type, "greedy");

    tsp_env.time_limit = (double)TSP_DEF_TL/1000;
    struct timeval tv; gettimeofday(&tv, NULL);
    tsp_env.time_start = ((double)tv.tv_sec)+((double)tv.tv_usec/1e+6);
    tsp_env.time_total = .0;

    tsp_env.tmp_choice = 0;
    
    tsp_env.g2opt_swap_pol = 0;
    tsp_env.tabu_tenure = TSP_DEF_TABU_TENURE;
    tsp_env.tabu_tenure_a = 0;
    tsp_env.tabu_tenure_f = .02;
    tsp_env.vns_fvns = 0;
    tsp_env.cplex_mipstart = 0;
    tsp_env.cplex_benders = 0;
    tsp_env.cplex_patching = 0;
    tsp_env.cplex_can_cb = 0;
    tsp_env.cplex_rel_cb = 0;

}

/**
 * @brief Initializes the problem instance
*/
void tsp_init_inst() {

    tsp_inst.nnodes = TSP_DEF_NNODES;
    tsp_inst.best_cost = INFINITY;
    tsp_inst.best_time = 0;
    tsp_inst.ncomp = 1;

}

void tsp_init_defs() {

    tsp_init_env();
    tsp_init_inst();
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

    char prefix[150], solution_file_name[500];

    if (tsp_env.seed > 0) 
        snprintf(prefix, sizeof(char)*150, "%lu_%d_%s", tsp_env.seed, tsp_inst.nnodes, tsp_env.alg_type);
    else
        snprintf(prefix, sizeof(char)*150, "%s_%s", tsp_env.file_name, tsp_env.alg_type);
    snprintf(solution_file_name, sizeof(char)*500, "%s/%s_%s", TSP_SOL_FOLDER, prefix, TSP_SOLUTION_FILE);  //where to save the file

    solution_file = fopen(solution_file_name, "w");

    if (solution_file == NULL) raise_error("Error writing the file for the solution.");

    fprintf(solution_file, "Algorithm: %s\n", tsp_env.alg_type);
    fprintf(solution_file, "Cost: %15.4f\n", tsp_inst.best_cost);
    fprintf(solution_file, "Time: %15.4fs\n", tsp_inst.best_time);
    fprintf(solution_file, "Total execution time: %15.4fs\n", tsp_env.time_total);
    if (tsp_env.status == 1) fprintf(solution_file, "The algorithm has exceeded the time limit and has been stopped.\n");
    else if (tsp_env.status == 2) fprintf(solution_file, "The algorithm has been terminated by the user.\n");
    else fprintf(solution_file, "The algorithm hasn't exceeded the time limit.\n");
    fprintf(solution_file, "--------------------\n");

    if (tsp_inst.ncomp>1) {

        for (int i=0; i<tsp_inst.ncomp; i++) {
            fprintf(solution_file, "LOOP %d:\n", i+1);
            int start_node, current_node;
            for (start_node=0; start_node<tsp_inst.nnodes && tsp_inst.comp[start_node]!=i+1; start_node++);
            current_node = start_node;
            do {
                fprintf(solution_file, "%4d %15.4f %15.4f\n", current_node, tsp_inst.coords[current_node].x, tsp_inst.coords[current_node].y);
                current_node = tsp_inst.succ[current_node];
            } while (current_node!=start_node);
            fprintf(solution_file, "%4d %15.4f %15.4f\n", start_node, tsp_inst.coords[start_node].x, tsp_inst.coords[start_node].y);
        }

    } else {

        for (int i = 0; i < tsp_inst.nnodes; i++)
            fprintf(solution_file, "%4d %15.4f %15.4f\n", tsp_inst.best_solution[i], tsp_inst.coords[tsp_inst.best_solution[i]].x, tsp_inst.coords[tsp_inst.best_solution[i]].y);
        fprintf(solution_file, "%4d %15.4f %15.4f\n", tsp_inst.best_solution[0], tsp_inst.coords[tsp_inst.best_solution[0]].x, tsp_inst.coords[tsp_inst.best_solution[0]].y);

    }

    fclose(solution_file);

}

void tsp_plot_solution() {

    int rows_read = 0, coord_files_number = 0;
    FILE *solution_file, *command_file;
    char plot_file_name[500], solution_file_name[500], solution_contents[100], gnuplot_command[500], prefix[150], gnuplot_title[1000];

    if (tsp_env.seed > 0) 
        snprintf(prefix, sizeof(char)*150, "%lu_%d_%s", tsp_env.seed, tsp_inst.nnodes, tsp_env.alg_type);
    else
        snprintf(prefix, sizeof(char)*150, "%s_%s", tsp_env.file_name, tsp_env.alg_type);

    snprintf(plot_file_name, sizeof(char)*500, "%s/%s_%s", TSP_SOL_FOLDER, prefix, TSP_PLOT_FILE);  //where to save the plot
    snprintf(solution_file_name, sizeof(char)*500, "%s/%s_%s", TSP_SOL_FOLDER, prefix, TSP_SOLUTION_FILE);  //where to read the file from

    solution_file = fopen(solution_file_name, "r");
    command_file = fopen(TSP_COMMAND_FILE, "w");

    if (solution_file == NULL || command_file == NULL) raise_error("Error with a file used to plot the solution.");

    // skip through the rows with the solution info
    while (rows_read < 6) if (fgetc(solution_file) =='\n') rows_read++;

    // build plot title with solution info
    snprintf(gnuplot_title, 1000, "Algorithm: %s; %d nodes; cost: %.4f; time: %.4fs", tsp_env.alg_type, tsp_inst.nnodes, tsp_inst.best_cost, tsp_inst.best_time);

    // copy nodes coordinates into coords_file
    if (tsp_inst.ncomp>1) {

        char coord_file_name[50];
        FILE *current_coord_file;
        while (fgets(solution_contents, 5, solution_file)) {
            if (!strcmp(solution_contents, "LOOP")) {
                if (coord_files_number) fclose(current_coord_file);
                sprintf(coord_file_name, "%d_%s", ++coord_files_number, TSP_COORDS_FILE);
                current_coord_file = fopen(coord_file_name, "w");
                while (fgetc(solution_file)!='\n');
            }
            else fprintf(current_coord_file, "%s", solution_contents);
        }
        fclose(current_coord_file);

    }
    else {

        FILE *coords_file = fopen(TSP_COORDS_FILE, "w");
        while (fgets(solution_contents, 100, solution_file)) fprintf(coords_file, "%s", solution_contents);
        fclose(coords_file);

    }
    
    // builds commands for gnuplot
    fprintf(command_file, "set term png\n");
    fprintf(command_file, "set output '%s'\n", plot_file_name);
    fprintf(command_file, "x=0.; y=0.\n");
    fprintf(command_file, "set title '%s'\n", gnuplot_title);
    if (tsp_inst.ncomp>1) {
        fprintf(command_file, "plot ");
        for (int i=0; i<coord_files_number; i++) {
            fprintf(command_file, "'%d_%s' u (x=$2):(y=$3) w lp lc rgb 'blue' title ''", i+1, TSP_COORDS_FILE);
            if (i!=coord_files_number-1) fprintf(command_file, ", ");
        }
        fprintf(command_file, "\n");
    }
    else {
        fprintf(command_file, "set xlabel 'Starting node highlighted as black point'\n");
        fprintf(command_file, "set label at %f, %f point pointtype 7 pointsize 2\n", tsp_inst.coords[tsp_inst.best_solution[0]].x, tsp_inst.coords[tsp_inst.best_solution[0]].y);
        fprintf(command_file, "plot '%s' u (x=$2):(y=$3) w lp lc rgb 'blue' title ''\n", TSP_COORDS_FILE);
    }

    fclose(solution_file);
    fclose(command_file);

    // execute gnuplot commands and remove all intermediate files
    snprintf(gnuplot_command, sizeof(char)*500, "gnuplot %s", TSP_COMMAND_FILE);
    system(gnuplot_command);
    remove(TSP_COORDS_FILE);
    for (int i=0; i<coord_files_number; i++) {
        char file[50];
        sprintf(file, "%d_%s", i+1, TSP_COORDS_FILE);
        remove(file);
    }
    remove(TSP_COMMAND_FILE);
    
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
    printf("Edge weight type: %s\n", TSP_DEF_EDGE_W_TYPE);

    printf("--------------------\n");

    printf("Algorithm: %s.\n", tsp_env.alg_type);

    if (tsp_env.g2opt_swap_pol) printf("Swap policy: %s.\n", ((tsp_env.g2opt_swap_pol == 1) ? "first swap" : "best swap"));
    if (!strcmp(tsp_env.alg_type, TSP_PARSING_TABU)) printf("Tabu tenure: %4d.\nTabu variability: %4d.\nTabu variability frequency: %10.4f.\n", tsp_env.tabu_tenure, tsp_env.tabu_tenure_a, tsp_env.tabu_tenure_f);
    if (tsp_env.vns_fvns) printf("Fast vns enabled.\n");
    if (tsp_env.cplex_mipstart) printf("Using a mipstart.\n");
    if (tsp_env.cplex_benders) printf("Using bender's loop.\n");
    if (tsp_env.cplex_patching) printf("Using patching.\n");
    if (tsp_env.cplex_can_cb) printf("Using candidate callback.\n");
    if (tsp_env.cplex_rel_cb) printf("Using relaxation callback.\n");
    if (tsp_env.tmp_choice) printf("Added temporary option.\n");
    
    printf("--------------------\n");

    if (tsp_verbose >= 100) {
        printf("Integrity checks enabled.\n");
        printf("--------------------\n");
    }

    if (tsp_verbose < 1000) return;

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
    printf("Cost: %15.4f\n", tsp_inst.best_cost);
    printf("Time:\t%15.4fs\n", tsp_inst.best_time);
    printf("Execution time: %8.4fs\n", tsp_env.time_total);
    if (tsp_verbose >= 500) {
        for (int i = 0; i < tsp_inst.nnodes; i++) printf("%d -> ", tsp_inst.best_solution[i]);
        printf("%d\n", tsp_inst.best_solution[0]);
    }
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
        case 5:
            break;
        case -5:
            printf("No solution has been found within the time limit.\n");
            break;
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

    if (visited != NULL) { free(visited); visited = NULL; }

    double c_cost = tsp_compute_path_cost(path);
    if (error == 0 && fabs(c_cost - cost) > TSP_EPSILON) error = 3;

    if (error >= 1) {
        print_warn("INTEGRITY COMPROMISED - error_code: %d ----- %s\n", error, message);
        if (error == 1) print_warn("Non-existent node in path.\n");
        else if (error == 2) print_warn("Double node in path.\n");
        else if (error == 3) print_warn("Cost: %.10f, Checked cost: %.10f, Difference: %.10f, Threshold: %.10f\n", cost, c_cost, fabs(c_cost - cost), TSP_EPSILON);
        else print_warn("Unknown error.\n");
        exit(1);
    }

    if (tsp_verbose >= 1000) print_info("Integrity check passed.\n");

}

#pragma endregion


#pragma region USEFUL METHODS

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