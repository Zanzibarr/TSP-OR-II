#include "../include/exact.h"

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

    int cpxerror = 0;

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
			cpxerror = CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname);
            if (cpxerror) raise_error("ERROR: wrong CPXnewcols on x var.s (%d).\n", cpxerror);
    		cpxerror = CPXgetnumcols(env, lp)-1 != tsp_cplex_coords_to_xpos(i,j);
            if (cpxerror) raise_error("ERROR: wrong position for x var.s (%d).\n", cpxerror);

		}

	} 

	// add degree constr.s 
	int* index = (int*)malloc(tsp_inst.nnodes * sizeof(int));
	double* value = (double*)malloc(tsp_inst.nnodes * sizeof(double));  
	
	// add the degree constraints
	for (int h = 0; h < tsp_inst.nnodes; h++) { // degree constraints

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
		
		cpxerror = CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]);
        if (cpxerror) raise_error("ERROR: wrong CPXaddrows [degree] (%d).\n", cpxerror);

	} 

    safe_free(value);
    safe_free(index);

	if (tsp_verbose >= 100) CPXwriteprob(env, lp, "cplex_outputs/lp/model.lp", NULL);

	if (cname[0] != NULL) free(cname[0]);
	safe_free(cname);

}

/**
 * @brief apply the "normal" patching to the cplex solution
 * 
 * @param ncomp number of components
 * @param comp list containing the component index of each node
 * @param succ successors type list containing the solution found by cplex
 * @param cost the cost of the solution
*/
void tsp_cplex_patch(int* ncomp, int* comp, int* succ, double* cost) {

    int starts[*ncomp];
    for (int i=0; i<(*ncomp); i++) starts[i]=-1;

    // glue components together until exceeding time limit or we have only one component left
    while ((*ncomp)!=1) {

        // edge to be removed expressed by first node in succ order
        int best_k1 = 0, best_k2 = 0, best_edge_k1 = 0, best_edge_k2 = 0;
        double best_delta = -INFINITY, delta_N, delta_R;
        // move stores how to glue components: given edges (i,succ[i]) and (j, succ[j]), replace them with:
        // - move='R' (Reverse): (i,j) and (succ[j],succ[i]) -> reverse k2 afterwards
        // - move='N' (Not reverse): (i,succ[j]) and (j,succ[i]) -> do not reverse k2 afterwards
        char best_move;

        // look for best patching move
        for (int k1=1; k1<=*ncomp; k1++) {

            int start_k1 = starts[k1-1];
            if (start_k1==-1) {
                for (start_k1=0; comp[start_k1]!=k1; start_k1++);
                starts[k1-1] = start_k1;
            }
            for (int k2=k1+1; k2<=*ncomp; k2++) {

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
        for (int i=best_k2-1; i<*ncomp-1; i++) starts[i] = starts[i+1];
        starts[*ncomp-1]=-1;

    }
    
    tsp_inst.ncomp = 1;

    //TODO: Update cost
    //TODO(maybe): Improve with 2opt

}

/**
 * @brief apply the "greedy" patching to the cplex solution
 * 
 * @param xstar cplex's xstar solution
 * @param ncomp number of components
 * @param comp list containing the component index of each node
 * @param succ successors type list containing the solution found by cplex
 * @param cost the cost of the solution
*/
void tsp_cplex_patch_greedy(const double* xstar, int* ncomp, int* comp, int* succ, double* cost) {

    int* patched = (int*) malloc(tsp_inst.nnodes * sizeof(int));

    for (int i = 0; i < tsp_inst.nnodes; i++) { patched[i] = -1; succ[i] = -1; comp[i] = 1; }
    *ncomp = 1;

    int start_node = 0;
    patched[0] = start_node;
    int current_node = start_node;

    // get patched solution with fractionary greedy
    for (int i = 1; i < tsp_inst.nnodes; i++) {

        double min = INFINITY;
        int min_j = -1;
        
        for (int j = 0; j < tsp_inst.nnodes; j++) {

            if (current_node == j) continue;
            if (succ[j] >= 0) continue;

            int xpos = tsp_cplex_coords_to_xpos(current_node, j);
            double edge_weight = tsp_get_edge_cost(current_node, j) * (1 - xstar[xpos]);

            if (edge_weight < min) {
                min = edge_weight;
                min_j = j;
            }

        }

        patched[i] = min_j;
        succ[current_node] = min_j;
        current_node = min_j;

    }

    *cost = tsp_compute_path_cost(patched);

    // Integrity check
    if (tsp_verbose >= 100) tsp_check_integrity(patched, *cost, "tsp.c: tsp_cplex_patch_greedy - 1");

    // improve patched solution with 2opt
    tsp_2opt(patched, cost, tsp_find_2opt_best_swap);

    // save patched solution
    current_node = patched[tsp_inst.nnodes - 1];
    for (int i = 0; i < tsp_inst.nnodes; i++) {
        succ[current_node] = patched[i];
        current_node = patched[i];
    }

    safe_free(patched);

    // Integrity check
    if (tsp_verbose >= 100) {
        int* check = (int*) calloc(tsp_inst.nnodes, sizeof(int));
        for (int i = 0; i < tsp_inst.nnodes; i++) {
            if (succ[i] < 0 || succ[i] >= tsp_inst.nnodes || check[succ[i]]) raise_error("Double node in (succ) patched solution with greedy (i: %d, succ[i]: %d, check[succ[i]]: %d).\n", i, succ[i], check[succ[i]]);
            check[succ[i]] = 1;
        }
        safe_free(check);
    }

}

/**
 * @brief function to add sec to cplex using concorde
 * 
 * @param cut_value concorde's cut value
 * @param cut_nnodes number of nodes in the cut
 * @param cut_index indexes of nodes in the cut
 * @param userhandle user data
 * 
 * @return concorde error code
*/
int tsp_concorde_callback_add_cplex_sec(double cut_value, int cut_nnodes, int* cut_index, void* userhandle) {

    CPXCALLBACKCONTEXTptr context = *(CPXCALLBACKCONTEXTptr*)userhandle;
    int cut_nedges = cut_nnodes * (cut_nnodes - 1) / 2;

    int cpxerror = 0;

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

    cpxerror = 0; cpxerror = CPXcallbackaddusercuts(context, 1, cut_nedges, &rhs, &sense, &izero, index, value, &purgeable, &local);
    if (cpxerror) raise_error("CPXcallbackaddusercuts() error (%d).\n", cpxerror);

	safe_free(value);
	safe_free(index);

	return 0;

}

void tsp_cplex_init(CPXENVptr* env, CPXLPptr* lp, int* cpxerror) {

    *env = CPXopenCPLEX(cpxerror);
	*lp = CPXcreateprob(*env, cpxerror, "TSP");

    if (*cpxerror) raise_error("CPX env or lp error (%d).\n", *cpxerror);

    // set cplex log file
    CPXsetdblparam(*env, CPX_PARAM_SCRIND, CPX_OFF);
    char cplex_log_file[100];
    sprintf(cplex_log_file, "%s/%lu-%d-%s.log", TSP_CPLEX_LOG_FOLDER, tsp_env.seed, tsp_inst.nnodes, tsp_env.alg_type);
    remove(cplex_log_file);
    *cpxerror = CPXsetlogfilename(*env, cplex_log_file, "w");
    if (*cpxerror) raise_error("CPXsetlogfilename error (%d).\n", *cpxerror);

    // build cplex model
    tsp_cplex_build_model(*env, *lp);

    // give cplex terminate condition
    CPXsetterminate(*env, &tsp_cplex_terminate);

    // create lp file from cplex model
    char cplex_lp_file[100];
    sprintf(cplex_lp_file, "%s/%lu-%d-%s.lp", TSP_CPLEX_LP_FOLDER, tsp_env.seed, tsp_inst.nnodes, tsp_env.alg_type);
    *cpxerror = CPXwriteprob(*env, *lp, cplex_lp_file, NULL);
    if (*cpxerror) raise_error("CPXwriteprob error (%d).\n", *cpxerror);

}

double tsp_cplex_compute_xstar_cost(const double* xstar) {

    double cost = 0;

    for ( int i = 0; i < tsp_inst.nnodes; i++ ) for ( int j = i+1; j < tsp_inst.nnodes; j++ )
        if ( xstar[tsp_cplex_coords_to_xpos(i,j)] > TSP_CPLEX_ZERO_THRESHOLD ) cost += tsp_get_edge_cost(i, j);

    return cost;

}

void tsp_cplex_path_to_ind_val(const int ncols, const int* path, int* indexes, double* values) {

    int k = 0;

    for (int i = 0; i < tsp_inst.nnodes-1; i++) {

        int xpos = tsp_cplex_coords_to_xpos(path[i], path[i+1]);
        indexes[k] = xpos;
        values[k++] = 1.0;

    }

    indexes[k] = tsp_cplex_coords_to_xpos(path[tsp_inst.nnodes - 1], path[0]);
    indexes[k++] = 1;

    // Integrity check
    if (k != tsp_inst.nnodes) raise_error("Something went wrong inside tsp_cplex_path_to_ind_val (k: %d).\n", k);

}

void tsp_cplex_check_best_sol(const int ncomp, const int* comp, const int* succ, const double cost) {

    double time = time_elapsed();
    
    if (ncomp == 1) {

        int* solution = (int*) malloc(tsp_inst.nnodes * sizeof(int));
        tsp_succ_to_path(succ, solution);

        tsp_check_best_sol(solution, cost, time);

        safe_free(solution);

    } else {

        int check = 0;
        for (int k = 1; k <= ncomp; k++)
            for (int i = 0; i < tsp_inst.nnodes; i++) {
                if(comp[i] != k) continue;
                
                int start = i;
                int from = start;
                int to = succ[from];
                while ( to != start ) {
                    from = to;
                    to = succ[to];
                    check++;
                }
                check++;

                break;

            }
            
        // Integrity check
        if (tsp_verbose >= 100) if (check != tsp_inst.nnodes) raise_error("Error while computing the cost of the multitour solution (%d).\n", check);

        if (cost >= tsp_inst.mt_cost - TSP_EPSILON) return;

        pthread_mutex_lock(&tsp_mutex_update_sol);

        if (cost < tsp_inst.mt_cost) {

            tsp_inst.mt_cost = cost;
            tsp_inst.best_time = time;

            tsp_inst.ncomp = ncomp;
            if (tsp_inst.comp == NULL) tsp_inst.comp = (int*)calloc(tsp_inst.nnodes, sizeof(int));
            if (tsp_inst.succ == NULL) tsp_inst.succ = (int*)calloc(tsp_inst.nnodes, sizeof(int));
            
            for (int i = 0; i < tsp_inst.nnodes; i++) { tsp_inst.comp[i] = comp[i]; tsp_inst.succ[i] = succ[i]; }

        }

        pthread_mutex_unlock(&tsp_mutex_update_sol);

    }

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

void tsp_cplex_add_sec(CPXENVptr env, CPXLPptr lp, const int* ncomp, const int* comp, const int* succ) {

    if ((*ncomp)==1) raise_error("ERROR: add_sec() error");

    int cpxerror = 0;

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
        cpxerror = CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]);
        if (cpxerror) raise_error("ERROR: wrong CPXaddrows [degree] (%d).\n", cpxerror);

        safe_free(value);
        safe_free(index);
        safe_free(comp_nodes);

    }

	if (cname[0] != NULL) free(cname[0]);
    safe_free(cname);

}

void tsp_cplex_patching(const double* xstar, int* ncomp, int* comp, int* succ, double* cost) {

    if (cost == NULL) cost = 0;

    switch (tsp_env.cplex_patching) {
        case 1:
            tsp_cplex_patch(ncomp, comp, succ, cost);
            break;
        case 2:
            tsp_cplex_patch_greedy(xstar, ncomp, comp, succ, cost);
            break;
        default:
            raise_error("Error chosing patching function.\n");
    }

    tsp_cplex_check_best_sol(*ncomp, comp, succ, *cost);

}

int tsp_cplex_callback_candidate(CPXCALLBACKCONTEXTptr context, const int nnodes) {

    int cpxerror = 0;

    // get candidate point
    int ncols = nnodes * (nnodes - 1) / 2;
    double* xstar = (double*) malloc(ncols * sizeof(double));
    double objval = CPX_INFBOUND;
    cpxerror = CPXcallbackgetcandidatepoint(context, xstar, 0, ncols-1, &objval);
    if (cpxerror) raise_error("CPXcallbackgetcandidatepoint() error (%d).\n", cpxerror);
    
    if (objval == CPX_INFBOUND) raise_error("CPXcallbackgetcandidatepoint() error, no candidate objval returned.\n");

    // space for data structures
    int* succ = (int*) calloc(tsp_inst.nnodes, sizeof(int));
    int* comp = (int*) calloc(tsp_inst.nnodes, sizeof(int));
    int ncomp = 0;

    // convert xstar to comp and succ
    tsp_cplex_decompose_xstar(xstar, comp, succ, &ncomp);

    if (ncomp == 1) {   // feasible solution (only one tour)
    
        safe_free(comp);
        safe_free(succ);
        safe_free(xstar);

        // user info
        double lower_bound = CPX_INFBOUND; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_BND, &lower_bound);
        double incumbent = CPX_INFBOUND; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent);
        double gap = (1 - lower_bound/incumbent) * 100;

        if (tsp_verbose >= 100) print_info("found feasible solution   -   lower_bound: %15.4f   -   incumbent: %15.4f   -   gap: %6.2f%c.\n", lower_bound, incumbent, gap, '%');

        return 0;

    }

    if (tsp_verbose >= 100) print_info("adding SEC for candidate     -   number of SEC: %d.\n", ncomp);

    // add as many SEC as connected components
    const char sense = 'L';
    const int izero = 0;
    int* index = (int*) calloc(ncols, sizeof(int));
    double* value = (double*) calloc(ncols, sizeof(double));

    for(int k = 1; k <= ncomp; k++) {

        int nnz = 0;
        double rhs = -1.0;

        for(int i = 0; i < tsp_inst.nnodes; i++) {

            if(comp[i] != k) continue;
            rhs++;
            for(int j = i+1; j < tsp_inst.nnodes; j++){
                if(comp[j] != k) continue;
                index[nnz]=tsp_cplex_coords_to_xpos(i,j);
                value[nnz]=1.0;
                nnz++;
            }

        }

        // reject candidate and add SEC
        cpxerror = CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, value);
        if (cpxerror) raise_error("CPXcallbackrejectcandidate() error (%d).\n", cpxerror);

    }

    //TODO: Try patching with some percentage
    if (ncomp != 1 && tsp_env.cplex_patching) {
        tsp_cplex_patching(xstar, &ncomp, comp, succ, NULL);
        //TODO: Give to cplex the patched solution
    }

    // free the memory
    safe_free(value);
    safe_free(index);
    safe_free(comp);
    safe_free(succ);
    safe_free(xstar);

    return 0;

}

int tsp_cplex_callback_relaxation(CPXCALLBACKCONTEXTptr context, const int nnodes) {

    int nodeuid = -1; 
    int cpxerror = 0;
    cpxerror = CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODEUID, &nodeuid);
    if (cpxerror) raise_error("CPXcallbackgetinfoint() error (%d).\n", cpxerror);

    if (nodeuid % 10) return 0; //TODO: Perfprof for different percentages

    if (tsp_verbose >= 1000) print_info("nodeuid: %d.\n", nodeuid);

    int ncols = (nnodes * (nnodes - 1) / 2);
	double* xstar = (double*) malloc(ncols * sizeof(double));  
	double objval = CPX_INFBOUND; 
    cpxerror = CPXcallbackgetrelaxationpoint(context, xstar, 0, ncols-1, &objval);
    if (cpxerror) raise_error("CPXcallbackgetcandidatepoint() error (%d).\n", cpxerror);

    if (objval == CPX_INFBOUND) raise_error("CPXcallbackgetcandidatepoint() error, no candidate objval returned.\n");

    int* elist = (int*) calloc(2 * ncols, sizeof(int));
    int k = 0;
    
    for (int i = 0; i < nnodes; i++) for (int j = i+1; j < nnodes; j++) {
        elist[k++] = i;
        elist[k++] = j;
    }

    int ncomp = -1;
    int* comps = NULL;
    int* compscount = NULL;
    cpxerror = CCcut_connect_components(nnodes, ncols, elist, xstar, &ncomp, &compscount, &comps);
    if (cpxerror) raise_error("CCcut_connect_components() error (%d).\n", cpxerror);

    if (ncomp == 1) {
        
        if (tsp_verbose >= 100) print_info("Adding SEC for relaxation    -   number of SEC: %d.\n", 1);
        int ccerror = CCcut_violated_cuts(nnodes, ncols, elist, xstar, 1.9, tsp_concorde_callback_add_cplex_sec, (void*)&context);
        if (ccerror) raise_error("CCcut_violated_cuts() error (%d).\n", ccerror);
    
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

            safe_free(subtour);

        }

        //TODO: Percentage
        if (tsp_env.cplex_patching) {
            
            int* succ = (int*)malloc(tsp_inst.nnodes * sizeof(int));
            int* comp = (int*)malloc(tsp_inst.nnodes * sizeof(int));
            tsp_cplex_decompose_xstar(xstar, comp, succ, &ncomp);

            tsp_cplex_patching(xstar, &ncomp, comp, succ, NULL);

            //TODO: Give to cplex the patched solution

            safe_free(comp);
            safe_free(succ);

        }

    }

    safe_free(compscount);
    safe_free(comps);

    safe_free(elist);
    safe_free(xstar);

    return 0;
    
}

void tsp_cplex_close(CPXENVptr env, CPXLPptr lp, double* xstar, int* comp, int* succ) {

    // free memory
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
    safe_free(comp);
    safe_free(succ);
    safe_free(xstar);

    // remove "clone<x>.log" files generated by cplex
    int file_number = 1;
    char clone_file[50];
    sprintf(clone_file, "clone0.log");
    remove(clone_file);
    do { sprintf(clone_file, "clone%d.log", file_number++); } while (!remove(clone_file));

}