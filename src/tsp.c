#include "../include/tsp.h"


#pragma region GLOBALS DEFINITIONS

tsp_instance tsp_inst;

double tsp_initial_time;
double tsp_total_time;
double tsp_time_limit;

int tsp_over_time;
int tsp_forced_termination;

char tsp_algorithms[TSP_ALG_NUMBER][50] = {"greedy", "g2opt", "g2opt-best", "tabu", "vns", "fvns", "cplex-base", "benders"};

uint64_t tsp_seed;

char tsp_alg_type[20];
char tsp_file_name[100];

int tsp_tabu_tenure;
int tsp_tabu_tenure_a;
double tsp_tabu_tenure_f;

char tsp_intermediate_costs_files[N_THREADS][30];

tsp_tabu tsp_tabu_tables[N_THREADS];

CPXENVptr tsp_cplex_env;
CPXLPptr tsp_cplex_lp;
double* tsp_cplex_solution;
double tsp_cplex_solution_cost;

/*
char tsp_test_flag = 0;
char tsp_test_run_file[100];
*/
#pragma endregion


#pragma region PRIVATE_METHODS
/**
 * @brief Compute the euclidian distance between the points i and j
 * 
 * @param i First node
 * @param j Second node
 * 
 * @return The euclidian distance between the two points
*/
double tsp_compute_distance(const int i, const int j) {
    
    return sqrt(pow(tsp_inst.coords[i].x - tsp_inst.coords[j].x, 2) + pow(tsp_inst.coords[i].y - tsp_inst.coords[j].y, 2));
    
}

/**
 * @brief Computes the dinamic tenure based on the counter
 * 
 * @param counter Counter used to calculate the tenure
 * @return The new tenure
 */
double tsp_dinamic_tenue(int counter) {

    return tsp_tabu_tenure_a * sin((double)counter * tsp_tabu_tenure_f) + tsp_tabu_tenure;

}

/**
 * @brief Debugging tool used to check whether the sorting of the sort_edges list is correct
*/
void tsp_check_sort_edges_integrity() {
    
    for (int i = 0; i < tsp_inst.nnodes; i++) {
        double min_cost = 0;
        for (int j = 0; j < tsp_inst.nnodes - 1; j++) {
            double checked_cost = tsp_inst.costs[i * tsp_inst.nnodes + tsp_inst.sort_edges[i * (tsp_inst.nnodes - 1) + j]];

            if (checked_cost < min_cost - TSP_EPSILON) {
                printf("SORT_EDGES INTEGRITY COMPROMISED\n");
                exit(1);
            }
            
            min_cost = checked_cost;
        }
    }

    #if TSP_VERBOSE >= 1000
    printf("sort_edges integrity check passed.\n");
    #endif

}
#pragma endregion


#pragma region MEMORY MANAGEMENT
/**
 * @brief Dinamically allocate the space for the costs list
*/
void tsp_allocate_costs_space() {

    if (!tsp_inst.nnodes) { printf("The nnodes variable hasn't been assigned yet."); exit(1); }

    tsp_inst.costs = (double*)calloc(tsp_inst.nnodes * tsp_inst.nnodes, sizeof(double));

}

/**
 * @brief Dinamically allocate the space for the sort_edges list
*/
void tsp_allocate_sort_edges_space() {

    if (!tsp_inst.nnodes) { printf("The nnodes variable hasn't been assigned yet."); exit(1); }

    tsp_inst.sort_edges = (int*)calloc(tsp_inst.nnodes * (tsp_inst.nnodes - 1), sizeof(int));

}

/**
 * @brief Dinamically allocate the space for the best solution
*/
void tsp_allocate_best_sol_space() {

    if (!tsp_inst.nnodes) { printf("The nnodes variable hasn't been assigned yet."); exit(1); }

    tsp_inst.best_solution = (int*)calloc(tsp_inst.nnodes, sizeof(int));
    
}

/**
 * @brief Dinamically allocate the space for the tabu list
*/
void tsp_allocate_tabu_space() {

    for (int thread = 0; thread < N_THREADS; thread++) {
        tsp_tabu_tables[thread].list = (tsp_tabu_entry*)calloc(tsp_inst.nnodes, sizeof(tsp_tabu_entry));
        for (int i = 0; i < tsp_inst.nnodes; i++) { tsp_tabu_tables[thread].list[i] = (tsp_tabu_entry){-1, -1, -1, -1}; }
    }

}

void tsp_allocate_coords_space() {

    if (!tsp_inst.nnodes) { printf("The nnodes variable hasn't been assigned yet."); exit(1); }

    tsp_inst.coords = (tsp_pair*)calloc(tsp_inst.nnodes, sizeof(tsp_pair));

}

void tsp_free_instance() {

    if (tsp_inst.coords != NULL) { free(tsp_inst.coords); tsp_inst.coords = NULL; }
    if (tsp_inst.costs != NULL) { free(tsp_inst.costs); tsp_inst.costs = NULL; }
    if (tsp_inst.sort_edges != NULL) { free(tsp_inst.sort_edges); tsp_inst.sort_edges = NULL; }
    if (tsp_inst.best_solution != NULL) { free(tsp_inst.best_solution); tsp_inst.best_solution = NULL; }

    for (int thread = 0; thread < N_THREADS; thread++)
        if (tsp_tabu_tables[thread].list != NULL) { free(tsp_tabu_tables[thread].list); tsp_tabu_tables[thread].list = NULL; }

    pthread_mutex_destroy(&tsp_mutex_update_sol);

}
#pragma endregion


#pragma region PRECOMPUTING
void tsp_precompute_sort_edges() {

    tsp_allocate_sort_edges_space();
    tsp_entry* list = (tsp_entry*)calloc(tsp_inst.nnodes, sizeof(tsp_entry));

    for (int i = 0; i < tsp_inst.nnodes; i++) {

        for (int j = 0; j < tsp_inst.nnodes; j++) {  //saving the entries to be sorted
            list[j].key = j;
            list[j].value = tsp_inst.costs[i * tsp_inst.nnodes + j];    //considering only the costs of the edges leaving node i
        }

        qsort((void*)list, (size_t)tsp_inst.nnodes, sizeof(tsp_entry), compare_tsp_entries); //sort by cost of the edge
        
        for (int j = 1; j < tsp_inst.nnodes; j++)
            tsp_inst.sort_edges[i * (tsp_inst.nnodes - 1) + j-1] = list[j].key;    //populate the ith row with the nodes ordered by increasing distance

    }

    if (list != NULL) { free(list); list = NULL; }

    #if TSP_VERBOSE >= 100
    tsp_check_sort_edges_integrity();
    #endif

}

void tsp_precompute_costs() {

    tsp_allocate_costs_space();

    for (int i = 0; i < tsp_inst.nnodes; i++) for (int j = 0; j < tsp_inst.nnodes; j++)
        tsp_inst.costs[i * tsp_inst.nnodes + j] = tsp_compute_distance(i, j);

}
#pragma endregion


#pragma region ALGORITHMS TOOLS
int compare_tsp_entries(const void* arg1, const void* arg2) {

    double diff = ((tsp_entry*)arg1)->value - ((tsp_entry*)arg2)->value;
    return (fabs(diff) < TSP_EPSILON) ? 0 : ((diff < 0) ? -1 : 1);

}

void tsp_check_best_sol(const int* path, const double cost, const double time) {

    if (cost > tsp_inst.best_cost + TSP_EPSILON) return;

    pthread_mutex_lock(&tsp_mutex_update_sol);

    if (cost < tsp_inst.best_cost - TSP_EPSILON) {

        for (int i = 0; i < tsp_inst.nnodes; i++) tsp_inst.best_solution[i] = path[i];
        tsp_inst.best_cost = cost;
        tsp_inst.best_time = time;

        #if TSP_VERBOSE >= 10
        printf("Time: %10.4f  -  New best solution : %15.4f\n", tsp_time_elapsed(), cost);
        #endif

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

int tsp_check_tabu(int t_index, int from, int to) {

    if (from > to) return tsp_check_tabu(t_index, to, from);

    return 
        tsp_tabu_tables[t_index].list[from].counter_1 != -1 &&  //if I have a tabu saved from the "from" node
        (
            tsp_tabu_tables[t_index].list[from].node_1 == to && //if the node "to" creates a tabu "from"-"to" (first option)...
            tsp_tabu_tables[t_index].counter - tsp_tabu_tables[t_index].list[from].counter_1 < tsp_dinamic_tenue(tsp_tabu_tables[t_index].counter)   //... and I still remember that as a tabu 
            ||
            tsp_tabu_tables[t_index].list[from].node_2 == to && //if the node "to" creates a tabu "from"-"to" (second option)...
            tsp_tabu_tables[t_index].counter - tsp_tabu_tables[t_index].list[from].counter_2 < tsp_dinamic_tenue(tsp_tabu_tables[t_index].counter)   //... and I still remember that as a tabu 
        );

}

void tsp_add_tabu(int t_index, int from, int to) {

    if (from > to) return tsp_add_tabu(t_index, to, from);

    if (tsp_tabu_tables[t_index].list[from].counter_1 != -1 && tsp_tabu_tables[t_index].counter - tsp_tabu_tables[t_index].list[from].counter_1 < tsp_dinamic_tenue(tsp_tabu_tables[t_index].counter)) {

        tsp_tabu_tables[t_index].list[from].node_2 = to;
        tsp_tabu_tables[t_index].list[from].counter_2 = tsp_tabu_tables[t_index].counter++;

    } else {

        tsp_tabu_tables[t_index].list[from].node_1 = to;
        tsp_tabu_tables[t_index].list[from].counter_1 = tsp_tabu_tables[t_index].counter++;

    }

}

double tsp_compute_path_cost(const int* path) {

    double cost = 0;

    for (int i = 0; i < tsp_inst.nnodes - 1; i++) cost += tsp_inst.costs[path[i] * tsp_inst.nnodes + path[i+1]];
    cost += tsp_inst.costs[path[tsp_inst.nnodes - 1] * tsp_inst.nnodes + path[0]];

    return cost;

}
#pragma endregion


#pragma region INITIALIZATIONS
void tsp_init_defs() {

    struct timeval tv;
    gettimeofday(&tv, NULL);
    tsp_initial_time = ((double)tv.tv_sec)+((double)tv.tv_usec/1e+6);

    tsp_seed = 0;
    tsp_time_limit = (double)TSP_DEF_TL/1000;
    
    tsp_over_time = 0;
    tsp_forced_termination = 0;

    tsp_tabu_tenure = TSP_DEF_TABU_TENURE;
    tsp_tabu_tenure_a = 0;
    tsp_tabu_tenure_f = .02;

    strcpy(tsp_file_name, "NONE");
    strcpy(tsp_alg_type, "greedy");

    tsp_inst.nnodes = TSP_DEF_NNODES;

    for (int i = 0; i < N_THREADS; i++) {
        snprintf(tsp_intermediate_costs_files[i], 30*sizeof(char), "%s/int_costs_%d.txt", TSP_PLOT_FOLDER, i);
        FILE* file = fopen(tsp_intermediate_costs_files[i], "w");
        fclose(file);
    }

    tsp_init_threads();

}

void tsp_init_solution() {

    tsp_allocate_best_sol_space();
    tsp_inst.best_cost = INFINITY;
    tsp_inst.best_time = 0;

    tsp_over_time = 0;
    tsp_forced_termination = 0;

    tsp_allocate_tabu_space();

}
#pragma endregion


#pragma region SAVING FILES
void tsp_save_solution() {
    
    FILE *solution_file;

    char prefix[150], solution_file_name[500];

    if (tsp_seed > 0) 
        snprintf(prefix, sizeof(char)*150, "%llu_%d_%s", tsp_seed, tsp_inst.nnodes, tsp_alg_type);
    else
        snprintf(prefix, sizeof(char)*150, "%s_%s", tsp_file_name, tsp_alg_type);
    snprintf(solution_file_name, sizeof(char)*500, "%s/%s_%s", TSP_SOL_FOLDER, prefix, TSP_SOLUTION_FILE);  //where to save the file

    solution_file = fopen(solution_file_name, "w");

    if (solution_file == NULL) {
        printf("Error writing the file for the solution.");
        exit(1);
    }

    fprintf(solution_file, "Algorithm: %s\n", tsp_alg_type);
    fprintf(solution_file, "Cost: %15.4f\n", tsp_inst.best_cost);
    fprintf(solution_file, "Time: %15.4fs\n", tsp_inst.best_time);
    fprintf(solution_file, "Total execution time: %15.4fs\n", tsp_total_time);
    if (tsp_over_time) fprintf(solution_file, "The algorithm has exceeded the time limit and has been stopped.\n");
    else if (tsp_forced_termination) fprintf(solution_file, "The algorithm has been terminated by the user.\n");
    else fprintf(solution_file, "The algorithm hasn't exceeded the time limit.\n");
    fprintf(solution_file, "--------------------\n");

    if (!strncmp(tsp_alg_type, "cplex", 5)) {

        int ncomp = 0;
        char* read_nodes_flags = (char*) calloc(tsp_inst.nnodes, 1);
        int start_node, current_node, read_nodes=0;
        while (read_nodes<tsp_inst.nnodes) {
            fprintf(solution_file, "LOOP %d:\n", ++ncomp);
            start_node = -1;
            while (read_nodes_flags[++start_node] && start_node<tsp_inst.nnodes);
            current_node = start_node;
            do {
                fprintf(solution_file, "%4d %15.4f %15.4f\n", current_node, tsp_inst.coords[current_node].x, tsp_inst.coords[current_node].y);
                read_nodes_flags[current_node] = 1;
                read_nodes++;
                current_node = tsp_inst.best_solution[current_node];
            } while (current_node!=start_node);
            fprintf(solution_file, "%4d %15.4f %15.4f\n", start_node, tsp_inst.coords[start_node].x, tsp_inst.coords[start_node].y);
        }

    }
    else {

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

    if (tsp_seed > 0) 
        snprintf(prefix, sizeof(char)*150, "%llu_%d_%s", tsp_seed, tsp_inst.nnodes, tsp_alg_type);
    else
        snprintf(prefix, sizeof(char)*150, "%s_%s", tsp_file_name, tsp_alg_type);

    snprintf(plot_file_name, sizeof(char)*500, "%s/%s_%s", TSP_SOL_FOLDER, prefix, TSP_PLOT_FILE);  //where to save the plot
    snprintf(solution_file_name, sizeof(char)*500, "%s/%s_%s", TSP_SOL_FOLDER, prefix, TSP_SOLUTION_FILE);  //where to read the file from

    solution_file = fopen(solution_file_name, "r");
    command_file = fopen(TSP_COMMAND_FILE, "w");

    if (solution_file == NULL || command_file == NULL) {
        printf("Error with a file used to plot the solution.");
        exit(1);
    }

    // skip through the rows with the solution info
    while (rows_read < 6) if (fgetc(solution_file) =='\n') rows_read++;

    // build plot title with solution info
    snprintf(gnuplot_title, 1000, "Algorithm: %s; number of nodes: %d; cost: %.4f; time: %.4fs", tsp_alg_type, tsp_inst.nnodes, tsp_inst.best_cost, tsp_inst.best_time);

    // copy nodes coordinates into coords_file
    if (strncmp(tsp_alg_type, "cplex", 5)) {
        FILE *coords_file = fopen(TSP_COORDS_FILE, "w");
        while (fgets(solution_contents, 100, solution_file)) fprintf(coords_file, "%s", solution_contents);
        fclose(coords_file);
    }
    else {
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
    
    // builds commands for gnuplot
    fprintf(command_file, "set term png\n");
    fprintf(command_file, "set output '%s'\n", plot_file_name);
    fprintf(command_file, "x=0.; y=0.\n");
    fprintf(command_file, "set title '%s'\n", gnuplot_title);
    fprintf(command_file, "set xlabel 'Starting node highlighted as black point'\n");
    if (strncmp(tsp_alg_type, "cplex", 5)) {
        fprintf(command_file, "set label at %f, %f point pointtype 7 pointsize 2\n", tsp_inst.coords[tsp_inst.best_solution[0]].x, tsp_inst.coords[tsp_inst.best_solution[0]].y);
        fprintf(command_file, "plot '%s' u (x=$2):(y=$3) w lp lc rgb 'blue' title ''\n", TSP_COORDS_FILE);
    }
    else {
        fprintf(command_file, "plot ");
        for (int i=0; i<coord_files_number; i++) {
            fprintf(command_file, "'%d_%s' u (x=$2):(y=$3) w lp lc rgb 'blue' title ''", i+1, TSP_COORDS_FILE);
            if (i!=coord_files_number-1) fprintf(command_file, ", ");
        }
        fprintf(command_file, "\n");
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

void tsp_save_intermediate_cost(const int list_index, const double cost) {
    
    FILE* file = fopen(tsp_intermediate_costs_files[list_index], "a");

    fprintf(file, "%f\n", cost);

    fclose(file);

}
#pragma endregion


#pragma region CPLEX

int tsp_cplex_coords_to_xpos(const int i, const int j) {

	if ( i == j ) { printf("ERROR: i == j in xpos"); exit(1); }
	if ( i > j ) return tsp_cplex_coords_to_xpos(j,i);

	return i * tsp_inst.nnodes + j - (( i + 1 ) * ( i + 2 )) / 2;

}

void tsp_cplex_build_model() {

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
			if ( CPXnewcols(tsp_cplex_env, tsp_cplex_lp, 1, &obj, &lb, &ub, &binary, cname) ) { printf("ERROR: wrong CPXnewcols on x var.s"); exit(1); }
    		if ( CPXgetnumcols(tsp_cplex_env, tsp_cplex_lp)-1 != tsp_cplex_coords_to_xpos(i,j) ) { printf("ERROR: wrong position for x var.s"); exit(1); }

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
		
		if ( CPXaddrows(tsp_cplex_env, tsp_cplex_lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) ) { printf("ERROR: wrong CPXaddrows [degree]"); exit(1); }

	} 

    if (value != NULL) { free(value); value = NULL; }
    if (index != NULL) { free(index); index = NULL; }	

	#if TSP_VERBOSE >= 100
    CPXwriteprob(tsp_cplex_env, tsp_cplex_lp, "model.tsp_cplex_lp", NULL);
    #endif

	if (cname[0] != NULL) { free(cname[0]); cname[0] = NULL; }
	if (cname != NULL) { free(cname); cname = NULL; }

}

void tsp_cplex_build_solution(int *ncomp, int *comp, int *succ) {

    #if TSP_VERBOSE >= 100

        int *degree = (int *) calloc(tsp_inst.nnodes, sizeof(int));
        for ( int i = 0; i < tsp_inst.nnodes; i++ ) {
            for ( int j = i+1; j < tsp_inst.nnodes; j++ ) {
                int k = tsp_cplex_coords_to_xpos(i,j);
                if ( fabs(tsp_cplex_solution[k]) > TSP_EPSILON && fabs(tsp_cplex_solution[k]-1.0) > TSP_EPSILON )
                    { printf(" wrong tsp_cplex_solution in build_sol()"); exit(1); }
                if ( tsp_cplex_solution[k] > 0.5 ) {
                    ++degree[i];
                    ++degree[j];
                }
            }
        }
        for ( int i = 0; i < tsp_inst.nnodes; i++ ) {
            if ( degree[i] != 2 ) { printf("wrong degree in build_sol()"); exit(1); }
        }	
        free(degree);

    #endif

	*ncomp = 0;
	for ( int i = 0; i < tsp_inst.nnodes; i++ ) {
		succ[i] = -1;
		comp[i] = -1;
	}
	
	for ( int start = 0; start < tsp_inst.nnodes; start++ ) {
		
        if ( comp[start] >= 0 ) continue;  // node "start" was already visited, just skip it

		// a new component is found
		(*ncomp)++;
		int i = start;
		int done = 0;
		while ( !done ) { // go and visit the current component
			comp[i] = *ncomp;
			done = 1;
			for ( int j = 0; j < tsp_inst.nnodes; j++ ) {
				if ( i != j && tsp_cplex_solution[tsp_cplex_coords_to_xpos(i,j)] > TSP_CPLEX_ZERO_THRESHOLD && comp[j] == -1 ) {
                    // the edge [i,j] is selected in tsp_cplex_solution and j was not visited before 
					succ[i] = j;
					i = j;
					done = 0;
					break;
				}
			}
		}	
		succ[i] = start;  // last arc to close the cycle
		
		// go to the next component...

	}

}

void tsp_cplex_save_solution() {
	
    //CPXgetx get the solution x* of our instance
	if ( CPXgetx(tsp_cplex_env, tsp_cplex_lp, tsp_cplex_solution, 0, CPXgetnumcols(tsp_cplex_env, tsp_cplex_lp)-1) )
        { printf("ERROR: CPXgetx() error"); exit(1); }	

    #if TSP_VERBOSE>0

	for ( int i = 0; i < tsp_inst.nnodes; i++ ) {

		for ( int j = i+1; j < tsp_inst.nnodes; j++ ) {

			if ( tsp_cplex_solution[tsp_cplex_coords_to_xpos(i,j)] > TSP_CPLEX_ZERO_THRESHOLD )
                printf("x(%3d,%3d) = 1\n", i,j);

		}

	}

    #endif

}

void tsp_cplex_convert_solution() {

    tsp_inst.best_cost = tsp_cplex_solution_cost;

    for (int i=0; i<tsp_inst.nnodes; i++) tsp_inst.best_solution[i] = -1;
    for (int i=0; i<tsp_inst.nnodes; i++) {
        if (tsp_inst.best_solution[i]!=-1) continue;
        int start_node = i, current_node = i;
        do {
            for (int j=0; j<tsp_inst.nnodes && tsp_inst.best_solution[current_node]==-1; j++) {
                if (j==current_node || tsp_inst.best_solution[j]==current_node) continue;
                if (tsp_cplex_solution[tsp_cplex_coords_to_xpos(current_node, j)] > TSP_CPLEX_ZERO_THRESHOLD) {
                    tsp_inst.best_solution[current_node] = j;
                    current_node = j;
                }
            }
        } while (current_node!=start_node);
    }
    // TODO: scrivere tempo in soluzione

}

void tsp_cplex_add_sec(int *ncomp, int* comp) {

    if (*ncomp==1) { printf("ERROR: add_sec() error"); exit(1); }

    int* index = (int*) calloc(CPXgetnumcols(tsp_cplex_env, tsp_cplex_lp), sizeof(double));
    double* value = (double*) calloc(CPXgetnumcols(tsp_cplex_env, tsp_cplex_lp), sizeof(double));
    int izero = 0;
    char** cname = (char**)calloc(1, sizeof(char *));	// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));

    for (int k=0; k<*ncomp; k++) {
        int nnz = 0;
        char sense = 'L';
        double rhs = -1.0;
        for (int i=0; i<tsp_inst.nnodes; i++) {
            if (comp[i]!=k+1) continue;
            rhs += 1.0;
            for (int j=0; j<tsp_inst.nnodes; j++) {
                if (comp[j]!=k+1) continue;
                index[nnz] = tsp_cplex_coords_to_xpos(i, j);
                value[nnz] = 1.0;
                nnz++;
                if ( CPXaddrows(tsp_cplex_env, tsp_cplex_lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) )
                    { printf("ERROR: wrong CPXaddrows [degree]"); exit(1); }
            }
        }
    }

}

#pragma endregion


#pragma region DEBUGGING TOOLS
void tsp_instance_info() {

    printf("--------------------\n");
    printf("Type of Instance: %s\n", ((tsp_seed == 0) ? "from file" : "random"));
    if (tsp_seed == 0) printf("File name: %s\n", tsp_file_name);
    else printf("Seed: %llu\n", tsp_seed);
    printf("Time limit: %10.4fs\n", tsp_time_limit);
    printf("Number of nodes: %4d\n", tsp_inst.nnodes);
    printf("Edge weight type: ATT\n");
    printf("--------------------\n");
    printf("Algorithm: %s\n", tsp_alg_type);
    if (!strcmp(tsp_alg_type, "tabu")) { printf("Fixed tenure: %d\nTenure variability: %d\nTenure frequency: %f\n", tsp_tabu_tenure, tsp_tabu_tenure_a, tsp_tabu_tenure_f); }
    printf("--------------------\n");

    #if TSP_VERBOSE >= 100
    printf("Integrity checks enabled.\n");
    printf("--------------------\n");
    #endif

    #if TSP_VERBOSE < 1000
    return;
    #endif

    printf("NODES:\n");
    for (int i = 0; i < tsp_inst.nnodes; i++) printf("node[%4d]: (%15.4f, %15.4f)\n", i, tsp_inst.coords[i].x, tsp_inst.coords[i].y);
    printf("--------------------\n");
    printf("COSTS\n");
    for (int i = 0; i < tsp_inst.nnodes; i++) for (int j = 0; j < tsp_inst.nnodes; j++) printf("v%4d - v%4d : %15.4f\n", i, j, tsp_inst.costs[i * tsp_inst.nnodes + j]);
    printf("--------------------\n");
    printf("MIN EDGES\n");
    for (int i = 0; i < tsp_inst.nnodes; i++) {
        for (int j = 0; j < tsp_inst.nnodes - 1; j++)
            printf("%4d(%15.4f) ", tsp_inst.sort_edges[i * (tsp_inst.nnodes - 1) + j], tsp_inst.costs[i * (tsp_inst.nnodes) + tsp_inst.sort_edges[i * (tsp_inst.nnodes - 1) + j]]);
        printf("\n");
    }
    printf("--------------------\n");

}

void tsp_print_solution() {

    printf("--------------------\nBEST SOLUTION:\n");
    printf("Cost: %15.4f\n", tsp_inst.best_cost);
    printf("Time:\t%15.4fs\n", tsp_inst.best_time);
    printf("Execution time %8.4fs\n", tsp_total_time);
    #if TSP_VERBOSE >= 500
    for (int i = 0; i < tsp_inst.nnodes; i++) printf("%d -> ", tsp_inst.best_solution[i]);
    printf("%d\n", tsp_inst.best_solution[0]);
    #endif
    if (tsp_over_time) printf("The algorithm exceeded the time limit and has been stopped.\n");
    else if (tsp_forced_termination) printf("The algorithm has been terminated by the user.\n");
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
        printf("\n\nINTEGRITY COMPROMISED - error_code: %d\n----- %s\n", error, message);
        if (error == 1) printf("\nNon-existent node in path.\n");
        else if (error == 2) printf("\nDouble node in path.\n");
        else if (error == 3) printf("\nCost: %.10f, Checked cost: %.10f, Difference: %.10f, Threshold: %.10f\n", cost, c_cost, fabs(c_cost - cost), TSP_EPSILON);
        else printf("\nUnknown error.\n");
        exit(1);
    }

    #if TSP_VERBOSE >= 1000
    printf("Integrity check passed.\n");
    #endif

}
#pragma endregion


#pragma region USEFUL METHODS
void tsp_init_rand() { for (int i = 0; i < 100; i++) rand(); }

double tsp_rnd_coord() { return (double)rand()/RAND_MAX*TSP_GRID_SIZE; }

double tsp_time_elapsed() {

    struct timeval tv;

    gettimeofday(&tv, NULL);
    return ((double)tv.tv_sec)+((double)tv.tv_usec/1e+6) - tsp_initial_time;

}
#pragma endregion