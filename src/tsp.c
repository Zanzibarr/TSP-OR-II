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

    safe_free(tsp_inst.coords);
    safe_free(tsp_inst.costs);
    safe_free(tsp_inst.sort_edges);
    safe_free(tsp_inst.best_solution);

    for (int thread = 0; thread < N_THREADS; thread++) safe_free(tsp_env.tabu_tables[thread].list);

    pthread_mutex_destroy(&tsp_mutex_update_sol);

    safe_free(tsp_inst.comp);
    safe_free(tsp_inst.succ);

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

    safe_free(list);

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

void tsp_check_best_sol(const int* path, const double cost, const double time) {

    // Integrity check
    if (tsp_verbose >= 100) tsp_check_integrity(path, cost, "tsp.c - tsp_check_best_sol - 1");

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

double tsp_get_edge_cost(const int i, const int j) {

    return tsp_inst.costs[i * tsp_inst.nnodes + j];

}

void tsp_succ_to_path(const int* succ, int* path) {

    int* visited = (int*)calloc(tsp_inst.nnodes, sizeof(int));

    safe_free(visited);

    for (int i=0, current_node=0; i<tsp_inst.nnodes; i++, current_node=succ[current_node]) path[i] = current_node;

    // Integrity check
    if (tsp_verbose >= 100) tsp_check_integrity(path, tsp_compute_path_cost(path), "tsp.c - tsp_succ_to_path.\n");

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

void tsp_2opt(int* path, double* cost, const int (*swap_function)(int*, double*)) {

    while ((*swap_function)(path, cost) > 0 && time_elapsed() < tsp_env.time_limit); {//repeat until I can't find new swaps that improve the cost of my solution
        // Integrity check
        if (tsp_verbose >= 100) tsp_check_integrity(path, *cost, "tsp.c - tsp_2opt - 1");
    }

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

}

/**
 * @brief Initializes the problem instance
*/
void tsp_init_inst() {

    tsp_inst.nnodes = TSP_DEF_NNODES;
    tsp_inst.best_cost = INFINITY;
    tsp_inst.best_time = 0;
    tsp_inst.ncomp = 1;
    tsp_inst.mt_cost = INFINITY;

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

int tsp_save_solution() {
    
    FILE *solution_file;

    char prefix[150], solution_file_name[500];

    struct timeval tv; gettimeofday(&tv, NULL);
    int timestamp = (int)tv.tv_sec;

    if (tsp_env.seed > 0) 
        snprintf(prefix, sizeof(char)*150, "%lu_%d_%s_%d", tsp_env.seed, tsp_inst.nnodes, tsp_env.alg_type, timestamp);
    else
        snprintf(prefix, sizeof(char)*150, "%s_%s_%d", tsp_env.file_name, tsp_env.alg_type, timestamp);
    snprintf(solution_file_name, sizeof(char)*500, "%s/%s_%s", TSP_SOL_FOLDER, prefix, TSP_SOLUTION_FILE);  //where to save the file

    solution_file = fopen(solution_file_name, "w");

    if (solution_file == NULL) raise_error("Error writing the file for the solution.");

    fprintf(solution_file, "Algorithm: %s\n", tsp_env.alg_type);

    if (tsp_env.g2opt_swap_pol) fprintf(solution_file, "Swap policy: %s.\n", ((tsp_env.g2opt_swap_pol == 1) ? "first swap" : "best swap"));
    if (!strcmp(tsp_env.alg_type, TSP_PARSING_TABU)) fprintf(solution_file, "Tabu tenure: %4d.\nTabu variability: %4d.\nTabu variability frequency: %10.4f.\n", tsp_env.tabu_tenure, tsp_env.tabu_tenure_a, tsp_env.tabu_tenure_f);
    if (tsp_env.vns_fvns) fprintf(solution_file, "Fast vns enabled.\n");
    if (tsp_env.cplex_mipstart) fprintf(solution_file, "Using a mipstart.\n");
    if (tsp_env.cplex_benders) fprintf(solution_file, "Using bender's loop.\n");
    switch (tsp_env.cplex_patching) {
        case 1:
            fprintf(solution_file, "Using normal patching.\n");
            break;
        case 2:
            fprintf(solution_file, "Using greedy patching.\n");
            break;
        case 0: break;
        default:
            raise_error("Error choosing the patching function.\n");
    }
    if (tsp_env.cplex_can_cb) fprintf(solution_file, "Using candidate callback.\n");
    if (tsp_env.cplex_rel_cb) fprintf(solution_file, "Using relaxation callback.\n");
    if (tsp_env.tmp_choice) fprintf(solution_file, "Added temporary option.\n");

    fprintf(solution_file, "Cost: %15.4f\n", (tsp_inst.ncomp == 1) ? tsp_inst.best_cost : tsp_inst.mt_cost);
    fprintf(solution_file, "Time: %15.4fs\n", tsp_inst.best_time);
    fprintf(solution_file, "Total execution time: %15.4fs\n", tsp_env.time_total);

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
        default: raise_error("Unexpected status.\n");
    }
    
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

    return timestamp;

}

void tsp_plot_solution(const int unique) {

    //FIXME: Move this to python

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
    switch (tsp_env.cplex_patching) {
        case 1:
            printf("Using normal patching.\n");
            break;
        case 2:
            printf("Using greedy patching.\n");
            break;
        case 0: break;
        default:
            raise_error("Error choosing the patching function.\n");
    }
    if (tsp_env.cplex_can_cb) printf("Using candidate callback.\n");
    if (tsp_env.cplex_rel_cb) printf("Using relaxation callback.\n");
    if (tsp_env.tmp_choice) printf("Added temporary option.\n");
    
    printf("--------------------\n");

    if (tsp_verbose >= 100) {
        printf("Integrity checks enabled (all kind).\n");
        print_warn("Execution might be slowed down.\n");
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
    printf("Cost: %15.4f\n", (tsp_inst.ncomp == 1) ? tsp_inst.best_cost : tsp_inst.mt_cost);
    printf("Time:\t%15.4fs\n", tsp_inst.best_time);
    printf("Execution time: %8.4fs\n", tsp_env.time_total);
    if (tsp_verbose >= 500) {
        //FIXME: Check for ncomp
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
        case 0: break;
        default: raise_error("Unexpected status.\n");
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
        if (error == 1) raise_error("Non-existent node in path.\n");
        else if (error == 2) raise_error("Double node in path.\n");
        else if (error == 3) raise_error("Cost: %.10f, Checked cost: %.10f, Difference: %.10f, Threshold: %.10f\n", cost, c_cost, fabs(c_cost - cost), TSP_EPSILON);
        else raise_error("Unknown error.\n");
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