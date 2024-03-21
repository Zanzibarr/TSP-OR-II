#include "../include/tsp.h"


#pragma region GLOBALS DEFINITIONS

tsp_instance tsp_inst;

double tsp_initial_time;
double tsp_total_time;
double tsp_time_limit;

int tsp_over_time;
int tsp_forced_termination; 

uint64_t tsp_seed;

char tsp_alg_type[20];
char tsp_file_name[100];

int tsp_tabu_tenure;
int tsp_tabu_tenure_a;
double tsp_tabu_tenure_f;

char tsp_intermediate_costs_files[N_THREADS][30];

tsp_tabu tsp_tabu_tables[N_THREADS];

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
 * @brief Comparator used by the qsort method to sort the sort_edges list
 * 
 * @param arg1 The first element to compare (casted as tsp_entry)
 * @param arg2 The second element to compare (casted as tsp_entry)
 * 
 * @return -1 if arg1<arg2, 0 if arg1 == arg2, 1 if arg1 > arg2
*/
int compare_tsp_entries(const void* arg1, const void* arg2) {

    double diff = ((tsp_entry*)arg1)->value - ((tsp_entry*)arg2)->value;
    return (fabs(diff) < TSP_EPSILON) ? 0 : ((diff < 0) ? -1 : 1);

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
 * @brief Debugging tool used to check wether the sorting of the sort_edges list is correct
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

    tsp_allocate_sort_edges_space(tsp_inst);
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
    tsp_check_sort_edges_integrity(tsp_inst);
    #endif

}

void tsp_precompute_costs() {

    tsp_allocate_costs_space(tsp_inst);

    for (int i = 0; i < tsp_inst.nnodes; i++) for (int j = 0; j < tsp_inst.nnodes; j++)
        tsp_inst.costs[i * tsp_inst.nnodes + j] = tsp_compute_distance(i, j);

}
#pragma endregion


#pragma region ALGORITHMS TOOLS
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

    if (from > to) { int c = from; from = to; to = c;}

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

    if (from > to) { int c = from; from = to; to = c;}

    if (tsp_tabu_tables[t_index].list[from].counter_1 != -1 && tsp_tabu_tables[t_index].counter - tsp_tabu_tables[t_index].list[from].counter_1 < tsp_dinamic_tenue(tsp_tabu_tables[t_index].counter)) {

        tsp_tabu_tables[t_index].list[from].node_2 = to;
        tsp_tabu_tables[t_index].list[from].counter_2 = tsp_tabu_tables[t_index].counter++;

    } else {

        tsp_tabu_tables[t_index].list[from].node_1 = to;
        tsp_tabu_tables[t_index].list[from].counter_1 = tsp_tabu_tables[t_index].counter++;

    }

}
#pragma endregion


#pragma region INIZIALIZATIONS
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
        snprintf(prefix, sizeof(char)*150, "%ld_%d_%s", tsp_seed, tsp_inst.nnodes, tsp_alg_type);
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
    for (int i = 0; i < tsp_inst.nnodes; i++)
        fprintf(solution_file, "%4d %15.4f %15.4f\n", tsp_inst.best_solution[i], tsp_inst.coords[tsp_inst.best_solution[i]].x, tsp_inst.coords[tsp_inst.best_solution[i]].y);
    fprintf(solution_file, "%4d %15.4f %15.4f\n", tsp_inst.best_solution[0], tsp_inst.coords[tsp_inst.best_solution[0]].x, tsp_inst.coords[tsp_inst.best_solution[0]].y);

    fclose(solution_file);

}

void tsp_plot_solution() {

    int rows_read = 0;
    FILE *solution_file, *coords_file, *command_file;
    char plot_file_name[500], solution_file_name[500], solution_contents[100], gnuplot_command[500], prefix[150], gnuplot_title[1000];

    if (tsp_seed > 0) 
        snprintf(prefix, sizeof(char)*150, "%ld_%d_%s", tsp_seed, tsp_inst.nnodes, tsp_alg_type);
    else
        snprintf(prefix, sizeof(char)*150, "%s_%s", tsp_file_name, tsp_alg_type);

    snprintf(plot_file_name, sizeof(char)*500, "%s/%s_%s", TSP_SOL_FOLDER, prefix, TSP_PLOT_FILE);  //where to save the plot
    snprintf(solution_file_name, sizeof(char)*500, "%s/%s_%s", TSP_SOL_FOLDER, prefix, TSP_SOLUTION_FILE);  //where to read the file from

    solution_file = fopen(solution_file_name, "r");
    coords_file = fopen(TSP_COORDS_FILE, "w");
    command_file = fopen(TSP_COMMAND_FILE, "w");

    if (solution_file == NULL || coords_file == NULL || command_file == NULL) {
        printf("Error with a file used to plot the solution.");
        exit(1);
    }

    // skip through the rows with the solution info
    while (rows_read < 6) if (fgetc(solution_file) =='\n') rows_read++;

    // build plot title with solution info
    snprintf(gnuplot_title, 1000, "Algorithm: %s; number of nodes: %d; cost: %.4f; time: %.4fs", tsp_alg_type, tsp_inst.nnodes, tsp_inst.best_cost, tsp_inst.best_time);

    // copy nodes coordinates into coords_file
    while (fgets(solution_contents, 100, solution_file)) fprintf(coords_file, "%s", solution_contents);
    // builds commands for gnuplot
    fprintf(command_file, "set term png\n");
    fprintf(command_file, "set output '%s'\n", plot_file_name);
    fprintf(command_file, "x=0.; y=0.\n");
    fprintf(command_file, "set title '%s'\n", gnuplot_title);
    fprintf(command_file, "set xlabel 'Starting node highlighted as black point'\n");
    fprintf(command_file, "set label at %f, %f point pointtype 7 pointsize 2\n", tsp_inst.coords[tsp_inst.best_solution[0]].x, tsp_inst.coords[tsp_inst.best_solution[0]].y);
    fprintf(command_file, "plot '%s' u (x=$2):(y=$3) w lp lc rgb 'blue' title ''", TSP_COORDS_FILE);

    fclose(solution_file);
    fclose(coords_file);
    fclose(command_file);

    // execute gnuplot commands and remove all intermediate files
    snprintf(gnuplot_command, sizeof(char)*500, "gnuplot %s", TSP_COMMAND_FILE);
    system(gnuplot_command);
    remove(TSP_COORDS_FILE);
    remove(TSP_COMMAND_FILE);
    
}

void tsp_save_intermediate_cost(const int list_index, const double cost) {
    
    FILE* file = fopen(tsp_intermediate_costs_files[list_index], "a");

    fprintf(file, "%f\n", cost);

    fclose(file);

}
#pragma endregion


#pragma region DEBUGGING TOOLS
void tsp_instance_info() {

    printf("--------------------\n");
    printf("Type of Instance: %s\n", ((tsp_seed == 0) ? "from file" : "random"));
    if (tsp_seed == 0) printf("File name: %s\n", tsp_file_name);
    else printf("Seed: %ld\n", tsp_seed);
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

void tsp_check_integrity(const int* path, const double cost) {

    double c_cost = 0;
    int first = path[0], error = 0;
    int* visited = (int*)calloc(tsp_inst.nnodes, sizeof(int));

    for (int i = 1; i < tsp_inst.nnodes; i++) {
        if (path[i] < 0 || path[i] >= tsp_inst.nnodes) { error = 1; break; }
        if (visited[path[i]]) { error = 2; break; }
        visited[path[i]] += 1;
        c_cost += tsp_inst.costs[path[i-1] * tsp_inst.nnodes + path[i]];
    }

    if (visited != NULL) { free(visited); visited = NULL; }

    c_cost += tsp_inst.costs[path[tsp_inst.nnodes - 1] * tsp_inst.nnodes + path[0]];

    if (error == 0 && fabs(c_cost - cost) > TSP_EPSILON) error = 3;

    if (error >= 1) {
        printf("\n\nINTEGRITY COMPROMISED - error_code: %d\n", error);
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