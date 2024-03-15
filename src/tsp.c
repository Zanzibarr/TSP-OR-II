#include "../include/tsp.h"

#pragma region GLOBALS DEFINITIONS
double tsp_initial_time = 0;
double tsp_total_time = 0;
int tsp_over_time = 0;
time_t tsp_time_limit = 0;

uint64_t tsp_seed = 0;

char tsp_alg_type[20] = "";
char tsp_file_name[100] = "";

tabu tsp_tabu_table;

pthread_t tsp_threads[N_THREADS];
int tsp_available_threads[N_THREADS];

int tsp_stoplight_update_sol = 0;
int tsp_mt_choice = 0;
#pragma endregion

#pragma region PRECOMPUTING
void tsp_precompute_sort_edges(tsp_instance* inst) { // precomputes the sort_edges "matrix"

    inst -> sort_edges = (int*)calloc(inst -> nnodes * (inst -> nnodes - 1), sizeof(int));
    tsp_entry* list = (tsp_entry*)calloc(inst -> nnodes, sizeof(tsp_entry));

    for (int i = 0; i < inst -> nnodes; i++) {

        for (int j = 0; j < inst -> nnodes; j++) {  //saving the entries to be sorted
            list[j].key = j;
            list[j].value = inst -> costs[i * inst -> nnodes + j];    //considering only the costs of the edges leaving node i
        }

        qsort((void*)list, (size_t)inst -> nnodes, sizeof(tsp_entry), compare_tsp_entries); //sort by cost of the edge
        
        for (int j = 1; j < inst -> nnodes; j++)
            inst -> sort_edges[i * (inst -> nnodes - 1) + j-1] = list[j].key;    //populate the ith row with the nodes ordered by increasing distance

    }

    if (list != NULL) { free(list); list = NULL; }

    #if TSP_VERBOSE >= 100
    tsp_check_sort_edges_integrity(inst);
    #endif

}

void tsp_precompute_costs(tsp_instance* inst) { //precomputes the costs "matrix"

    tsp_allocate_costs_space(inst);

    for (int i = 0; i < inst -> nnodes; i++) for (int j = 0; j < inst -> nnodes; j++)
        inst -> costs[i * inst -> nnodes + j] = tsp_compute_distance(inst, i, j);

}

double tsp_compute_distance(const tsp_instance* inst, int i, int j) { //euclidian distance between two points in the instance
    
    return sqrt(pow(inst -> coords[i].x - inst -> coords[j].x, 2) + pow(inst -> coords[i].y - inst -> coords[j].y, 2));
    
}

int compare_tsp_entries( const void* arg1, const void* arg2) { //compare function for the qsort

    double diff = ((tsp_entry*)arg1)->value - ((tsp_entry*)arg2)->value;
    return (fabs(diff) < TSP_EPSILON) ? 0 : ((diff < 0) ? -1 : 1);

}
#pragma endregion

#pragma region MULTITHREADING
void tsp_init_threads() {

    #if TSP_VERBOSE == 5
    printf("Initializing threads.\n");
    #endif

    for (int i = 0; i < N_THREADS; i++) tsp_available_threads[i] = 1;

}

int tsp_wait_for_thread() {

    #if TSP_VERBOSE == 5
    int rnd_index = (int)tsp_rnd_coord();
    printf("- %4d - Waiting for thread.\n", rnd_index);
    #endif

    while (1)
        for (int i = 0; i < N_THREADS; i++)
            if (tsp_available_threads[i]) {
                #if TSP_VERBOSE == 5
                printf("- %4d - Thread %d available.\n", rnd_index, i);
                #endif
                tsp_available_threads[i] = 0;
                return i;
            }

}

void tsp_free_thread(int index) {

    #if TSP_VERBOSE == 5
    printf("Freeing thread %d.\n", index);
    #endif

    pthread_join(tsp_threads[index], NULL);
    tsp_available_threads[index] = 1;

}

void tsp_wait_all_threads() {

    #if TSP_VERBOSE == 5
    printf("Waiting for all threads to finish.\n");
    #endif

    int free = 0;
    while (!free) {
        free = 1;
        for (int i = 0; i < N_THREADS; i++)
            if (!tsp_available_threads[i]) free = 0;
    }

    #if TSP_VERBOSE == 5
    printf("All threads finished.\n");
    #endif

}
#pragma endregion

#pragma region ALGORITHMS TOOLS
void tsp_check_best_sol(tsp_instance* inst, int* path, double cost, double time) { //update the best solution found so far

    if (cost > inst -> best_cost + TSP_EPSILON) return;

    while(!tsp_stoplight_update_sol) {
        #if TSP_VERBOSE == 5
        printf("----- Waiting to update best solution. -----\n");
        #endif
    }
    
    tsp_stoplight_update_sol = 0;

    if (cost < inst -> best_cost - TSP_EPSILON) {

        for (int i = 0; i < inst -> nnodes; i++) inst -> best_solution[i] = path[i];
        inst -> best_cost = cost;
        inst -> best_time = time;

        #if TSP_VERBOSE >= 10
        printf("New best solution\n");
        #endif

    }

    tsp_stoplight_update_sol = 1;

}

void tsp_reverse(int* path, int start, int end) { //reverse the array specified between two specified indexes

    int c;

    while (start < end) {
        c = path[start];
        path[start] = path[end];
        path[end] = c;
        start++; end--;
    }

}

int tsp_check_tabu(int node_1, int node_2) {

    for (int i = 0; i < TSP_TABU_SIZE; i++) if (tsp_tabu_table.table_1[i] == node_1 && tsp_tabu_table.table_2[i] == node_2) return 1;

    return 0;

}

void tsp_add_tabu(int node_1, int node_2) {

    tsp_tabu_table.table_1[tsp_tabu_table.pointer++] = node_1;
    tsp_tabu_table.table_2[tsp_tabu_table.pointer++] = node_2;

    tsp_tabu_table.pointer = tsp_tabu_table.pointer % TSP_TABU_SIZE;

}
#pragma endregion

#pragma region INIZIALIZATIONS
void tsp_init_defs(tsp_instance* inst) { //default values

    tsp_seed = 0;
    tsp_time_limit = (double)TSP_DEF_TL/1000;

    strcpy(tsp_file_name, "NONE");
    strcpy(tsp_alg_type, "greedy");

    inst -> nnodes = TSP_DEF_NNODES;

    struct timeval tv;
    gettimeofday(&tv, NULL);
    tsp_initial_time = ((double)tv.tv_sec)+((double)tv.tv_usec/1e+6);

    tsp_stoplight_update_sol = 1;

    tsp_tabu_table.pointer = 0;
    for (int i = 0; i < TSP_TABU_SIZE; i++) { tsp_tabu_table.table_1[i] = -1; tsp_tabu_table.table_2[i] = -1; }

    tsp_init_threads();

}

void tsp_init_solution(tsp_instance* inst) { //initialize the best solution

    tsp_allocate_best_sol_space(inst);
    inst -> best_cost = INFINITY;
    inst -> best_time = 0;

    tsp_over_time = 0;

}
#pragma endregion

#pragma region SAVING FILES

void tsp_save_solution(const tsp_instance* inst) {//save the best solution found in a file
    
    FILE *solution_file;

    char prefix[150], solution_file_name[500];

    if (tsp_seed > 0) 
        snprintf(prefix, sizeof(char)*150, "%ld_%d_%s", tsp_seed, inst -> nnodes, tsp_alg_type);
    else
        snprintf(prefix, sizeof(char)*150, "%s_%s", tsp_file_name, tsp_alg_type);
    snprintf(solution_file_name, sizeof(char)*500, "%s/%s_%s", TSP_SOL_FOLDER, prefix, TSP_SOLUTION_FILE);  //where to save the file

    solution_file = fopen(solution_file_name, "w");

    if (solution_file == NULL) {
        printf("Error writing the file for the solution.");
        exit(1);
    }

    fprintf(solution_file, "Algorithm: %s\n", tsp_alg_type);
    fprintf(solution_file, "Cost: %15.4f\n", inst -> best_cost);
    fprintf(solution_file, "Time: %15.4fs\n", inst -> best_time);
    fprintf(solution_file, "Total execution time: %15.4fs\n", tsp_total_time);
    fprintf(solution_file, "The algorithm %s exceeded the time limit%s\n", (tsp_over_time ? "has" : "hasn't"), (tsp_over_time ? " and has been stopped." : "."));
    fprintf(solution_file, "--------------------\n");
    for (int i = 0; i < inst -> nnodes; i++)
        fprintf(solution_file, "%4d %15.4f %15.4f\n", inst -> best_solution[i], inst->coords[inst -> best_solution[i]].x, inst->coords[inst -> best_solution[i]].y);
    fprintf(solution_file, "%4d %15.4f %15.4f\n", inst -> best_solution[0], inst ->coords[inst -> best_solution[0]].x, inst->coords[inst -> best_solution[0]].y);

    fclose(solution_file);

}

void tsp_plot_solution(const tsp_instance* inst) { //plot the best solution found

    int rows_read = 0;
    FILE *solution_file, *coords_file, *command_file;
    char plot_file_name[500], solution_file_name[500], solution_contents[100], gnuplot_command[500], prefix[150], gnuplot_title[1000];

    if (tsp_seed > 0) 
        snprintf(prefix, sizeof(char)*150, "%ld_%d_%s", tsp_seed, inst -> nnodes, tsp_alg_type);
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
    snprintf(gnuplot_title, 1000, "Algorithm: %s; number of nodes: %d; cost: %.4f; time: %.4fs", tsp_alg_type, inst -> nnodes, inst -> best_cost, inst -> best_time);

    // copy nodes coordinates into coords_file
    while (fgets(solution_contents, 100, solution_file)) fprintf(coords_file, "%s", solution_contents);
    // builds commands for gnuplot
    fprintf(command_file, "set term png\n");
    fprintf(command_file, "set output '%s'\n", plot_file_name);
    fprintf(command_file, "x=0.; y=0.\n");
    fprintf(command_file, "set title '%s'\n", gnuplot_title);
    fprintf(command_file, "set xlabel 'Starting node highlighted as black point'\n");
    fprintf(command_file, "set label at %f, %f point pointtype 7 pointsize 2\n", inst -> coords[inst -> best_solution[0]].x, inst -> coords[inst -> best_solution[0]].y);
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
#pragma endregion

#pragma region DEBUGGING TOOLS
void tsp_instance_info(const tsp_instance* inst) { //prints the instance info and problem's parameters

    printf("--------------------\n");
    printf("Type of Instance: %s\n", ((tsp_seed == 0) ? "from file" : "random"));
    if (tsp_seed == 0) printf("File name: %s\n", tsp_file_name);
    else printf("Seed: %ld\n", tsp_seed);
    printf("Time limit: %lds\n", tsp_time_limit);
    printf("Number of nodes: %4d\n", inst -> nnodes);
    printf("Edge weight type: ATT\n");
    printf("--------------------\n");
    printf("Algorithm: %s\n", tsp_alg_type);
    printf("--------------------\n");

    #if TSP_VERBOSE < 1000
    return;
    #endif

    printf("NODES:\n");
    for (int i = 0; i < inst -> nnodes; i++) printf("node[%4d]: (%15.4f, %15.4f)\n", i, inst -> coords[i].x, inst -> coords[i].y);
    printf("--------------------\n");
    printf("COSTS\n");
    for (int i = 0; i < inst -> nnodes; i++) for (int j = 0; j < inst -> nnodes; j++) printf("v%4d - v%4d : %15.4f\n", i, j, inst -> costs[i * inst -> nnodes + j]);
    printf("--------------------\n");
    printf("MIN EDGES\n");
    for (int i = 0; i < inst -> nnodes; i++) {
        for (int j = 0; j < inst -> nnodes - 1; j++)
            printf("%4d(%15.4f) ", inst -> sort_edges[i * (inst -> nnodes - 1) + j], inst -> costs[i * (inst -> nnodes) + inst -> sort_edges[i * (inst -> nnodes - 1) + j]]);
        printf("\n");
    }
    printf("--------------------\n");

}

void tsp_print_solution(const tsp_instance* inst) { //prints the solution to stdout

    printf("--------------------\nBEST SOLUTION:\n");
    printf("Cost: %15.4f\n", inst -> best_cost);
    printf("Time:\t%15.4fs\n", inst -> best_time);
    printf("Execution time %8.4fs\n", tsp_total_time);
    #if TSP_VERBOSE >= 500
    for (int i = 0; i < inst -> nnodes; i++) printf("%d -> ", inst -> best_solution[i]);
    printf("%d\n", inst -> best_solution[0]);
    #endif
    if (tsp_over_time) printf("The algorithm exceeded the time limit and has been stopped.\n");
    printf("--------------------\n");

}

void tsp_check_sort_edges_integrity(const tsp_instance* inst) { //debugging tool
    
    for (int i = 0; i < inst -> nnodes; i++) {
        double min_cost = 0;
        for (int j = 0; j < inst -> nnodes - 1; j++) {
            double checked_cost = inst -> costs[i * inst -> nnodes + inst -> sort_edges[i * (inst -> nnodes - 1) + j]];

            if (checked_cost < min_cost - TSP_EPSILON) {
                printf("SORT_EDGES INTEGRITY COMPROMISED\n");
                exit(1);
            }
            
            min_cost = checked_cost;
        }
    }

    #if TSP_VERBOSE >= 100
    printf("sort_edges integrity check passed.\n");
    #endif

}

void tsp_check_integrity(const tsp_instance* inst, double cost, int* path) { //debugging tool

    double c_cost = 0;
    int first = path[0], error = 0;
    int* visited = (int*)calloc(inst -> nnodes, sizeof(int));

    for (int i = 1; i < inst -> nnodes; i++) {
        if (path[i] < 0 || path[i] >= inst -> nnodes) { error = 1; break; }
        if (visited[path[i]]) { error = 2; break; }
        visited[path[i]] += 1;
        c_cost += inst -> costs[path[i-1] * inst -> nnodes + path[i]];
    }

    if (visited != NULL) { free(visited); visited = NULL; }

    c_cost += inst -> costs[path[inst -> nnodes - 1] * inst -> nnodes + path[0]];

    if (error == 0 && fabs(c_cost - cost) > TSP_EPSILON) error = 3;

    if (error >= 1) {
        printf("\n\nINTEGRITY COMPROMISED - error_code: %d\n", error);
        if (error == 1) printf("\nNon-existent node in path.\n");
        if (error == 2) printf("\nDouble node in path.\n");
        if (error == 3) printf("\nCost: %.10f, Checked cost: %.10f, Difference: %.10f, Threshold: %.10f\n", cost, c_cost, fabs(c_cost - cost), TSP_EPSILON);
        exit(1);
    }

}
#pragma endregion

#pragma region MEMORY MANAGEMENT
void tsp_allocate_coords_space(tsp_instance* inst) { //dinamically allocate the space for the coords list

    if (!inst -> nnodes) { printf("The nnodes variable hasn't been assigned yet."); exit(1); }

    inst -> coords = (tsp_pair*)calloc(inst -> nnodes, sizeof(tsp_pair));

}

void tsp_allocate_costs_space(tsp_instance* inst) { //dinamicallu allocate the space for the costs matrix

    if (!inst -> nnodes) { printf("The nnodes variable hasn't been assigned yet."); exit(1); }

    inst -> costs = (double*)calloc(inst -> nnodes * inst -> nnodes, sizeof(double));

}

void tsp_allocate_best_sol_space(tsp_instance* inst) { //dinamically allocate the space for the best solution list

    if (!inst -> nnodes) { printf("The nnodes variable hasn't been assigned yet."); exit(1); }

    inst -> best_solution = (int*)calloc(inst -> nnodes, sizeof(int));
    
}

void tsp_free_instance(tsp_instance* inst) { //frees the dinamically allocated memory

    if (inst -> coords != NULL) { free(inst -> coords); inst -> coords = NULL; }
    if (inst -> costs != NULL) { free(inst -> costs); inst -> costs = NULL; }
    if (inst -> sort_edges != NULL) { free(inst -> sort_edges); inst -> sort_edges = NULL; }
    if (inst -> best_solution != NULL) { free(inst -> best_solution); inst -> best_solution = NULL; }

}
#pragma endregion

#pragma region USEFUL METHODS
void tsp_init_rand() { for (int i = 0; i < 100; i++) rand(); }  //fixing first random values being small

double tsp_rnd_coord() { return (double)rand()/RAND_MAX*TSP_GRID_SIZE; }  //generate a random number between 0 and GRID_SIZE

double tsp_time_elapsed() { //time elapsed since the beginning of the execution of the chosen algorithm

    struct timeval tv;

    gettimeofday(&tv, NULL);
    return ((double)tv.tv_sec)+((double)tv.tv_usec/1e+6) - tsp_initial_time;

}
#pragma endregion