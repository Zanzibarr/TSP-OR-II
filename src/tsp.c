#include "../include/tsp.h"

#pragma region GLOBALS DEFINITIONS
double tsp_initial_time = 0;
double tsp_total_time = 0;
int tsp_over_time = 0;
double tsp_time_limit = 0;

uint64_t tsp_seed = 0;

char tsp_alg_type[20] = "";
char tsp_file_name[100] = "";

tsp_tabu tsp_tabu_tables[N_THREADS];
//tsp_tabu_2 tsp_tabu_tables_2[N_THREADS];

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

        #if TSP_VERBOSE == 9
        printf("%10.4f New best solution : %15.4f\n", tsp_time_elapsed(), cost);
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

double tsp_dinamic_tenue(int counter) {

    return (1 * sin((double)counter * 1/50) + 1) * TSP_TABU_TENURE;
    //return sin((double)counter/20) * TSP_TABU_TENURE / 1.1 + TSP_TABU_TENURE;

}

int tsp_check_tabu(int t_index, int from, int to) {

    if (from > to) { int c = from; from = to; to = c;}

    //return (tsp_tabu_table[t_index].table[node_1 * tsp_tabu_table[t_index].size + node_2] >= 0) && (tsp_tabu_table[t_index].counter - tsp_tabu_table[t_index].table[node_1 * tsp_tabu_table[t_index].size + node_2] < sin((double)tsp_tabu_table[t_index].counter/50) * TSP_TABU_TENURE / 2 + TSP_TABU_TENURE);
    //return (tsp_tabu_table[t_index].table[node_1 * tsp_tabu_table[t_index].size + node_2] >= 0) && (tsp_tabu_table[t_index].counter - tsp_tabu_table[t_index].table[node_1 * tsp_tabu_table[t_index].size + node_2] < TSP_TABU_TENURE);

    //printf("from %d, to %d, list[from].node_1 %d, list[from].counter_1 %d, list[from].node_2 %d, list[from].counter_2 %d\n", from, to, tsp_tabu_tables[t_index].list[from].node_1, tsp_tabu_tables[t_index].list[from].counter_1, tsp_tabu_tables[t_index].list[from].node_2, tsp_tabu_tables[t_index].list[from].counter_2);
    //printf("Counter: %d", tsp_tabu_tables[t_index].counter);

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

    //printf("TABU: tenue: %10.4f counter: %4d\n", sin((double)tsp_tabu_table.counter/10) * TSP_TABU_TENURE / 2 + TSP_TABU_TENURE, tsp_tabu_table.counter);

    //tsp_tabu_tables[t_index].table[node_1 * tsp_tabu_tables[t_index].size + node_2] = tsp_tabu_tables[t_index].counter++;

    if (tsp_tabu_tables[t_index].list[from].counter_1 != -1 && tsp_tabu_tables[t_index].counter - tsp_tabu_tables[t_index].list[from].counter_1 < tsp_dinamic_tenue(tsp_tabu_tables[t_index].counter)) {

        tsp_tabu_tables[t_index].list[from].node_2 = to;
        tsp_tabu_tables[t_index].list[from].counter_2 = tsp_tabu_tables[t_index].counter++;

    } else {

        tsp_tabu_tables[t_index].list[from].node_1 = to;
        tsp_tabu_tables[t_index].list[from].counter_1 = tsp_tabu_tables[t_index].counter++;

    }

}
/**
int tsp_check_tabu_2(int t_index, int node_1, int node_2, int node_3, int node_4) {

    if (node_1 > node_3) { int temp = node_1; node_1 = node_3; node_3 = temp; temp = node_2; node_2 = node_4; node_4 = temp; }
    if (node_1 > node_2) { int temp = node_1; node_1 = node_2; node_2 = temp; }
    if (node_3 > node_4) { int temp = node_3; node_3 = node_4; node_4 = temp; }

    for (int i = 0; i < tsp_tabu_tables_2[t_index].list[node_1].size; i++) {
        if (
            tsp_tabu_tables_2[t_index].list[node_1].list[i].node_2 == node_2 &&
            tsp_tabu_tables_2[t_index].list[node_1].list[i].node_3 == node_3 &&
            tsp_tabu_tables_2[t_index].list[node_1].list[i].node_4 == node_4 &&
            tsp_tabu_tables_2[t_index].counter - tsp_tabu_tables_2[t_index].list[node_1].list[i].counter < tsp_dinamic_tenue(tsp_tabu_tables_2[t_index].counter)
        ) return 1;
    }

    return 0;

}

void tsp_add_tabu_2(int t_index, int node_1, int node_2, int node_3, int node_4) {

    if (node_1 > node_3) { int temp = node_1; node_1 = node_3; node_3 = temp; temp = node_2; node_2 = node_4; node_4 = temp; }
    if (node_1 > node_2) { int temp = node_1; node_1 = node_2; node_2 = temp; }
    if (node_3 > node_4) { int temp = node_3; node_3 = node_4; node_4 = temp; }

    //printf("New tabu %d: %d %d %d %d\n", tsp_tabu_tables_2[t_index].counter, node_1, node_2, node_3, node_4);

    tsp_tabu_tables_2[t_index].list[node_1].list[tsp_tabu_tables_2[t_index].list[node_1].size-1] = (tsp_tabu_quadruple){tsp_tabu_tables_2[t_index].counter++, node_2, node_3, node_4};

    tsp_tabu_quadruple* tmp = (tsp_tabu_quadruple*)calloc(++tsp_tabu_tables_2[t_index].list[node_1].size, sizeof(tsp_tabu_quadruple));
    for (int i = 0; i < tsp_tabu_tables_2[t_index].list[node_1].size-1; i++) {
        //printf("tabu_list[%d].list[%d] = {%d, %d, %d, %d}\n", node_1, i, tsp_tabu_tables_2[t_index].list[node_1].list[i].counter, tsp_tabu_tables_2[t_index].list[node_1].list[i].node_2, tsp_tabu_tables_2[t_index].list[node_1].list[i].node_3, tsp_tabu_tables_2[t_index].list[node_1].list[i].node_4);
        //printf("tabu_list[%d].size = %d\n", node_1, tsp_tabu_tables_2[t_index].list[node_1].size);
        tmp[i].counter = tsp_tabu_tables_2[t_index].list[node_1].list[i].counter;
        tmp[i].node_2 = tsp_tabu_tables_2[t_index].list[node_1].list[i].node_2;
        tmp[i].node_3 = tsp_tabu_tables_2[t_index].list[node_1].list[i].node_3;
        tmp[i].node_4 = tsp_tabu_tables_2[t_index].list[node_1].list[i].node_4;
    }
    free(tsp_tabu_tables_2[t_index].list[node_1].list);
    tsp_tabu_tables_2[t_index].list[node_1].list = tmp;

}*/
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

    tsp_init_threads();

}

void tsp_init_solution(tsp_instance* inst) { //initialize the best solution

    tsp_allocate_best_sol_space(inst);
    inst -> best_cost = INFINITY;
    inst -> best_time = 0;

    tsp_over_time = 0;

    for (int thread = 0; thread < N_THREADS; thread++) {
        tsp_tabu_tables[thread].list = (tsp_tabu_entry*)calloc(inst -> nnodes, sizeof(tsp_tabu_entry));
        for (int i = 0; i < inst -> nnodes; i++) { tsp_tabu_tables[thread].list[i] = (tsp_tabu_entry){-1, -1, -1, -1}; }
    }
    /**
    for (int thread = 0; thread < N_THREADS; thread++) {
        tsp_tabu_tables_2[thread].list = (tsp_tabu_entry_2*)calloc(inst -> nnodes, sizeof(tsp_tabu_entry_2));
        for (int i = 0; i < inst -> nnodes; i++) { tsp_tabu_tables_2[thread].list[i].size = 1; tsp_tabu_tables_2[thread].list[i].list = (tsp_tabu_quadruple*)calloc(1, sizeof(tsp_tabu_quadruple)); }
    }*/

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
    snprintf(solution_file_name, sizeof(char)*500, "%s/%s%s_%s", TSP_SOL_FOLDER, prefix, (tsp_mt_choice) ? "_mt" : "", TSP_SOLUTION_FILE);  //where to save the file

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

    snprintf(plot_file_name, sizeof(char)*500, "%s/%s%s_%s", TSP_SOL_FOLDER, prefix, (tsp_mt_choice) ? "_mt" : "", TSP_PLOT_FILE);  //where to save the plot
    snprintf(solution_file_name, sizeof(char)*500, "%s/%s%s_%s", TSP_SOL_FOLDER, prefix, (tsp_mt_choice) ? "_mt" : "", TSP_SOLUTION_FILE);  //where to read the file from

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
    printf("Time limit: %10.4f\n", tsp_time_limit);
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

    for (int thread = 0; thread < N_THREADS; thread++)
        if (tsp_tabu_tables[thread].list != NULL) { free(tsp_tabu_tables[thread].list); tsp_tabu_tables[thread].list = NULL; }
    /**
    for (int thread = 0; thread < N_THREADS; thread++)
        if (tsp_tabu_tables_2[thread].list != NULL) {
            for (int i = 0; i < inst -> nnodes; i++) {
                if (tsp_tabu_tables_2[thread].list[i].list != NULL) { free(tsp_tabu_tables_2[thread].list[i].list); tsp_tabu_tables_2[thread].list[i].list = NULL; }
            }
            free(tsp_tabu_tables_2[thread].list); tsp_tabu_tables_2[thread].list = NULL;
        }
    */
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