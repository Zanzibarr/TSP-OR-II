#include "tsp.h"
#include "algorithms.h"
#include "inst_gen.h"

void tsp_init_defs(tsp_instance* inst);
void tsp_parse_cmd(const int argc, const char** argv, tsp_instance* inst);
void tsp_instance_info(const tsp_instance* inst);

void tsp_solve(tsp_instance* inst);

void tsp_print_solution(const tsp_instance* inst);
void tsp_save_solution(const tsp_instance* inst);
void tsp_plot_solution(const tsp_instance* inst);

void tsp_help();

int main(int argc, const char** argv) {

    tsp_instance inst;
    tsp_init_defs(&inst);

    tsp_parse_cmd(argc, argv, &inst);
    if (tsp_verbose >= 0) tsp_instance_info(&inst);

    if (tsp_verbose > 0) printf("Starting the execution of the %s algorithm\n", tsp_alg_type);

    tsp_solve(&inst);   //algorithm to find optimal(ish) solutions
    
    if (tsp_verbose >= 0) tsp_print_solution(&inst);
    tsp_save_solution(&inst);
    tsp_plot_solution(&inst);
    
    tsp_free_instance(&inst);   //frees the dinamically allocated memory

}

void tsp_init_defs(tsp_instance* inst) {  //default values

    tsp_seed = 0;
    tsp_time_limit = TSP_DEF_TL;
    tsp_verbose = 0;

    strcpy(tsp_file_name, "NONE");
    strcpy(tsp_edge_weight_type, "ATT");

    strcpy(tsp_alg_type, "greedy");

    inst -> nnodes = TSP_DEF_NNODES;

}

void tsp_parse_cmd(const int argc, const char** argv, tsp_instance* inst) { //parse the command line arguments to prepare the instance and the problem's parameters

    int check = -1;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], TSP_FILE_P)) {

            if (check == 1) {
                printf("Cannot parse both a seed and a file_name.\n");
                exit(1);
            }
            
            strcpy(tsp_file_name, argv[++i]);
            check = 0;
            
        }
        else if (!strcmp(argv[i], TSP_SEED)) {

            if (check == 0) {
                printf("Cannot parse both a seed and a file_name.\n");
                exit(1);
            }

            tsp_seed = atoi(argv[++i]);
            srand(tsp_seed);
            tsp_init_rand();
            check = 1;
            
        }
        else if (!strcmp(argv[i], TSP_TIME_LIMIT)) { tsp_time_limit = atoi(argv[++i]); }
        else if (!strcmp(argv[i], TSP_NNODES)) {

            if (check == 0) {
                printf("Cannot parse the number of nodes if a file_name is specified.\n");
                exit(1);
            }
            
            inst -> nnodes = atoi(argv[++i]);
            check = 1;
            
        }
        else if (!strcmp(argv[i], TSP_HELP)) { tsp_help(); }
        else if (!strcmp(argv[i], TSP_QUIET)) { tsp_verbose = -1; }
        else if (!strcmp(argv[i], TSP_VERBOSE)) { tsp_verbose = 1; }
        else if (!strcmp(argv[i], TSP_ALGORITHM)) { strcpy(tsp_alg_type, argv[++i]); }
        else { printf("Error parsing the command line arguments; use %s to view the command line options.", TSP_HELP); exit(1); }
    }

    if (tsp_seed > 0)   //if the seed is not at 0 (default value), then a seed has been specified -> generate instance randomly
        tsp_gen_random_instance(inst);
    else    //no seed specified: generating instance from filename given
        tsp_gen_instance_from_file(inst);

    tsp_precompute_costs(inst);
    tsp_precompute_min_edges(inst);

}

void tsp_instance_info(const tsp_instance* inst) {  //prints the instance info and problem's parameters

    printf("--------------------\n");
    printf("Type of Instance: %s\n", ((tsp_seed == 0) ? "from file" : "random"));
    if (tsp_seed == 0) printf("File name: %s\n", tsp_file_name);
    else printf("Seed: %ld\n", tsp_seed);
    printf("Time limit: %lds\n", tsp_time_limit);
    printf("Number of nodes: %d\n", inst -> nnodes);
    printf("Edge weight type: ATT\n");
    printf("--------------------\n");
    printf("Algorithm: %s\n", tsp_alg_type);
    printf("--------------------\n");

    if (tsp_verbose == 0) return;

    printf("NODES:\n");
    for (int i = 0; i < inst -> nnodes; i++) printf("node[%d]: (%f, %f)\n", i, inst -> coords[i].x, inst -> coords[i].y);
    printf("--------------------\n");
    printf("COSTS\n");
    for (int i = 0; i < inst -> nnodes; i++) for (int j = 0; j < inst -> nnodes; j++) printf("v%d - v%d : %f\n", i, j, inst -> costs[i][j]);
    printf("--------------------\n");
    printf("MIN EDGES\n");
    for (int i = 0; i < inst -> nnodes; i++) {
        for (int j = 0; j < inst -> nnodes - 1; j++)
            printf("%d ", inst -> min_edges[i][j]);
        printf("\n");
    }
    printf("--------------------\n");

}

void tsp_solve(tsp_instance* inst) {  //solve the instance based on the type of the algorithm specified

    tsp_init_solution(inst);

    int result = 0;

    if (!strcmp(tsp_alg_type, "greedy")) result = tsp_solve_greedy(inst, 0);
    else if(!strcmp(tsp_alg_type, "g2opt")) result = tsp_solve_greedy(inst, 1);
    
    else {
        printf("Error choosing the algorithm to use.");
        exit(1);
    }

    tsp_total_time = tsp_time_elapsed();

    if (result == -1) {
        printf("The algorithm has been stopped since it exceeded the time limit.\n");
    }

}

void tsp_print_solution(const tsp_instance* inst) {

    printf("--------------------\nBEST SOLUTION:\n");
    printf("Cost: %f\n", inst -> tsp_best_cost);
    printf("Time: %fs\n", inst -> tsp_best_time);
    for (int i = 0; i < inst -> nnodes; i++) printf("%d->", inst -> tsp_best_solution[i]);
    printf("%d\n", inst -> tsp_best_solution[0]);
    printf("Total execution time: %fs\n", tsp_total_time);
    printf("--------------------\n");

}

void tsp_save_solution(const tsp_instance* inst) {  //save the best solution found in a file
    
    FILE *solution_file;

    solution_file = fopen("solution_file.txt", "w");
    fprintf(solution_file, "BEST SOLUTION:\n");
    fprintf(solution_file, "Cost: %f\n", inst->tsp_best_cost);
    fprintf(solution_file, "Time for best solution: %fs\n", inst->tsp_best_time);
    fprintf(solution_file, "Total execution time: %fs\n", tsp_total_time);
    fprintf(solution_file, "Algorithm used: %s\n", tsp_alg_type);
    for (int i = 0; i < inst -> nnodes; i++)
        fprintf(solution_file, "%d %f %f\n", inst -> tsp_best_solution[i], inst->coords[inst -> tsp_best_solution[i]].x, inst->coords[inst -> tsp_best_solution[i]].y);
    fprintf(solution_file, "%d %f %f\n", inst -> tsp_best_solution[0], inst ->coords[inst -> tsp_best_solution[0]].x, inst->coords[inst -> tsp_best_solution[0]].y);

    fclose(solution_file);

}

void tsp_plot_solution(const tsp_instance* inst) {  //plot the best solution found

    int rows_read, character, remove_success;
    FILE *command_file, *solution_file, *coords_file;
    char *solution_contents, *gnuplot_command;

    remove(TSP_PLOT_FILE);
    solution_file = fopen(TSP_SOLUTION_FILE, "r");
    coords_file = fopen(TSP_COORDS_FILE, "w");
    command_file = fopen(TSP_COMMAND_FILE, "w");
    solution_contents = malloc(100);
    gnuplot_command = malloc(100);

    // skip through the rows with the solution info
    rows_read = 0;
    while (rows_read < 5) {
        character = fgetc(solution_file);
        if (character=='\n') rows_read++;
    }

    // copy nodes coordinates into coords_file
    while (fgets(solution_contents, 100, solution_file)) fprintf(coords_file, "%s", solution_contents);
    // builds commands for gnuplot
    fprintf(command_file, "set term png\nset output '%s'\nx=0.; y=0.\nplot 'coords_file.txt' u (x=$2):(y=$3) w lp", TSP_PLOT_FILE);

    free(solution_contents);
    fclose(solution_file);
    fclose(coords_file);
    fclose(command_file);

    // execute gnuplot commands and remove all intermediate files
    strcpy(gnuplot_command, "gnuplot ");
    strcpy(gnuplot_command+8, TSP_COMMAND_FILE);
    system(gnuplot_command);
    remove(TSP_COORDS_FILE);
    remove(TSP_COMMAND_FILE);
    free(gnuplot_command);
    
}

void tsp_help() {   //instructions to use the program

    printf("Use:\n");
    printf("%s <file_name> : to specify a file to obtain the TPS values from.\n", TSP_FILE_P);
    printf("%s <int> : specify the time limit in seconds. (default: %ds)\n", TSP_TIME_LIMIT, (int)TSP_DEF_TL);
    printf("%s <int> : (default type of instance) specify the seed to use to create random TPS data (the seed 0 cannot be used due to implementation choices).\n", TSP_SEED);
    printf("%s <int> : specity the number of nodes in the problem (default: %d).\n", TSP_NNODES, TSP_DEF_NNODES);
    printf("%s <str> : Type of algorithm to use ([greedy, g2opt]), (default: greedy).\n", TSP_ALGORITHM);
    printf("%s : Logging option (quiet).\n", TSP_QUIET);
    printf("%s : Logging option (verbose).\n", TSP_VERBOSE);

    exit(0);

}