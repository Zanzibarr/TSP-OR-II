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
    if (TSP_VERBOSE >= 0) tsp_instance_info(&inst);

    if (TSP_VERBOSE > 0) printf("Starting the execution of the %s algorithm\n", tsp_alg_type);

    tsp_solve(&inst);   //algorithm to find optimal(ish) solutions
    
    if (TSP_VERBOSE >= 0) tsp_print_solution(&inst);
    tsp_save_solution(&inst);
    tsp_plot_solution(&inst);
    
    tsp_free_instance(&inst);   //frees the dinamically allocated memory

}

void tsp_init_defs(tsp_instance* inst) {  //default values

    tsp_seed = 0;
    tsp_time_limit = TSP_DEF_TL;

    strcpy(tsp_file_name, "NONE");
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

    if (TSP_VERBOSE == 0) return;

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

    int result = 0;
    tsp_init_solution(inst);

    if (!strcmp(tsp_alg_type, "greedy")) result = tsp_solve_greedy(inst, 0);
    else if(!strcmp(tsp_alg_type, "g2opt")) result = tsp_solve_greedy(inst, 1);
    
    else {
        printf("Error choosing the algorithm to use.");
        exit(1);
    }

    tsp_total_time = tsp_time_elapsed();

    if (result == -1)
        tsp_over_time = 1;

}

void tsp_print_solution(const tsp_instance* inst) {

    printf("--------------------\nBEST SOLUTION:\n");
    printf("Cost: %f\n", inst -> tsp_best_cost);
    printf("Time: %fs\n", inst -> tsp_best_time);
    for (int i = 0; i < inst -> nnodes; i++) printf("%d->", inst -> tsp_best_solution[i]);
    printf("%d\n", inst -> tsp_best_solution[0]);
    printf("Total execution time: %fs\n", tsp_total_time);
    if (tsp_over_time) printf("The algorithm exceeded the time limit and has been stopped.\n");
    printf("--------------------\n");

}

void tsp_save_solution(const tsp_instance* inst) {  //save the best solution found in a file
    
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
    fprintf(solution_file, "Cost: %f\n", inst -> tsp_best_cost);
    fprintf(solution_file, "Time: %fs\n", inst -> tsp_best_time);
    fprintf(solution_file, "Total execution time: %fs\n", tsp_total_time);
    fprintf(solution_file, "The algorithm %s exceeded the time limit%s\n", (tsp_over_time ? "has" : "hasn't"), (tsp_over_time ? " and has been stopped." : "."));
    fprintf(solution_file, "--------------------\n");
    for (int i = 0; i < inst -> nnodes; i++)
        fprintf(solution_file, "%d %f %f\n", inst -> tsp_best_solution[i], inst->coords[inst -> tsp_best_solution[i]].x, inst->coords[inst -> tsp_best_solution[i]].y);
    fprintf(solution_file, "%d %f %f\n", inst -> tsp_best_solution[0], inst ->coords[inst -> tsp_best_solution[0]].x, inst->coords[inst -> tsp_best_solution[0]].y);

    fclose(solution_file);

}

void tsp_plot_solution(const tsp_instance* inst) {  //plot the best solution found

    int rows_read = 0;
    FILE *solution_file, *coords_file, *command_file;
    char plot_file_name[500], solution_file_name[500], solution_contents[100], gnuplot_command[500], prefix[150];

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

    // copy nodes coordinates into coords_file
    while (fgets(solution_contents, 100, solution_file)) fprintf(coords_file, "%s", solution_contents);
    // builds commands for gnuplot
    fprintf(command_file, "set term png\nset output '%s'\nx=0.; y=0.\nplot '%s' u (x=$2):(y=$3) w lp", plot_file_name, TSP_COORDS_FILE);

    fclose(solution_file);
    fclose(coords_file);
    fclose(command_file);

    // execute gnuplot commands and remove all intermediate files
    snprintf(gnuplot_command, sizeof(char)*500, "gnuplot %s", TSP_COMMAND_FILE);
    system(gnuplot_command);
    remove(TSP_COORDS_FILE);
    remove(TSP_COMMAND_FILE);
    
}

void tsp_help() {   //instructions to use the program

    printf("Use:\n");
    printf("%s <file_name> : to specify a file to obtain the TPS values from.\n", TSP_FILE_P);
    printf("%s <int> : specify the time limit in seconds. (default: %ds)\n", TSP_TIME_LIMIT, (int)TSP_DEF_TL);
    printf("%s <int> : (default type of instance) specify the seed to use to create random TPS data (the seed 0 cannot be used due to implementation choices).\n", TSP_SEED);
    printf("%s <int> : specity the number of nodes in the problem (default: %d).\n", TSP_NNODES, TSP_DEF_NNODES);
    printf("%s <str> : Type of algorithm to use ([greedy, g2opt]), (default: greedy).\n", TSP_ALGORITHM);

    exit(0);

}