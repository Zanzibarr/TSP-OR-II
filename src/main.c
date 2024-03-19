#include "../include/inst_gen.h"
#include "../include/algorithms.h"

/**
 * Instructions to use the program
*/
void tsp_help() {

    printf("Use:\n");
    printf("%s <file_name> : to specify a file to obtain the TPS values from.\n", TSP_FILE_P);
    printf("%s <int> : specify the time limit in seconds. (default: %4ds)\n", TSP_TIME_LIMIT, (int)TSP_DEF_TL);
    printf("%s <int> : (default type of instance) specify the seed to use to create random TPS data (the seed 0 cannot be used due to implementation choices).\n", TSP_SEED);
    printf("%s <int> : specity the number of nodes in the problem (default: %4d).\n", TSP_NNODES, TSP_DEF_NNODES);
    printf("%s <str> : Type of algorithm to use ([greedy, g2opt, g2opt-best, tabu]), (default: greedy).\n", TSP_ALGORITHM);

    exit(0);

}

/**
 * Parse the command line arguments to prepare the instance and the problem's parameters
*/
void tsp_parse_cmd(const int argc, const char** argv, tsp_instance* inst) {

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
        else if (!strcmp(argv[i], TSP_TIME_LIMIT)) { tsp_time_limit = atof(argv[++i]); }
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
    tsp_precompute_sort_edges(inst);

}

/**
 * Solve the instance based on the type of the algorithm specified
*/
void tsp_solve(tsp_instance* inst) {

    int result = 0;
    tsp_init_solution(inst);

    if (!strcmp(tsp_alg_type, "greedy")) result = tsp_solve_greedy(inst, NULL);                             //greedy
    else if(!strcmp(tsp_alg_type, "g2opt")) result = tsp_solve_greedy(inst, tsp_find_2opt_swap);            //greedy with 2opt (first swap policy)
    else if(!strcmp(tsp_alg_type, "g2opt-best")) result = tsp_solve_greedy(inst, tsp_find_2opt_best_swap);  //greedy with 2opt (best  swap policy)
    else if(!strcmp(tsp_alg_type, "tabu")) result = tsp_solve_tabu(inst);                                   //tabu
    
    else {
        printf("Error choosing the algorithm to use.");
        exit(1);
    }

    tsp_total_time = tsp_time_elapsed();

    if (result == -1)
        tsp_over_time = 1;

}

int main(int argc, const char** argv) {

    tsp_instance inst;
    tsp_init_defs(&inst);

    tsp_parse_cmd(argc, argv, &inst);

    #if TSP_VERBOSE >= 0
    tsp_instance_info(&inst);
    #endif

    tsp_solve(&inst);   //algorithm to find optimal(ish) solutions
    
    #if TSP_VERBOSE >= 0
    tsp_print_solution(&inst);
    #endif
    tsp_save_solution(&inst);
    tsp_plot_solution(&inst);
    
    tsp_free_instance(&inst);   //frees the dinamically allocated memory

}
