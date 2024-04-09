#include "../include/inst_gen.h"
#include "../include/algorithms.h"

/**
 * @brief Capture the Ctrl+C signal and terminate peacefully the program
 * 
 * @param signum The signal captured
*/
void signal_callback_handler(const int signum) {

    printf("\n\n---------------------------------------------");
    printf("\n- Caught ctrl+C signal, exiting peacefully. -");
    printf("\n---------------------------------------------\n");
    
    tsp_check_best_sol(NULL, INFINITY, 0);    //wait for eventual best solution updates in the queue
    tsp_time_limit = 0; //signal all running events to stop

    printf("\n---------------------------------------------");
    printf("\n-  Waiting for all processes to terminate.  -");
    printf("\n---------------------------------------------\n\n");

    tsp_wait_all_threads(); //wait for all process to stop properly

    printf("\n\n---------------------------------------------");
    printf("\n-      All processes ended peacefully.      -");
    printf("\n---------------------------------------------\n\n");

    tsp_total_time = tsp_time_elapsed();
    tsp_forced_termination = 1;
    
    if (tsp_verbose >= 0) tsp_print_solution();
    tsp_save_solution();
    tsp_plot_solution();
    
    tsp_free_instance(); //frees the dinamically allocated memory and other finishing operations
    
    exit(0);

}

/**
 * @brief Instructions to use the program
*/
void tsp_help() {

    printf("Use:\n");
    printf("%s <file_name> : to specify a file to obtain the TPS values from.\n", TSP_PARSING_FILE);
    printf("%s <int> : specify the seed to use to create random TPS data (the seed 0 cannot be used due to implementation choices).\n", TSP_PARSING_SEED);
    printf("%s <int> : specity the number of nodes in the problem (default: %4d).\n", TSP_PARSING_NNODES, TSP_DEF_NNODES);
    printf("%s <str> : Type of algorithm to use ([", TSP_PARSING_ALGORITHM);
    for (int i=0; i<TSP_ALG_NUMBER-1; i++) printf("%s, ", tsp_algorithms[i]);
    printf("%s]), (default: %s)\n", tsp_algorithms[TSP_ALG_NUMBER-1], tsp_algorithms[0]);
    printf("%s <int> : Verbose parameter to control the amount of output (default: 100).\n", TSP_PARSING_VERBOSE);
    printf("%s <int> : Tenure for the tabu algorithm (default: %4d).\n", TSP_PARSING_TENURE, TSP_DEF_TABU_TENURE);
    printf("%s <int> : Amplitude parameter for the dinamic tenure (default: %d).\n", TSP_PARSING_TENURE_A, tsp_tabu_tenure_a);
    printf("%s <double> : Frequency parameter for the dinamic tenure (default: %10.4f).\n", TSP_PARSING_TENURE_F, tsp_tabu_tenure_f);
    printf("%s <int> : specify the time limit in seconds (default: %4ds).\n", TSP_PARSING_TIME_LIMIT, (int)TSP_DEF_TL);

    exit(0);

}

/**
 * @brief Parse the command line arguments to prepare the instance and the problem's parameters
 * 
 * @param argv List of strings
 * @param argc Length of the list of strings
*/
void tsp_parse_cmd(const char** argv, const int argc) {

    int check = -1;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], TSP_PARSING_FILE)) {

            if (check == 1) {
                printf("Cannot parse both a seed and a file_name.\n");
                exit(1);
            }
            
            strcpy(tsp_file_name, argv[++i]);
            check = 0;
            
        }
        else if (!strcmp(argv[i], TSP_PARSING_SEED)) {

            if (check == 0) {
                printf("Cannot parse both a seed and a file_name.\n");
                exit(1);
            }

            tsp_seed = atoi(argv[++i]);
            srand(tsp_seed);
            tsp_init_rand();
            check = 1;
            
        }
        else if (!strcmp(argv[i], TSP_PARSING_TIME_LIMIT)) { tsp_time_limit = atof(argv[++i]); }
        else if (!strcmp(argv[i], TSP_PARSING_NNODES)) {

            if (check == 0) {
                printf("Cannot parse the number of nodes if a file_name is specified.\n");
                exit(1);
            }
            
            tsp_inst.nnodes = atoi(argv[++i]);
            check = 1;
            
        }
        else if (!strcmp(argv[i], TSP_PARSING_HELP)) { tsp_help(); }
        else if (!strcmp(argv[i], TSP_PARSING_ALGORITHM)) {

            if (!strcmp(argv[i+1], "benders")) snprintf(tsp_alg_type, 20, "cplex-benders");
            else if (!strcmp(argv[i+1], "cplex")) snprintf(tsp_alg_type, 20, "cplex-base");
            else strcpy(tsp_alg_type, argv[i+1]);

            i++;

        }
        else if (!strcmp(argv[i], TSP_PARSING_TENURE)) { tsp_tabu_tenure = atoi(argv[++i]); }
        else if (!strcmp(argv[i], TSP_PARSING_TENURE_A)) { tsp_tabu_tenure_a = atoi(argv[++i]); }
        else if (!strcmp(argv[i], TSP_PARSING_TENURE_F)) { tsp_tabu_tenure_f = atof(argv[++i]); }
        else if (!strcmp(argv[i], TSP_PARSING_VERBOSE)) { tsp_verbose = atoi(argv[++i]); }

        else { printf("Error parsing %s from the command line arguments; use %s to view the command line options.", argv[i], TSP_PARSING_HELP); exit(1); }
    }

    if (tsp_seed > 0)   //if the seed is not at 0 (default value), then a seed has been specified -> generate instance randomly
        tsp_gen_random_instance();
    else    //no seed specified: generating instance from filename given
        tsp_gen_instance_from_file();

    tsp_precompute_costs();
    tsp_precompute_sort_edges();

}

/**
 * @brief Solve the instance based odwdn the type of the algorithm specified, saves the solution and frees the instance
*/
void tsp_solve() {

    int result = 0;
    tsp_init_solution();

    // HEURISTICS

    if (!strcmp(tsp_alg_type, "greedy")) result = tsp_solve_greedy(NULL);                             	//greedy
    else if(!strcmp(tsp_alg_type, "g2opt")) result = tsp_solve_greedy(tsp_find_2opt_swap);            	//greedy with 2opt (first swap policy)
    else if(!strcmp(tsp_alg_type, "g2opt-best")) result = tsp_solve_greedy(tsp_find_2opt_best_swap);  	//greedy with 2opt (best  swap policy)
    else if(!strcmp(tsp_alg_type, "tabu")) result = tsp_solve_tabu();                                   //tabu
    else if(!strcmp(tsp_alg_type, "vns")) result = tsp_solve_vns();                                     //vns
    else if(!strcmp(tsp_alg_type, "fvns")) result = tsp_solve_fvns();                                   //fvns

    // CPLEX ALGORITHMS (base cplex without SECs, Benders's loop)
    
    else if (!strncmp(tsp_alg_type, "cplex", 5)) result = tsp_cplex_solve();
    
    else {
        printf("Error choosing the algorithm to use.");
        exit(1);
    }

    tsp_total_time = tsp_time_elapsed();

    if (result == -1)
        tsp_over_time = 1;

}

int main(int argc, const char** argv) {

    signal(SIGINT, signal_callback_handler);

    tsp_init_defs();

    tsp_parse_cmd(argv, argc);

    if (tsp_verbose >= 0) tsp_instance_info();

    tsp_solve();   //algorithm to find optimal(ish) solutions
    
    if (tsp_verbose >= 0) tsp_print_solution();
    tsp_save_solution();
    tsp_plot_solution();
    
    tsp_free_instance();   //frees the dinamically allocated memory and other finishing operations

}
