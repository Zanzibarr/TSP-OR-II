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

            if (check == 1) tsp_raise_error("Cannot parse both a seed and a file_name.\n");
            
            strcpy(tsp_file_name, argv[++i]);
            check = 0;
            
        }
        else if (!strcmp(argv[i], TSP_PARSING_SEED)) {

            if (check == 0) tsp_raise_error("Cannot parse both a seed and a file_name.\n");

            tsp_seed = atoi(argv[++i]);
            srand(tsp_seed);
            tsp_init_rand();
            check = 1;
            
        }
        else if (!strcmp(argv[i], TSP_PARSING_TIME_LIMIT)) { tsp_time_limit = atof(argv[++i]); }
        else if (!strcmp(argv[i], TSP_PARSING_NNODES)) {

            if (check == 0) tsp_raise_error("Cannot parse the number of nodes if a file_name is specified.\n");
            
            tsp_inst.nnodes = atoi(argv[++i]);
            check = 1;
            
        }
        else if (!strcmp(argv[i], TSP_PARSING_HELP)) { tsp_help(); }
        else if (!strcmp(argv[i], TSP_PARSING_ALGORITHM)) { strcpy(tsp_alg_type, argv[++i]); }
        else if (!strcmp(argv[i], TSP_PARSING_TENURE)) { tsp_tabu_tenure = atoi(argv[++i]); }
        else if (!strcmp(argv[i], TSP_PARSING_TENURE_A)) { tsp_tabu_tenure_a = atoi(argv[++i]); }
        else if (!strcmp(argv[i], TSP_PARSING_TENURE_F)) { tsp_tabu_tenure_f = atof(argv[++i]); }
        else if (!strcmp(argv[i], TSP_PARSING_VERBOSE)) { tsp_verbose = atoi(argv[++i]); }

        else tsp_raise_error("Error parsing %s from the command line arguments; use %s to view the command line options.", argv[i], TSP_PARSING_HELP);
    }

}

/**
 * @brief Solve the instance based odwdn the type of the algorithm specified, saves the solution and frees the instance
*/
void tsp_solve() {

    int result = 0;
    tsp_init_solution();

    //user info
    if (tsp_verbose >= 0) tsp_instance_info();

    //use algorithm selected
    switch(tsp_find_alg(tsp_alg_type)) {
        case 0:     // greedy
            result = tsp_solve_greedy(NULL);
            break;
        case 1:     // greedy with 2opt (first swap policy)
            result = tsp_solve_greedy(tsp_find_2opt_swap);
            break;
        case 2:     // greedy with 2opt (best  swap policy)
            result = tsp_solve_greedy(tsp_find_2opt_best_swap);
            break;
        case 3:     // tabu
            result = tsp_solve_tabu();
            break;
        case 4:     // vns
            result = tsp_solve_vns();
            break;
        case 5:     // fvns
            result = tsp_solve_fvns();
            break;
        case 6:     // cplex-base
            result = tsp_solve_cplex();
            break;
        case 7:     // cplex-benders
            result = tsp_solve_cplex();
            break;
        case 8:     // cplex-benders-patching
            result = tsp_solve_cplex();
            break;
        case 9:
            result = tsp_solve_cplex_bnc();
            break;
        default:    // algorithm not found
            tsp_raise_error("Error choosing the algorithm to use.");
    }

    tsp_total_time = tsp_time_elapsed();

    if (result) tsp_over_time = result;
    
    if (tsp_verbose >= 0) tsp_print_solution();
    tsp_save_solution();
    tsp_plot_solution();

}

int main(int argc, const char** argv) {

    signal(SIGINT, signal_callback_handler);
    tsp_init_defs();

    tsp_parse_cmd(argv, argc);
    tsp_solve();

    tsp_free_instance();

}