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
    tsp_env.time_limit = 0; //signal all running events to stop

    printf("\n---------------------------------------------");
    printf("\n-  Waiting for all processes to terminate.  -");
    printf("\n---------------------------------------------\n\n");

    tsp_wait_all_threads(); //wait for all process to stop properly

    printf("\n\n---------------------------------------------");
    printf("\n-      All processes ended peacefully.      -");
    printf("\n---------------------------------------------\n\n");

    tsp_env.time_total = time_elapsed();
    tsp_env.status = 2;
    
    if (tsp_verbose >= 0) tsp_print_solution();
    tsp_save_solution();
    tsp_plot_solution();
    
    tsp_free_all(); //frees the dinamically allocated memory and other finishing operations
    
    exit(0);

}

/**
 * @brief Instructions to use the program
*/
void tsp_help() {

    //TODO
    /*
    printf("Use:\n");
    printf("%s <file_name> : to specify a file to obtain the TPS values from.\n", TSP_PARSING_FILE);
    printf("%s <int> : specify the seed to use to create random TPS data (the seed 0 cannot be used due to implementation choices).\n", TSP_PARSING_SEED);
    printf("%s <int> : specity the number of nodes in the problem.\n", TSP_PARSING_NNODES);
    printf("%s <str> : Type of algorithm to use.\n");
    printf("%s <int> : Verbose parameter to control the amount of output.\n", TSP_PARSING_VERBOSE);
    printf("%s <int> : Tenure for the tabu algorithm.\n", TSP_PARSING_TENURE, TSP_DEF_TABU_TENURE);
    printf("%s <int> : Amplitude parameter for the dinamic tenure.\n", TSP_PARSING_TENURE_A, tsp_tabu_tenure_a);
    printf("%s <double> : Frequency parameter for the dinamic tenure (default: %10.4f).\n", TSP_PARSING_TENURE_F, tsp_tabu_tenure_f);
    printf("%s <int> : specify the time limit in seconds (default: %4ds).\n", TSP_PARSING_TIME_LIMIT, (int)TSP_DEF_TL);
    */

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

            if (check == 1) raise_error("Cannot parse both a seed and a file_name.\n");
            
            strcpy(tsp_env.file_name, argv[++i]);
            check = 0;
            
        }
        else if (!strcmp(argv[i], TSP_PARSING_SEED)) {

            if (check == 0) raise_error("Cannot parse both a seed and a file_name.\n");

            tsp_env.seed = atoi(argv[++i]);
            srand(tsp_env.seed);
            init_rand();
            check = 1;
            
        }
        else if (!strcmp(argv[i], TSP_PARSING_TIME_LIMIT)) { tsp_env.time_limit = atof(argv[++i]); }
        else if (!strcmp(argv[i], TSP_PARSING_NNODES)) {

            if (check == 0) raise_error("Cannot parse the number of nodes if a file_name is specified.\n");
            
            tsp_inst.nnodes = atoi(argv[++i]);
            check = 1;
            
        }
        else if (!strcmp(argv[i], TSP_PARSING_HELP)) { tsp_help(); }
        else if (!strcmp(argv[i], TSP_PARSING_ALGORITHM)) {
        
            strcpy(tsp_env.alg_type, argv[++i]);
            if (!strcmp(tsp_env.alg_type, TSP_PARSING_G2OPT)) tsp_env.g2opt_swap_pol = 1;
        
        }
        else if (!strcmp(argv[i], TSP_PARSING_VERBOSE)) { tsp_verbose = atoi(argv[++i]); }

        else if (!strcmp(argv[i], TSP_PARSING_TMP_CHOICE)) { tsp_env.tmp_choice = atoi(argv[++i]); }

        else if (!strcmp(argv[i], TSP_PARSING_BEST_SWAP)) { tsp_env.g2opt_swap_pol = 2; }
        else if (!strcmp(argv[i], TSP_PARSING_TENURE)) { tsp_env.tabu_tenure = atoi(argv[++i]); }
        else if (!strcmp(argv[i], TSP_PARSING_TENURE_A)) { tsp_env.tabu_tenure_a = atoi(argv[++i]); }
        else if (!strcmp(argv[i], TSP_PARSING_TENURE_F)) { tsp_env.tabu_tenure_f = atof(argv[++i]); }
        else if (!strcmp(argv[i], TSP_PARSING_FVNS)) { tsp_env.vns_fvns = 1; }
        else if (!strcmp(argv[i], TSP_PARSING_MIPSTART)) { tsp_env.cplex_mipstart = 1; }
        else if (!strcmp(argv[i], TSP_PARSING_CPLEX_BENDERS)) { tsp_env.cplex_benders = 1; }
        else if (!strcmp(argv[i], TSP_PARSING_CPLEX_PATCHING)) { tsp_env.cplex_patching = 1; }
        else if (!strcmp(argv[i], TSP_PARSING_CPLEX_CANDIDATE)) { tsp_env.cplex_can_cb = 1; }
        else if (!strcmp(argv[i], TSP_PARSING_RELAX_CALLBACK)) { tsp_env.cplex_rel_cb = 1; }

        else raise_error("Error parsing %s from the command line arguments; use %s to view the command line options.", argv[i], TSP_PARSING_HELP);
    
    }

}

/**
 * @brief Solve the instance based odwdn the type of the algorithm specified, saves the solution and frees the instance
*/
void tsp_solve() {

    tsp_init_solution();

    //user info
    if (tsp_verbose >= 0) tsp_instance_info();

    //use algorithm selected
    switch(tsp_find_alg(tsp_env.alg_type)) {
        case 0:     // greedy
            tsp_solve_greedy();
            break;
        case 1:     // g2opt
            tsp_solve_g2opt();
            break;
        case 2:     // tabu
            tsp_solve_tabu();
            break;
        case 3:     // vns
            tsp_solve_vns();
            break;
        case 4:     // cplex
            tsp_solve_cplex();
            break;
        default:    // algorithm not found
            raise_error("Error choosing the algorithm to use.");
    }

    tsp_env.time_total = time_elapsed();
    
    if (tsp_verbose >= 0) tsp_print_solution();
    tsp_save_solution();
    tsp_plot_solution();

}

int main(int argc, const char** argv) {

    signal(SIGINT, signal_callback_handler);
    tsp_init_defs();

    tsp_parse_cmd(argv, argc);
    
    tsp_solve();

    tsp_free_all();

}