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

    //TODO: Can I tell cplex to stop and retrieve it's solution?
    
    if (tsp_verbose >= 0) tsp_print_solution();
    int unique = tsp_save_solution();
    //tsp_plot_solution(unique);
    
    tsp_free_all(); //frees the dinamically allocated memory and other finishing operations
    
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
        else if (!strcmp(argv[i], TSP_PARSING_CPLEX_CANDIDATE)) { tsp_env.cplex_can_cb = 1; }
        else if (!strcmp(argv[i], TSP_PARSING_RELAX_CALLBACK)) { tsp_env.cplex_rel_cb = 1; }
        else if (!strcmp(argv[i], TSP_PARSING_CPLEX_PATCHING)) { tsp_env.cplex_patching = 1; }
        else if (!strcmp(argv[i], TSP_PARSING_CPLEX_GREEDY_PATCHING)) { tsp_env.cplex_patching = 2; }

        else raise_error("Error parsing %s from the command line arguments. See the README.md to get instructions.\n", argv[i]);
    
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
    int unique = tsp_save_solution();
    //tsp_plot_solution(unique);

}

int main(int argc, const char** argv) {

    signal(SIGINT, signal_callback_handler);
    tsp_init_defs();

    tsp_parse_cmd(argv, argc);
    
    tsp_solve();

    tsp_free_all();

}