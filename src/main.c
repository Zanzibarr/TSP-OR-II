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
    
    #if TSP_VERBOSE >= 0
    tsp_print_solution();
    #endif
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
            if (!strcmp(argv[++i], "benders"))
                snprintf(tsp_alg_type, 20, "cplex_%s", argv[i]);
            else strcpy(tsp_alg_type, argv[i]);
        }
        else if (!strcmp(argv[i], TSP_PARSING_TENURE)) { tsp_tabu_tenure = atoi(argv[++i]); }
        else if (!strcmp(argv[i], TSP_PARSING_TENURE_A)) { tsp_tabu_tenure_a = atoi(argv[++i]); }
        else if (!strcmp(argv[i], TSP_PARSING_TENURE_F)) { tsp_tabu_tenure_f = atof(argv[++i]); }
        else { printf("Error parsing the command line arguments; use %s to view the command line options.", TSP_PARSING_HELP); exit(1); }
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

/**
 * @brief Performs a test run using the textual file specified in the cmd commands.
 * TODO Tratta il main come black box per fare i test, crea uno script esterno che chiama il main in loop in modo che anche il debugging sia più facile...
void tsp_perform_test_run() {

    // open the specified file for the test run
    char test_run_file_path[500];
    snprintf(test_run_file_path, 500, "%s/%s", TSP_TEST_RUN_FILES_FOLDER, tsp_test_run_file);
    FILE* test_run_file = fopen(test_run_file_path, "r");
    if (test_run_file == NULL) {
        printf("Error opening the file specified for the test run");
        exit(1);
    }

    // read the first line and look for a seed
    char line[100], seed_command[5], seed[50];
    fgets(line, 150, test_run_file);
    strncpy(seed_command, line, 4);
    if (!strcmp(seed_command, "seed")) {
        strncpy(seed, line+5, 50);
        tsp_seed = atoi(seed);
        if (tsp_seed==0) {
            printf("Invalid seed detected.\n");
            exit(1);
        }
        srand(tsp_seed);
        tsp_init_rand();
    }
    else {
        printf("No specified seed; it will be generated by the program: ");
        tsp_seed = (uint64_t) time(NULL);
        printf("%d.\n", tsp_seed);
        fclose(test_run_file);
        test_run_file = fopen(test_run_file_path, "r");
    }

    // read all lines and save the algorithms to use
    char **test_algorithms;
    tsp_test_hyperparameters* test_hyperparameters;
    int alg_number = 0;
    while (fgets(line, 150, test_run_file)) {

        test_algorithms = (char**) realloc(test_algorithms, sizeof(char*)*(++alg_number));
        test_hyperparameters = (tsp_test_hyperparameters*)
            realloc(test_hyperparameters, sizeof(tsp_test_hyperparameters)*(alg_number));
        char algorithm[50], hyperparameters[100];
        int last_pos;
        sscanf(line, "%s %n", algorithm, &last_pos);
        strcpy(hyperparameters, line+last_pos);
        hyperparameters[strlen(hyperparameters)-1]='\0';

        char alg_present = 0;
        for (int i=0; i<tsp_algorithms_number; i++) {
            if (!strcmp(algorithm, tsp_algorithms[i])) { alg_present = 1; break; }
        }
        if (!alg_present) {
            printf("Non-valid algorithm found: %s.\n", algorithm);
            exit(1);
        }
        test_algorithms[alg_number-1] = (char*) calloc(strlen(algorithm), 1);
        strcpy(test_algorithms[alg_number-1], algorithm);

        tsp_test_hyperparameters hyper = {0,0,0};
        if (!strcmp(algorithm, "tabu")) {
            char command[3];
            int string_pos = 0, hyper_length = strlen(hyperparameters);
            char t_read = 0, a_read = 0, f_read = 0, c;
            while (string_pos<hyper_length) {
                c = *(hyperparameters+string_pos);
                command[0] = c;
                command[1] = *(hyperparameters+string_pos+1);
                command[2] = '\0';
                if (!strcmp(command, "-t")) t_read++;
                else {
                    if (!strcmp(command, "-a")) a_read++;
                    else {
                        if (!strcmp(command, "-f")) f_read++;
                        else {
                            printf("Unrecognized specifier detected: %s.\n", command);
                            exit(1);
                        }
                    }
                }
                if (t_read>1 || a_read>1 || f_read>1) {
                    printf("Command %s detected more than once.\n", command);
                    exit(1);
                }
                string_pos += 3;
                int start_string_pos = string_pos, i=0;
                char point_read=0, value[50];
                while (*(hyperparameters+string_pos)!=' ' && *(hyperparameters+string_pos)) {
                    c = *(hyperparameters+string_pos);
                    char condition = c>='0' && c<='9';
                    if (!strcmp(command, "-f")) {
                        char condition2 = c=='.' && !point_read;
                        condition = condition || condition2;
                        if (condition2) point_read = 1;
                    }
                    if (!condition) {
                        printf("Invalid character detected after %s command.\n", command);
                        exit(1);
                    }
                    value[i++] = c;
                    string_pos++;
                }
                value[i] = '\0';
                if (!strcmp(command, "-t")) hyper.tabu_tenure = atoi(value);
                else if (!strcmp(command, "-a")) hyper.tabu_tenure_a = atoi(value);
                else if (!strcmp(command, "-f")) hyper.tabu_tenure_f = strtod(value, '\0');
                string_pos++;
            }
            if (t_read<1) {
                printf("No base tenure value detected.\n");
                exit(1);
            }
        }

        if (!strcmp(algorithm, "tabu")) test_hyperparameters[alg_number-1] = hyper;
        //else test_hyperparameters[alg_number-1] = NULL;


    }

    /*for (int i=0; i<alg_number; i++) {
        printf("Algoritmo numero %d: %s\n", i+1, test_algorithms[i]);
        if (test_hyperparameters[i]!=NULL) printf("Iperparametri per algoritmo numero %d: %s\n", i+1, test_hyperparameters[i]);
    }

    fclose(test_run_file);

    // generate one test instance at a time and solve it with every test algorithm

    char test_results_file_path[500];
    tsp_test_run_file[strlen(tsp_test_run_file)-4]='\0';
    snprintf(test_results_file_path, 500, "%s/%s_txtfile_results.txt", TSP_TEST_RUN_RESULTS_FOLDER, tsp_test_run_file);
    remove(test_results_file_path);
    FILE* test_results_file = fopen(test_results_file_path, "a");
    fprintf(test_results_file, "%d, ", alg_number);
    for (int i=0; i<alg_number; i++) {
        tsp_test_hyperparameters* h = &test_hyperparameters[i];
        fprintf(test_results_file, "%s", test_algorithms[i]);
        if (!strcmp(test_algorithms[i], "tabu")) {
            fprintf(test_results_file, "-t%d-a%d-f%f", h->tabu_tenure, h->tabu_tenure_a, h->tabu_tenure_f);
        }
        if (i!=alg_number-1) fprintf(test_results_file, ", ");
    }
    fprintf(test_results_file, "\n");

    tsp_inst.nnodes = TSP_TEST_NNODES;
    tsp_time_limit = TSP_TEST_TL;
    for (int i_inst=0; i_inst<TSP_TEST_NINST; i_inst++) {

        tsp_gen_random_instance();
        tsp_precompute_costs();
        tsp_precompute_sort_edges();

        fprintf(test_results_file, "rndinst_%d_%d, ", i_inst, tsp_seed);

        for (int i_alg=0; i_alg<alg_number; i_alg++) {

            tsp_init_solution();
            tsp_set_initial_time();

            printf("Starting algorithm %d: %s...\n", i_alg+1, test_algorithms[i_alg]);
            if (!strcmp(test_algorithms[i_alg], "greedy")) tsp_solve_greedy(NULL);
            else if(!strcmp(test_algorithms[i_alg], "g2opt")) tsp_solve_greedy(tsp_find_2opt_swap);
            else if(!strcmp(test_algorithms[i_alg], "g2opt-best")) tsp_solve_greedy(tsp_find_2opt_best_swap);
            else if(!strcmp(test_algorithms[i_alg], "tabu")) {
                //rndinst_<indice>_<seed>
                tsp_tabu_tenure = test_hyperparameters[i_alg].tabu_tenure;
                tsp_tabu_tenure_a = test_hyperparameters[i_alg].tabu_tenure_a;
                tsp_tabu_tenure_f = test_hyperparameters[i_alg].tabu_tenure_f;
                //printf("Hyperparameters for this tabu: tenure %d, amplitude %d, frequency %f\n", tsp_tabu_tenure, tsp_tabu_tenure_a, tsp_tabu_tenure_f);
                tsp_solve_tabu();
            }
            else if(!strcmp(test_algorithms[i_alg], "vns")) tsp_solve_vns();

            printf("Algorithm %d found its best solution in %fs.\n", i_alg+1, tsp_inst.best_time);
            fprintf(test_results_file, "%f", tsp_inst.best_time);
            if (i_alg!=alg_number-1) fprintf(test_results_file, ", ");

        }

        if (i_inst!=TSP_TEST_NINST-1) fprintf(test_results_file, "\n");

        tsp_free_instance();

    }

    fclose(test_results_file);

    // run the Python script to create the perfomance profile
    char test_plot_file_path[500];
    snprintf(test_plot_file_path, 500, "%s/%s_txtfile_plot.pdf", TSP_TEST_RUN_RESULTS_FOLDER, tsp_test_run_file);
    remove(test_plot_file_path);
    char python_command[100], script_path[1000];
    snprintf(python_command, 1000, "python perf_prof/perfprof.py -D , -T %f -S 2 -M 20 %s %s -P 'first plot' ", tsp_time_limit, test_results_file_path, test_plot_file_path);
    system(python_command);

}
*/

int main(int argc, const char** argv) {

    signal(SIGINT, signal_callback_handler);

    tsp_init_defs();

    tsp_parse_cmd(argv, argc);

    #if TSP_VERBOSE >= 0
    tsp_instance_info();
    #endif

    tsp_solve();   //algorithm to find optimal(ish) solutions
    
    #if TSP_VERBOSE >= 0
    tsp_print_solution();
    #endif
    tsp_save_solution();
    tsp_plot_solution();
    
    tsp_free_instance();   //frees the dinamically allocated memory and other finishing operations

}
