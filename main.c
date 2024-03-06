#include "tsp.h"

void tsp_init_defs(tsp_instance* inst);
void tsp_parse_cmd(const int argc, const char** argv, tsp_instance* inst);
void tsp_instance_info(const tsp_instance* inst);

void tsp_solve(tsp_instance* inst);
void tsp_solve_greedy(tsp_instance* inst);
void tsp_solve_g2opt(tsp_instance* inst);

void tsp_cmd_solution(const tsp_instance* inst);
void tsp_save_solution(const tsp_instance* inst);
void tsp_plot_solution(const tsp_instance* inst);

void tsp_gen_random_instance(tsp_instance* inst);
void tsp_gen_instance_from_file(tsp_instance* inst);
int tsp_process_file_line(char* line, tsp_instance* inst, int code);
void tsp_process_node_line(char* line, tsp_instance* inst);
int tsp_process_header_line(const char* line, tsp_instance* inst);

void tsp_help();

int main(int argc, const char** argv) {

    tsp_instance inst;
    tsp_init_defs(&inst);

    tsp_parse_cmd(argc, argv, &inst);
    if (tsp_verbose >= 0) tsp_instance_info(&inst);  //prints instance info

    tsp_solve(&inst);   //algorithm to find optimal(ish) solutions

    tsp_cmd_solution(&inst); 
    tsp_save_solution(&inst);
    tsp_plot_solution(&inst);
    
    tsp_free_instance(&inst);   //frees the dinamically allocated memory

}

void tsp_init_defs(tsp_instance* inst) {  //default values

    tsp_seed = TSP_DEF_SEED;
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

}


void tsp_solve(tsp_instance* inst) {  //solve the instance based on the type of the algorithm specified

    tsp_init_solution(inst);

    if (!strcmp(tsp_alg_type, "greedy")) tsp_solve_greedy(inst);
    else if(!strcmp(tsp_alg_type, "g2opt")) tsp_solve_g2opt(inst);
    
    else {
        printf("Error choosing the algorithm to use.");
        exit(1);
    }

}

void tsp_solve_greedy(tsp_instance* inst) {   //solve using greedy algorithm

    // TODO

    /* procedura:
    1-prendere nodo iniziale
    2-iniziare ciclo in cui:
        a-cerco nodi uscenti
        b-prendo minimo
        c-aggiungo arco a soluzione
        d-mi sposto su altro nodo toccato da arco 
    3-proseguo finchè non torno a nodo iniziale
    4-ritorno
    */

    int i, j, current_node, next_node;
    double best_cost, current_cost;
    int* sol;

    current_node = 0;
    sol = inst->tsp_best_solution;
    sol[0] = current_node;
    for (i=1; i<inst->nnodes; i++) {
        best_cost = DBL_MAX;
        for (j=0; next_node<inst->nnodes; next_node++) {
            if (j==current_cost) continue;
            current_cost = inst->costs[current_node][j];
            if (current_cost<best_cost) {
                best_cost = current_cost;
                next_node = j;
            }
        }
        sol[i] = next_node;
    }


}

void tsp_solve_g2opt(tsp_instance* inst) {    //solve using greedy + 2opt algorithm

    // TODO

}

void tsp_save_solution(const tsp_instance* inst) {  //save the best solution found in a file
    //TODO
}

void tsp_plot_solution(const tsp_instance* inst) {  //plot the best solution found
    //TODO
}


void tsp_gen_random_instance(tsp_instance* inst) {  //generates a random instance

    tsp_allocate_coords_space(inst);    //allocate the memory dinamically using the number of nodes saved in the instance
    
    for (int i = 0; i < inst -> nnodes; i++) {
        inst -> coords[i].x = tsp_rnd_coord();
        inst -> coords[i].y = tsp_rnd_coord();
    }

}

void tsp_gen_instance_from_file(tsp_instance* inst) {   //generates an instance from a TSP file

    FILE* fp;
    char line[200];

    fp = fopen(tsp_file_name, "r");
    if (fp == NULL) {   //If the file specified doesn't exist print the error message
        printf("Error reading the file.");
        exit(1);
    }

    int code = 0;   //used to understand what line I'm working in: 0 - I expect an header, 1 - I expect a node
    while (fgets(line, sizeof(line), fp) != NULL && code >= 0)  //read each line of the file
        code = tsp_process_file_line(line, inst, code); //process the line using the code as context to what I'm reading

}

int tsp_process_file_line(char* line, tsp_instance* inst, int code) {   //process a line from the TSP file

    if (code == 1)  //code == 1 -> I'm expecting a node

        if (isdigit(line[0])) { //checking if it's a node

            tsp_process_node_line(line, inst);  //process the node
            return 1;   //expecting another node next

        } else  //if it's not a node then I've read all nodes and I should process it as an header
            code = 0;

    if (code == 0)  //code == 0 -> I'm reading an header
        return tsp_process_header_line(line, inst); //process the header

}

void tsp_process_node_line(char* line, tsp_instance* inst) {    //process a node

    int coord[3];   //format of the node: <index> <x coord> <y coord>
    int counter = 0, old_c = 0, len=strlen(line);

    for (int i = 0; i < 3; i++) {   //read the 3 numbers
        for (counter; line[counter] != ' ' && counter < len; counter++);    //split using spaces
        if (counter < len) line[counter] = 0;

        coord[i] = atoi(line + old_c);  //read the number
        
        old_c = ++counter;

    }

    inst -> coords[coord[0]-1].x = coord[1];    //saving the node info
    inst -> coords[coord[0]-1].y = coord[2];

}

int tsp_process_header_line(const char* line, tsp_instance* inst) { //process a header (return: 0 if next line is a header, 1 if next line is expected to be a node, -1 if I'm done reading the file)

    if (!strncmp(line, "DIMENSION", strlen("DIMENSION"))) { //save the dimension and dinamically allocate the space
        
        inst -> nnodes = atoi(line+strlen("DIMENSION : "));
        tsp_allocate_coords_space(inst);

        return 0;
    
    }
    if (!strncmp(line, "EDGE_WEIGHT_TYPE", strlen("EDGE_WEIGHT_TYPE"))) {   //check the edge weight type

        if (strncmp(line+strlen("EDGE_WEIGHT_TYPE : "), tsp_edge_weight_type, strlen(tsp_edge_weight_type))) { printf("Unexpected weight type. ATT is the expected one."); exit(1); }

        return 0;

    }
    if (!strncmp(line, "NODE_COORD_SECTION", strlen("NODE_COORD_SECTION"))) return 1;   //next lines should be nodes
    if (!strncmp(line, "EOF", strlen("EOF"))) return -1;    //found the end of the line

    return 0;
    
}


void tsp_help() {   //instructions to use the program

    printf("Use:\n");
    printf("%s <file_name> : to specify a file to obtain the TPS values from.\n", TSP_FILE_P);
    printf("%s <int> : specify the time limit in seconds. (default: %ds)\n", TSP_TIME_LIMIT, (int)TSP_DEF_TL);
    printf("%s <int> : (default type of instance) specify the seed to use to create random TPS data (the seed 0 cannot be used due to implementation choices, default: %d).\n", TSP_SEED, TSP_DEF_SEED);
    printf("%s <int> : specity the number of nodes in the problem (default: %d).\n", TSP_NNODES, TSP_DEF_NNODES);
    printf("%s <str> : Type of algorithm to use ([greedy, g2opt]), (default: greedy).\n", TSP_ALGORITHM);
    printf("%s : Logging option (quiet).\n", TSP_QUIET);
    printf("%s : Logging option (verbose).\n", TSP_VERBOSE);

    exit(0);

}