#include "tsp.h"

void tsp_init_defs(tsp_instance* inst);
void tsp_parse_cmd(const int argc, const char** argv, tsp_instance* inst);
void tsp_help();
void tsp_instance_info(const tsp_instance* inst);
void tsp_gen_random_instance(tsp_instance* inst);
void tsp_gen_instance_from_file(tsp_instance* inst);
int tsp_process_file_line(char* line, tsp_instance* inst, int code);
void tsp_process_node_line(char* line, tsp_instance* inst);
int tsp_process_header_line(const char* line, tsp_instance* inst);
void tsp_save_solution(const tsp_instance* inst);
void tsp_plot_solution(const tsp_instance* inst);
void tsp_free_instance(tsp_instance* inst);

int main(int argc, const char** argv) {

    if (argc < 2) { //basics check for the correct use of the program
        printf("Too few arguments. You must at least specify the file or the seed to obtain the data.\nUse %s to get help for command line use.", TSP_HELP);
        exit(1);
    }

    tsp_instance inst;  //generate the instance
    tsp_init_defs(&inst); //initialize the instance

    tsp_parse_cmd(argc, argv, &inst);   //parse the command line arguments
    if (tsp_verbose > 0) tsp_instance_info(&inst);  //prints the info of the instance and the problem's parameters

    // SOLVE
    tsp_save_solution(&inst);   //save the current solution into a file
    tsp_plot_solution(&inst);   //plot the solution using gnuplot
    
    tsp_free_instance(&inst);   //frees the dinamically allocated memory

}

void tsp_init_defs(tsp_instance* inst) {  //default values

    tsp_seed = 0;
    tsp_time_limit = TSP_DEF_TL;
    tsp_verbose = 0;

    strcpy(tsp_file_name, "NONE");
    strcpy(tsp_edge_weight_type, "ATT");

    inst -> nnodes = TSP_DEF_NNODES;

}

void tsp_parse_cmd(const int argc, const char** argv, tsp_instance* inst) { //parse the command line arguments to prepare the instance and the problem's parameters

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], TSP_FILE_P)) { strcpy(tsp_file_name, argv[++i]); }
        else if (!strcmp(argv[i], TSP_SEED)) { tsp_seed = atoi(argv[++i]); srand(tsp_seed); tsp_init_rnd(); }
        else if (!strcmp(argv[i], TSP_TIME_LIMIT)) { tsp_time_limit = atoi(argv[++i]); }
        else if (!strcmp(argv[i], TSP_NNODES)) { inst -> nnodes = atoi(argv[++i]); }
        else if (!strcmp(argv[i], TSP_HELP)) { tsp_help(); }
        else if (!strcmp(argv[i], TSP_QUIET)) { tsp_verbose = -1; }
        else if (!strcmp(argv[i], TSP_VERBOSE)) { tsp_verbose = 1; }
        else { printf("Error parsing the command line arguments; use %s to view the command line options.", TSP_HELP); exit(1); }
    }

    if (tsp_seed > 0)   //if the seed is not at 0 (default value), then a seed has been specified -> generate instance randomly
        tsp_gen_random_instance(inst);
    else    //no seed specified: generating instance from filename given
        tsp_gen_instance_from_file(inst);

}

void tsp_help() {   //instructions to use the program

    printf("Use:\n");
    printf("%s <file_name> : to specify a file to obtain the TPS values from\n", TSP_FILE_P);
    printf("%s <int> : specify the time limit in seconds\n", TSP_TIME_LIMIT);
    printf("%s <int> : specify the seed to use to create random TPS data (the seed 0 cannot be used due to implementative choices)", TSP_SEED);
    printf("\n%s <int> : specity the number of nodes in the problem (needed only for the random generated grid)", TSP_NNODES);

    exit(0);

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

}

void tsp_gen_random_instance(tsp_instance* inst) {  //generates a random instance
    
    inst -> coords = calloc(inst -> nnodes, sizeof(tsp_pair));  //allocate the memory dinamically using the number of nodes saved in the instance

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
        inst -> coords = calloc(inst -> nnodes, sizeof(tsp_pair));
        
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

void tsp_save_solution(const tsp_instance* inst) {
    //TODO
}

void tsp_plot_solution(const tsp_instance* inst) {
    //TODO
}

void tsp_free_instance(tsp_instance* inst) {    //frees the dinamically allocated memory

    free (inst -> coords);

}