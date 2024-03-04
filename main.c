#include "tsp.h"

void initialize_defaults(instance* inst);
void parse_cmd(const int argc, const char** argv, instance* inst);
void help();
void instance_info(const instance* inst);
void gen_random_instance(instance* inst);
void gen_instance_from_file(instance* inst);
int process_file_line(char* line, instance* inst, int code);
void process_node_line(char* line, instance* inst);
int process_header_line(const char* line, instance* inst);
void save_solution(const instance* inst);
void plot_solution(const instance* inst);
void free_instance(instance* inst);

int main(int argc, const char** argv) {

    if (argc < 2) { //basics check for the correct use of the program
        printf("Too few arguments. You must at least specify the file or the seed to obtain the data.\nUse %s to get help for command line use.", HELP);
        exit(1);
    }

    instance inst;  //generate the instance
    initialize_defaults(&inst); //initialize the instance

    parse_cmd(argc, argv, &inst);   //parse the command line arguments
    if (verbose > 0) instance_info(&inst);  //prints the info of the instance and the problem's parameters

    // SOLVE
    save_solution(&inst);   //save the current solution into a file
    plot_solution(&inst);   //plot the solution using gnuplot
    
    free_instance(&inst);   //frees the dinamically allocated memory

}

void initialize_defaults(instance* inst) {  //default values

    seed = 0;
    time_limit = DEF_TL;
    verbose = 0;

    strcpy(file_name, "NONE");
    strcpy(edge_weight_type, "ATT");

    inst -> nnodes = DEF_NNODES;

}

void parse_cmd(const int argc, const char** argv, instance* inst) { //parse the command line arguments to prepare the instance and the problem's parameters

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], FILE_P)) { strcpy(file_name, argv[++i]); continue; }
        if (!strcmp(argv[i], SEED)) { seed = atoi(argv[++i]); srand(seed); continue; }
        if (!strcmp(argv[i], TIME_LIMIT)) { time_limit = atoi(argv[++i]); continue; }
        if (!strcmp(argv[i], NNODES)) { inst -> nnodes = atoi(argv[++i]); continue; }
        if (!strcmp(argv[i], HELP)) { help(); }
        if (!strcmp(argv[i], QUIET)) { verbose = -1; }
        if (!strcmp(argv[i], VERBOSE)) { verbose = 1; }
    }

    if (seed > 0)   //if the seed is not at 0 (default value), then a seed has been specified -> generate instance randomly
        gen_random_instance(inst);
    else    //no seed specified: generating instance from filename given
        gen_instance_from_file(inst);

}

void help() {   //instructions to use the program

    printf("Use:\n");
    printf("%s <file_name> : to specify a file to obtain the TPS values from\n", FILE_P);
    printf("%s <int> : specify the time limit in seconds\n", TIME_LIMIT);
    printf("%s <int> : specify the seed to use to create random TPS data (the seed 0 cannot be used due to implementative choices)", SEED);
    printf("\n%s <int> : specity the number of nodes in the problem (needed only for the random generated grid)", NNODES);

    exit(0);

}

void instance_info(const instance* inst) {  //prints the instance info and problem's parameters

    printf("--------------------\n");
    printf("Type of Instance: %s\n", ((seed == 0) ? "from file" : "random"));
    if (seed == 0) printf("File name: %s\n", file_name);
    else printf("Seed: %ld\n", seed);
    printf("Time limit: %lds\n", time_limit);
    printf("Number of nodes: %d\n", inst -> nnodes);
    printf("Edge weight type: %s\n", edge_weight_type);
    printf("--------------------\n");

}

void gen_random_instance(instance* inst) {  //generates a random instance
    
    inst -> coords = calloc(inst -> nnodes, sizeof(pair));  //allocate the memory dinamically using the number of nodes saved in the instance

    for (int i = 0; i < inst -> nnodes; i++) {
        inst -> coords[i].x = rnd();
        inst -> coords[i].y = rnd();
    }

}

void gen_instance_from_file(instance* inst) {   //generates an instance from a TSP file

    FILE* fp;
    char line[200];

    fp = fopen(file_name, "r");
    if (fp == NULL) {   //If the file specified doesn't exist print the error message
        printf("Error reading the file.");
        exit(1);
    }

    int code = 0;   //used to understand what line I'm working in: 0 - I expect an header, 1 - I expect a node
    while (fgets(line, sizeof(line), fp) != NULL && code >= 0)  //read each line of the file
        code = process_file_line(line, inst, code); //process the line using the code as context to what I'm reading

}

int process_file_line(char* line, instance* inst, int code) {   //process a line from the TSP file

    if (code == 1)  //code == 1 -> I'm expecting a node

        if (isdigit(line[0])) { //checking if it's a node

            process_node_line(line, inst);  //process the node
            return 1;   //expecting another node next

        } else  //if it's not a node then I've read all nodes and I should process it as an header
            code = 0;

    if (code == 0)  //code == 0 -> I'm reading an header
        return process_header_line(line, inst); //process the header

}

void process_node_line(char* line, instance* inst) {    //process a node

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

int process_header_line(const char* line, instance* inst) { //process a header (return: 0 if next line is a header, 1 if next line is expected to be a node, -1 if I'm done reading the file)

    if (!strncmp(line, "NAME", strlen("NAME"))) return 0;   //ignored
    if (!strncmp(line, "TYPE", strlen("TYPE"))) return 0;   //ignored
    if (!strncmp(line, "COMMENT", strlen("COMMENT"))) return 0; //ignored
    if (!strncmp(line, "DIMENSION", strlen("DIMENSION"))) { //save the dimension and dinamically allocate the space
        inst -> nnodes = atoi(line+strlen("DIMENSION : "));
        inst -> coords = calloc(inst -> nnodes, sizeof(pair));
        return 0;
    }
    if (!strncmp(line, "EDGE_WEIGHT_TYPE", strlen("EDGE_WEIGHT_TYPE"))) {   //save the edge weight type
        int k = 0;
        for (int i = strlen("EDGE_WEIGHT_TYPE : "); i < strlen(line) && line[i] != '\n'; i++, k++)
            edge_weight_type[k] = line[i];

        return 0;
    }
    if (!strncmp(line, "EDGE_WEIGHT_FORMAT", strlen("EDGE_WEIGHT_FORMAT"))) {
        //TODO
        return 0;
    }
    if (!strncmp(line, "EDGE_DATA_FORMAT", strlen("EDGE_DATA_FORMAT"))) {
        //TODO
        return 0;
    }
    if (!strncmp(line, "NODE_COORD_TYPE", strlen("NODE_COORD_TYPE"))) {
        //TODO
        return 0;
    }
    if (!strncmp(line, "DISPLAY_DATA_TYPE", strlen("DISPLAY_DATA_TYPE"))) {
        //TODO
        return 0;
    }
    if (!strncmp(line, "EOF", strlen("EOF"))) return -1;    //found the end of the line
    if (!strncmp(line, "NODE_COORD_SECTION", strlen("NODE_COORD_SECTION"))) return 1;   //next lines should be nodes

}

void save_solution(const instance* inst) {
    //TODO
}

void plot_solution(const instance* inst) {
    //TODO
}

void free_instance(instance* inst) {    //frees the dinamically allocated memory

    free (inst -> coords);

}