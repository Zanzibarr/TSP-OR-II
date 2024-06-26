#include "../include/inst_gen.h"

/**
 * @brief Process a line from a tsp file as a pair of coordinates
 * 
 * @param line The line to process
*/
void tsp_process_node_line(const char* line) {

    double coord[3];   //format of the node: <index> <x coord> <y coord>
    int counter = 0, old_c = 0, len=strlen(line);
    char line_t[strlen(line)];

    for (int i = 0; i < strlen(line); i++) line_t[i] = line[i];

    for (int i = 0; i < 3; i++) {   //read the 3 numbers

        for (; line_t[counter] != ' ' && counter < len; counter++);    //split using spaces
        if (counter < len) line_t[counter] = 0;   //if I wasn't at the end of the line, split the string using a terminator

        coord[i] = atof(line_t + old_c);  //read the number
        
        old_c = ++counter;

    }

    tsp_inst.coords[(int)coord[0]-1].x = coord[1];    //saving the node info
    tsp_inst.coords[(int)coord[0]-1].y = coord[2];

}

/**
 * @brief Process a line from a tsp file as a header
 * 
 * @param line The line to process
 * 
 * @return -1 : reached EOF, 0 : next line is a header, 1 : next line is SUPPOSED to be a node
*/
int tsp_process_header_line(const char* line) {

    if (!strncmp(line, "DIMENSION", strlen("DIMENSION"))) {
        
        tsp_inst.nnodes = atoi(line+strlen("DIMENSION : "));
        tsp_allocate_coords_space();

        return 0;
    
    }
    if (!strncmp(line, "EDGE_WEIGHT_TYPE", strlen("EDGE_WEIGHT_TYPE"))) {

        if (strncmp(line+strlen("EDGE_WEIGHT_TYPE : "), TSP_EDGE_W_TYPE, strlen(TSP_EDGE_W_TYPE))) raise_error("Error in tsp_process_header_line: unexpected weight type. %s is the expected one.", TSP_EDGE_W_TYPE);

        return 0;

    }
    if (!strncmp(line, "NODE_COORD_SECTION", strlen("NODE_COORD_SECTION"))) return 1;   //next lines should be nodes
    if (!strncmp(line, "EOF", strlen("EOF"))) return -1;

    return 0;
    
}

/**
 * @brief Process a line from a tsp file
 * 
 * @param line The line to process
 * @param code The expected line: 0 : reading a header, 1 : expecting a node
 * 
 * @return The code expected for the next line
*/
int tsp_process_file_line(const char* line, int code) {

    if (code == 1)  { //code == 1 -> I'm expecting a node

        if (isdigit(line[0])) { //checking if it's a node

            tsp_process_node_line(line);  //process the node
            return 1;   //expecting another node next

        }

    }

    //if it's not a node then I've read all nodes and I should process it as an header   ---   code == 0 -> I'm reading an header
    return tsp_process_header_line(line); //process the header

}

// GENERATING RANDOM INSTANCE
void tsp_gen_random_instance() {

    tsp_allocate_coords_space();
    
    for (int i = 0; i < tsp_inst.nnodes; i++) {
        tsp_inst.coords[i].x = rnd_coord();
        tsp_inst.coords[i].y = rnd_coord();
    }

}

// GENERATING INSTANCE FROM FILE
void tsp_gen_instance_from_file() {

    FILE* fp;
    char line[200], relative_file_name[120];
    int code = 0;

    sprintf(relative_file_name, "%s/%s", TSP_INST_FOLDER, tsp_env.file_name);    //where to read the file from

    fp = fopen(relative_file_name, "r");

    if (fp == NULL) raise_error("Error in tsp_gen_instance_from_file: reading the file used to generate the instance.");

    while (fgets(line, sizeof(line), fp) != NULL && code >= 0)
        code = tsp_process_file_line(line, code);   //code used to understand what line I'm working in: 0 - I expect an header, 1 - I expect a node

}
