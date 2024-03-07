#include "tsp.h"

void tsp_gen_random_instance(tsp_instance* inst);
void tsp_gen_instance_from_file(tsp_instance* inst);
int tsp_process_file_line(char* line, tsp_instance* inst, int code);
void tsp_process_node_line(char* line, tsp_instance* inst);
int tsp_process_header_line(const char* line, tsp_instance* inst);


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
