#include "tsp.h"

void initialize_defaults(instance* inst);
void parse_cmd(const int argc, const char** argv, instance* inst);
void help();
void gen_random_instance(instance* inst);
void gen_instance_from_file(instance* inst);
int process_file_line(char* line, instance* inst, int code);
void free_instance(instance* inst);

int main(int argc, const char** argv) {

    if (argc < 2) {
        printf("Too few arguments. You must at least specify the file or the seed to obtain the data.\nUse %s to get help for command line use.", HELP);
        exit(1);
    }

    instance inst;
    initialize_defaults(&inst);

    parse_cmd(argc, argv, &inst);

    if (verbose > 0) {
        printf("--------------------\n");
        printf("Type of Instance: %s\n", ((seed == 0) ? "from file" : "random"));
        if (seed == 0) printf("File name: %s\n", file_name);
        else printf("Seed: %ld\n", seed);
        printf("Time limit: %lds\n", time_limit);
        printf("Number of nodes: %d\n", inst.nnodes);
        printf("Edge weight type: %s\n", edge_weight_type);
        printf("--------------------\n");
    }

    free_instance(&inst);

    //plot_rnd(&inst);

}

void initialize_defaults(instance* inst) {

    seed = 0;
    time_limit = DEF_TL;
    verbose = 0;

    strcpy(file_name, "NONE");
    strcpy(edge_weight_type, "ATT");

    inst -> nnodes = DEF_NNODES;

}

void parse_cmd(const int argc, const char** argv, instance* inst) {

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], FILE_P)) { strcpy(file_name, argv[++i]); continue; }
        if (!strcmp(argv[i], SEED)) { seed = atoi(argv[++i]); srand(seed); continue; }
        if (!strcmp(argv[i], TIME_LIMIT)) { time_limit = atoi(argv[++i]); continue; }
        if (!strcmp(argv[i], NNODES)) { inst -> nnodes = atoi(argv[++i]); continue; }
        if (!strcmp(argv[i], HELP)) { help(); }
        if (!strcmp(argv[i], QUIET)) { verbose = -1; }
        if (!strcmp(argv[i], VERBOSE)) { verbose = 1; }
    }

    if (seed > 0)
        gen_random_instance(inst);
    else
        gen_instance_from_file(inst);

}

void help() {

    printf("Use:\n%s <file_name> : to specify a file to obtain the TPS values from\n%s <int> : specify the time limit in seconds\n%s <int> : specify the seed to use to create random TPS data (the seed 0 cannot be used due to implementative choices)\n%s <int> : specity the number of nodes in the problem (needed only for the random generated grid)", FILE_P, TIME_LIMIT, SEED, NNODES);

    exit(0);

}

void gen_random_instance(instance* inst) {
    
    inst -> coords = calloc(inst -> nnodes, sizeof(pair));

    for (int i = 0; i < inst -> nnodes; i++) {
        inst -> coords[i].x = rnd();
        inst -> coords[i].y = rnd();
    }

}

void gen_instance_from_file(instance* inst) {

    FILE* fp;
    int c;
    char line[200];

    fp = fopen(file_name, "r");
    if (fp == NULL) {
        printf("Error reading the file.");
        exit(1);
    }

    int code = 0;
    while (fgets(line, sizeof(line), fp) != NULL && code >= 0)
        code = process_file_line(line, inst, code);

}

int process_file_line(char* line, instance* inst, int code) {

    if (code == 1) {

        if (isdigit(line[0])) {
            
            int coord[3];
            int counter = 0, old_c = 0, len=strlen(line);

            for (int i = 0; i < 3; i++) {
                for (counter; line[counter] != ' ' && counter < len; counter++);
                if (counter < len) line[counter] = 0;

                coord[i] = atoi(line + old_c);
                
                old_c = ++counter;

            }

            inst -> coords[coord[0]-1].x = coord[1];
            inst -> coords[coord[0]-1].y = coord[2];

            return 1;

        } else
            code = 0;

    }

    if (code == 0) {
        
        if (!strncmp(line, "NAME", strlen("NAME"))) return 0;
        if (!strncmp(line, "TYPE", strlen("TYPE"))) return 0;
        if (!strncmp(line, "COMMENT", strlen("COMMENT"))) return 0;
        if (!strncmp(line, "DIMENSION", strlen("DIMENSION"))) {
            inst -> nnodes = atoi(line+strlen("DIMENSION : "));
            inst -> coords = calloc(inst -> nnodes, sizeof(pair));
            return 0;
        }
        if (!strncmp(line, "EDGE_WEIGHT_TYPE", strlen("EDGE_WEIGHT_TYPE"))) {
            int k = 0;
            for (int i = strlen("EDGE_WEIGHT_TYPE : "); i < strlen(line) && line[i] != '\n'; i++, k++) {
                edge_weight_type[k] = line[i];
            }
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
        if (!strncmp(line, "EOF", strlen("EOF"))) return -1;
        if (!strncmp(line, "NODE_COORD_SECTION", strlen("NODE_COORD_SECTION"))) return 1;
    
    }

}

void free_instance(instance* inst) {

    free (inst -> coords);

}