#include "tsp.h"

void initialize_defaults(instance* inst);
void parse_cmd(const int argc, const char** argv, instance* inst);
void help();
void gen_random_instance(instance* inst);
void gen_instance_from_file(const char* file_name, instance* inst);
void free_instance(instance* inst);

int main(int argc, const char** argv) {

    if (argc < 2) {
        printf("Too few arguments. You must at least specify the file or the seed to obtain the data.");
        exit(1);
    }

    instance inst;
    initialize_defaults(&inst);

    parse_cmd(argc, argv, &inst);

    free_instance(&inst);

    //plot_rnd(&inst);

}

void initialize_defaults(instance* inst) {

    time_limit = DEF_TL;
    seed = 0;

    inst -> nnodes = DEF_NNODES;

}

void parse_cmd(const int argc, const char** argv, instance* inst) {

    char file_name[1000];
    file_name[0] = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], FILE) == 0) { strcpy(file_name, argv[++i]); continue; }
        if (strcmp(argv[i], SEED) == 0) { srand(atoi(argv[++i])); continue; }
        if (strcmp(argv[i], TIME_LIMIT) == 0) { time_limit = atoi(argv[++i]); continue; }
        if (strcmp(argv[i], NNODES) == 0) { inst -> nnodes = atoi(argv[++i]); continue; }
        if (strcmp(argv[i], HELP) == 0) { help(); }
    }

    if (seed > 0) {
        printf("HERE\n");
        gen_random_instance(inst);
    } else {
        gen_instance_from_file(file_name, inst);
    }

}

void help() {

    printf("Use:\n");
    printf(FILE);
    printf(" <file_name> : to specify a file to obtain the TPS values from\n");
    printf(TIME_LIMIT);
    printf(" <int> : specify the time limit in seconds\n");
    printf(SEED);
    printf(" <int> : specify the seed to use to create random TPS data (the seed 0 cannot be used due to implementative choices)\n");
    printf(NNODES);
    printf(" <int> : specity the number of nodes in the problem (needed only for the random generated grid)");

    exit(0);

}

void gen_random_instance(instance* inst) {
    
    inst -> xcoords = malloc(sizeof(double) * inst -> nnodes);
    inst -> ycoords = malloc(sizeof(double) * inst -> nnodes);

    for (int i = 0; i < inst -> nnodes; i++) {
        inst -> xcoords[i] = rnd();
        inst -> ycoords[i] = rnd();
    }

    for (int i = 0; i < inst -> nnodes; i++) {
        printf("(%f;%f)", inst -> xcoords[i], inst -> ycoords[i]);
    }

}

void gen_instance_from_file(const char* file_name, instance* inst) {

    printf(file_name);

}

void free_instance(instance* inst) {

    free (inst -> xcoords);
    free (inst -> ycoords);

}