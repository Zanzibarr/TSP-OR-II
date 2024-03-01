#include "tsp.h"

void parse_cmd(const int argc, const char** argv, instance* inst);

int main(int argc, const char** argv) {

    if (argc < 3) {
        printf("Too few arguments. You must at least specify the file or the seed to obtain the data.");
        exit(1);
    }

    instance inst;
    initialize_defaults();

    parse_cmd(argc, argv, &inst);

    plot_rnd(&inst);

}

void initialize_defaults() {

    time_limit = DEF_TL;
    seed = -1;

}

void parse_cmd(const int argc, const char** argv, instance* inst) {

    char file_name[1000];

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], FILE) == 0) { strcpy(file_name, argv[++i]); continue; }
        if (strcmp(argv[i], SEED) == 0) { seed = atoi(argv[++i]); continue; }
        if (strcmp(argv[i], TIME_LIMIT) == 0) { time_limit = atof(argv[++i]); continue; }
    }

    if (seed >= 0) {
        gen_random_instance(seed, &inst);
    } else {
        gen_instance_from_file(file_name, &inst);
    }

}