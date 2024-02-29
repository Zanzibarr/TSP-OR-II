#include "tsp.h"

int main(int argc, char** argv) {

    if (argc < 3) {
        printf("Too few arguments. You must at least specify the file or the seed to obtain the data.");
        exit(1);
    }

    instance inst;

    parse_cmd(argc, argv, &inst);

    plot_rnd(&inst);

}

void parse_cmd(const int argc, const char** argv, instance* inst) {

    char* file_name;
    int seed = -1;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], FILE) == 0) { strcpy(file_name, argv[++i]); continue; }
        if (strcmp(argv[i], SEED) == 0) { seed = atoi(argv[++i]); continue; }
        if (strcmp(argv[i], TIME_LIMIT) == 0) { time_limit = atoi(argv[++i]); continue; }
    }

    if (seed >= 0) {
        gen_random_instance(seed, &inst);
    } else {
        gen_instance_from_file(file_name, &inst);
    }

}

// tsp -file ... -tl ...