#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void create_sample_solution(char* sample_file, int nodes_number, int lower_rand, int upper_rand);
void plot_solution(char* solution_file);

int main(int argc, const char** argv) {
    char* solution_file;
    int sample_nodes_number, lower_rand, upper_rand;

    solution_file = "sample_solution.txt";
    remove(solution_file);
    printf("Creating sample solution file. Enter the desired number of nodes: ");
    scanf("%u", &sample_nodes_number);
    printf("Enter range of desired node coordinates as two integers ordered in increasing order. ");
    scanf("%u %u", &lower_rand, &upper_rand);
    create_sample_solution(solution_file, sample_nodes_number, lower_rand, upper_rand);
    printf("Sample solution file created.");

    plot_solution(solution_file);
}

void create_sample_solution(char* sample_file, int nodes_number, int lower_rand, int upper_rand) {
    int rand_x, rand_y, real_upper, first_x, first_y;
    FILE *file;

    real_upper = upper_rand-lower_rand;
    file = fopen(sample_file, "a");
    for (int i=0; i<nodes_number; i++) {
        rand_x = (rand() % real_upper) + lower_rand;
        rand_y = (rand() % real_upper) + lower_rand;
        if (i==0) {
            first_x = rand_x;
            first_y = rand_y;
        }
        fprintf(file, "%d %d %d\n", i, rand_x, rand_y);
    }
    fprintf(file, "%d %d %d\n", nodes_number, first_x, first_y);
    fclose(file);
}

void plot_solution(char* solution_file) {
    int remove_success;
    FILE *file, *command_file;

    file = fopen(solution_file, "r");
    command_file = fopen("command.txt", "a");

    fprintf(command_file, "x=0.; y=0.\n");
    fprintf(command_file, "plot '%s' u (x=$2):(y=$3) w lp\n", solution_file);
    fclose(command_file);

    system("gnuplot -persistent command.txt");
    remove_success = -1;
    while (remove_success!=0) {
        remove_success = remove("command.txt");
    }
}