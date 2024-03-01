#include <gnuplot_c.h>

int main(int argc, const char** argv) {
    h_GPC_Plot* plotter;
    plotter = gpc_init_xy("Best solution", "X coord", "Y coord", GPC_AUTO_SCALE, GPC_KEY_DISABLE);
    fprintf(plotter->pipe, "help");
}