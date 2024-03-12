#ifndef _GEN_H
#define _GEN_H

#include "tsp.h"

// GENERATING RANDOM INSTANCE
void    tsp_gen_random_instance(tsp_instance* inst);

// GENERATING INSTANCE FROM FILE
void    tsp_gen_instance_from_file(tsp_instance* inst);

int     tsp_process_file_line(char* line, tsp_instance* inst, int code);
void    tsp_process_node_line(char* line, tsp_instance* inst);
int     tsp_process_header_line(const char* line, tsp_instance* inst);

#endif