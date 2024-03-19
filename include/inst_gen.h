#ifndef _GEN_H
#define _GEN_H

#include "tsp.h"

// GENERATING RANDOM INSTANCE
void    tsp_gen_random_instance(tsp_instance* inst);

// GENERATING INSTANCE FROM FILE
void    tsp_gen_instance_from_file(tsp_instance* inst);

#endif