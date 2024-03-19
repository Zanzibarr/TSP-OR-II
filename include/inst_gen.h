#ifndef _GEN_H
#define _GEN_H

#include "tsp.h"

// GENERATING RANDOM INSTANCE

/**
 * Generates a random instance
 *
 * @param inst Where to save the instance
 */
void tsp_gen_random_instance(tsp_instance *inst);

// GENERATING INSTANCE FROM FILE

/**
 * Generates an instance from a TSP file
 *
 * @param inst Where to save the instance
 */
void tsp_gen_instance_from_file(tsp_instance *inst);

#endif