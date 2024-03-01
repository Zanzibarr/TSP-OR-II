#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

// PARSING
#define FILE "-file"
#define TIME_LIMIT "-tl"
#define SEED "-seed"

// TIME MANAGEMENT
#define MS_SEC 1e+3
#define DEF_TL 8.64e+7

// NUMBERS
#define EPSYLON 1e-9

uint64_t time_limit;
uint64_t seed;

double random_value() {
    
}

#endif