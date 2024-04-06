# TSP-OR-II
University project for the OR2 course

## How to build
Open terminal in this folder.
In Linux, write:
```shell
make
```
or
```shell
make && clean
```
In Windows (mingw32 compiler), write:
```shell
mingw32-make -f Makefile_win
```

Then to run:
```shell
./tsp <options>
```

## Command line options 
 - "-file" to specify a file to obtain the TPS values from
 - "-seed" specify the seed to use to create random TPS data (the seed 0 cannot be used due to implementation choices)
 - "-nodes" specity the number of nodes in the problem
 - "-alg" Type of algorithm to use ([greedy, g2opt, g2opt_best, tabu, vns, fvns, base_cplex, benders])
 - "-tenure" Tenure for the tabu algorithm
 - "-tenure-a" Amplitude parameter for the dinamic tenure
 - "-tenure-f" Frequency parameter for the dinamic tenure
 - "-tl" specify the time limit in seconds

## Verbose options
The user needs to specify the TSP_VERBOSE value in the utils.h file before compiling each time they wish to change it.
Here is a list of the values:
 - $<0$ for quiet                               (nothing)
 - $[0, 10[$ for normal                         (basic info for final user)
 - $5$ for thread info                          (multithreading)
 - $[10, 50[$ for new best solutions            (visual info)
 - $[50, 100[$ to plot intermediate costs       (plotting)
 - $[100, 500[$ for integrity checks            (integrity checks enabled) (suggested while in development)
 - $[500, 100[$ to see the path in the solution (advanced debugging)
 - $\geq1000$ for super-verbose                 (full verbose)
