# TSP-OR-II
University project for the OR2 course

## How to build
Open terminal in this folder and write:
```shell
make
```
>I use >make && clean

Then to run:
```shell
./tsp <options>
```

### Command line options 
 - "-file" to specify a file to obtain the TPS values from
 - "-seed" specify the seed to use to create random TPS data (the seed 0 cannot be used due to implementation choices)
 - "-nodes" specity the number of nodes in the problem
 - "-alg" Type of algorithm to use ([greedy, g2opt, g2opt_best, tabu, vns, fvns, base_cplex, benders])
 - "-tenure" Tenure for the tabu algorithm
 - "-tenure-a" Amplitude parameter for the dinamic tenure
 - "-tenure-f" Frequency parameter for the dinamic tenure
 - "-tl" specify the time limit in seconds

## Test run file
The user needs to list in the textual file specified by the -test option what algorithms they want to use in the test run, along with their hyperparameters if needed.
The file must be saved in the perf_prof/test_spec directory. However, the user does not to specify the whole filepath, but only the name file;
Here is a list of the algorithms, along with their hyperparameters and how they need to be expressed in the file:
 - Greedy algorithm: "greedy", no hyperparameters;
 - Greedy + 2-opt algorithm: "g2opt" for the first swap policy, "g2opt-best" for the best swap policy, no hyperparameters in both cases;
 - Tabu search: "tabu -t <tenure value> -a <tenure amplitude value> -f <tenure frequency value>": -t, -a, -f specify the value for the dinamic tenure formula, respectively the base value, the amplitude and the frequency;
    > they can be written in any order;
    > both amplitude and frequency have 0 as default values, used if they are not specified;
    > the base value must be specified and cannot be 0;
    > there can be at most one value for all three hyperparameters;
 - VNS: "vns", no hyperparameters.

If the user wishes to specify the seed to be used to generate the test instances, the first row of the file must read "-seed <seed value>". If the first row does not read this, the seed will be generated by the program. 0 is not a valid value for the seed: it will be treated as the sign of an error and cause the program to stop.


## Verbose options
The user needs to specify the TSP_VERBOSE value in the utils.h file before compiling each time they wish to change it.
Here is a list of the values:
 - <0 for quiet                             (nothing)
 - [0, 10[ for normal                       (basic info for final user)
 - == 5 for thread info                     (multithreading)
 - >=10 for new best solutions              (visual info)
 - >=50 to plot intermediate costs          (plotting)
 - >=100 for integrity checks               (integrity checks enabled) (suggested while in development)
 - >=500 to see the path in the solution    (advanced debugging)
 - >=1000 for super-verbose                 (full verbose)
