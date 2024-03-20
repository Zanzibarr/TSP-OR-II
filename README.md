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
./main <options>
```

### Command line options
 - "-file" to specify a file to obtain the TPS values from
 - "-seed" specify the seed to use to create random TPS data (the seed 0 cannot be used due to implementation choices)
 - "-nodes" specity the number of nodes in the problem
 - "-alg" Type of algorithm to use ([greedy, g2opt, g2opt_best, tabu, vns])
 - "-tenure" Tenure for the tabu algorithm
 - "-tenure-a" Amplitude parameter for the dinamic tenure
 - "-tenure-f" Frequency parameter for the dinamic tenure
 - "-tl" specify the time limit in seconds


## Verbose options
Look at comment in the utils.h file under the TSP_VERBOSE definition.  

You need to specify the TSP_VERBOSE value before compiling each time you wish to change it.
