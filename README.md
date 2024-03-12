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
- "-file": specify the .txt file from which to read the TSP instance
- "-tl": specify the time limit for the solving algorithm (in seconds)
- "-seed": specify the seed to use for the random number generator
- "-nodes": specify the number of nodes for the randomly generated instance
- "-help": print on the terminal the instructions on how to use the program
- "-alg": specify the algorithm to use to solve the instance

## Verbose options
Look at comment in the utils.h file under the TSP_VERBOSE definition.  

You need to specify the TSP_VERBOSE value before compiling each time you wish to change it.
