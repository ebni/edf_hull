# Minimal EDF deadlines by convex hull

This repository contains the C code to prune the unnecessary deadlines of an EDF scheduled real-time task. The theory supporting the implemented method is described into a paper submitted to the [RTCSA19](https://rtcsa2019.github.io/index/). If interested in getting a copy, please send me an email.

## Quick and dirty run

If you are lucky and want to run it first then understand it later, type the following commands:

- `git clone https://github.com/ebni/edf_hull.git`
- `cd edf_hull`
- `make run`

Want to understand more? Just open the Makefile and go top-down. Or read next.

## Long story

This code detects the minimal set of constraints which are needed to
guarantee EDF schedulability. This action is made by the following steps:

1. the task set parameters are read from stdin
1. for each constraint, a vector is created
1. the set of constraints is transformed
1. the convex hull of the constraints is computed
1. the vertices of the convex hull are printed

### Reading the task parameters

The parameters of the tasks are read from `stdin` in the following order:

1. the number of tasks (integer)
1. for task period, deadline, and offset are read (double). If tasks have no offset, then it must be specified `0` as third task parameter.

Notice that the parameters are read as floating point numbers. Hence,
the least common multiple among periods (often called hyperperiod in
the literature) is computed with a tolerance: two floating point
numbers differing by less than a compile-time macro `EPS` are declared
equal. More details in the function `double lcm_double(double, double)`

### Creating the constraints

For each absolute deadline, the number of jobs per task within such a
deadline are computed, as requested by EDF exact schedulability
condition [[Baruah et al. 1990]](https://doi.org/10.1109/REAL.1990.128746). These numbers of jobs are treated as vectors to be processed.

Also the constraints of non-negative execution time are added.

### Translating the constraints

To achieve the minimal number of constraints, the constraints are
translated. A detailed explanation of this phase is given in the paper.

### Computing the convex hull

The convex hull of the vectors representing the constraints is computed by the [Qhull](https://github.com/qhull/qhull) library, available on github. The constraints are written to the file `data_in.txt` which is taken as input by the executable `qhull`. `qhull` then writes the set of the vertices of the convex hull in the file `data_out.txt`.

The format of the vectors written in `data_in.txt` is as follows:

1. size _n_ of the vectors (number of tasks in our case)
1. number _m_ of vectors
1. _m_ rows of space-separated vectors of size _n_

The format of the set of vertices written in `data_out.txt` is:

1. number _v_ of vertices of the convex hull
1. sequence of _v_ indices (one per row) of the vectors which are vertices

### Printing the result

Finally, the program displays what is the tightest set of constraints that needs to be checked, by reading the indices in `data_out.txt`.
