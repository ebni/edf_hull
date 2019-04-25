# Minimal EDF deadlines by convex hull

This repository contains the C code to prune the unnecessary deadlines
of an EDF scheduled real-time task. The theory supporting the
implemented method is described into a paper submitted to the
[RTCSA19](https://rtcsa2019.github.io/index/). If interested in getting a copy, please
send me an email.

## Quick and dirty run

If you are lucky and want to first run it then understand it then type
the following commands:
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
1. the result of the convex

If you are interestend in the theory behind it, please refer to the paper 



If you are 