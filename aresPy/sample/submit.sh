#!/bin/sh
CODEPATH=/home/xq/Workspace/ARES_PBC/src
mpirun -n 4 $CODEPATH/ares.pbc  #> ares.$1 2>&1
