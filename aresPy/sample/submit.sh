#!/bin/sh 
#SBATCH  --job-name=aresjob  
#SBATCH  --output=log.out.%j
#SBATCH  --error=log.err.%j
#SBATCH  --partition=debug
#SBATCH  --nodes=1
#SBATCH  --ntasks=2
#SBATCH  --ntasks-per-node=2
#SBATCH  --cpus-per-task=1
#SBATCH  --exclusive

source /work/env/intel2018

srun hostname | sort | uniq >> /tmp/nodefile.$$
NP=`srun hostname | wc -l`

CODEPATH=/work/home/xinyu/workplace/PhdProgram/aresPy/src
mpirun -n 1 $CODEPATH/ares.pbc  > ares.$1 2>&1
