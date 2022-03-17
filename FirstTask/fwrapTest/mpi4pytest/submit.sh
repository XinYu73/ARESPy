#!/bin/sh 
#SBATCH  --job-name=xinyujob  
#SBATCH  --output=log.out.%j
#SBATCH  --error=log.err.%j
#SBATCH  --partition=debug
#SBATCH  --nodes=1
#SBATCH  --ntasks=4
#SBATCH  --ntasks-per-node=4
#SBATCH  --cpus-per-task=1
##SBATCH  --exclusive

source /work/env/intel2018
source /work/software/miniconda2/bin/activate /work/home/xinyu/soft/XinYuEnv

srun hostname | sort | uniq >> /tmp/nodefile.$$
NP=`srun hostname | wc -l`


#ARESPATH=/work/home/xinyu/soft/ARES/ARES_PBC/src
#mpirun -genv I_MPI_FABRICS shm:tcp -machinefile /tmp/nodefile.$$ -n $NP $ARESPATH/ares.pbc  > aresTest.log 2>&1

mpiexec -n 4 python hello_mpi.py

#CODEPATH=/home/xq/Workspace/ARES_PBC/src
#mpirun -n 4 $CODEPATH/ares.pbc  #> ares.$1 2>&1
