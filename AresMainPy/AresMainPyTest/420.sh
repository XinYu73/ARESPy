#!/bin/sh 
#SBATCH  --job-name=vasp_job  
#SBATCH  --output=log.out.%j
#SBATCH  --error=log.err.%j
#SBATCH  --partition=debug
#SBATCH  --nodes=1
#SBATCH  --ntasks=4
#SBATCH  --ntasks-per-node=4
#SBATCH  --cpus-per-task=1
#SBATCH  --exclusive

source /work/env/intel2018

srun hostname | sort | uniq >> /tmp/nodefile.$$
NP=`srun hostname | wc -l`

#./calypso.x > caly.log 2>&1 
mpiexec -n 4 python test2.py >mpiexecOut
#mpirun -genv I_MPI_FABRICS shm:ofa -machinefile /tmp/nodefile.$$ -n $NP /work/software/vasp.6.1.0/vasp_std  > vasp.log 2>&1  

rm -rf /tmp/nodefile.$$