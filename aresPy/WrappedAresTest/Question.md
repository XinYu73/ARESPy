# Question

```python
import sys

sys.path.insert(0, '/work/home/xinyu/workplace/PhdProgram/aresPy/srcWrap/preprocessed')

import arespy
from mpi4py import MPI
comm = MPI.COMM_WORLD
arespy.Smpi_Math_Module.smpi_init()
rank = comm.Get_rank()
print("comm:    ",comm)
print("rank:    ",rank)
arespy.Smpi_Math_Module.start_time('[Total Time]',True)
arespy.Smpi_Math_Module.start_time('[SCF Time]',True)
arespy.Read_Module.read_file('ares.in')
arespy.Scalapack_Module.init_scala()
arespy.Begin_Module.initial_grid_pbc()
arespy.Potential_Module.vlpp()
arespy.Scf_Module.electronicscf()
# arespy.Smpi_Math_Module.end_time('[SCF Time]',True)
# arespy.Smpi_Math_Module.end_time('[Total Time]',True)
# arespy.Smpi_Math_Module.write_time('[Total Time]',True)
a = arespy.Constants()

print(a.dp)
```

```bash
comm:     <mpi4py.MPI.Intracomm object at 0x2b1ac07e4ef0>
rank:     0
Task name >>> OPTIMUS_PRIME                 
[The Atomic Spherical Cutting Radius is (Angs)]  -4.00000000000000     
Ecut=   940.075301025635      eV
The order of finite difference is:           8
Spin unpolarized
Max simulate steps is:           0
Max simulate steps is:           0
Chebyshev filter is used,order:          12
first Chebyshev filter order:         -24
Chebyshev:first step diag free
RayleighRitz needn''t OrthNorm step
Initialize the subspace by Pseudo orbitals
initial density by program
[The Adding Charge]  0.000000000000000E+000
Use Simple Mixing + (r)Pulay Mixing + kerker
parallel dims           1           1
[Reading PseudoPP FILE]C.pbe-mt-cpi.UPF              
[Reading PseudoPP FILE]H.pbe-mt-cpi.UPF              
[Ion and Charge Num.]           8   8.00000000000000     
[Eigen States Num.]          19
XC: PBE-GGA (Pewdew-Burke-Ernzerhof 1996)
communicate all grids
Total Charge# in Sphere is   7.99999999999999     
Non-Local pseudopotential has been used
R-space-GRIDS:          50          50          50
K-space-GRIDS:           1           1           1
Num of K-used:           1
>----SCF Iterations for Solving KS Equations----<
CheFSI: Dim. of Pseudo Subspace:           8
CheFSI: # of random samping states:          11
Segmentation fault
```

问题可能产生的位置:

1. 源代码preprocess出错
    1. 用preprocessed代码从新编译ares.pbd
    2. 用mpirun -n 2测试, 结果显示没问题
    3. 说明preprocessed的代码没问题
2. f90wrap的代码有问题  主要是duplicate define
