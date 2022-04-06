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