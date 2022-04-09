import sys

sys.path.insert(0, '/work/home/xinyu/workplace/PhdProgram/aresPy/srcWrap/preprocessed')
import arespy as ares
from mpi4py import MPI
comm = MPI.COMM_WORLD
ares.Smpi_Math_Module.smpi_init()
ares.Read_Module.read_file('ares.in')
# print(ares.Parameters().system_name.decode('utf8'))
# ares.Smpi_Math_Module.smpi_init()

ares.Smpi_Math_Module.start_time('[Total Time]',True)
ares.Smpi_Math_Module.start_time('[SCF Time]',True)
print(ares.Struct_Module.struct_type.pos)
ares.Scalapack_Module.init_scala()
ares.Begin_Module.initial_grid_pbc()
ares.Potential_Module.vlpp()
ares.Scf_Module.electronicscf()
ares.Smpi_Math_Module.end_time('[SCF Time]',True)
out = ares.Output_Module()
out.output()

ares.Smpi_Math_Module.end_time('[Total Time]',True)
ares.Smpi_Math_Module.write_time('[Total Time]',True)
ares.Smpi_Math_Module.smpi_exit()
