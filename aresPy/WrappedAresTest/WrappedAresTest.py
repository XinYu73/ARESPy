import sys

sys.path.insert(0, '/work/home/xinyu/workplace/PhdProgram/aresPy/srcWrap/preprocessed')
import arespy as ares
import numpy as np
from mpi4py import MPI
comm = MPI.COMM_WORLD



ares.Smpi_Math_Module.smpi_init()
ares.Smpi_Math_Module.start_time('[Total Time]',True)
ares.Smpi_Math_Module.start_time('[SCF Time]',True)
ares.Read_Module.read_file('ares.in')
aresParameters=ares.Parameters()
print("aresParameters.system_name.decode('utf-8')",aresParameters.system_name.decode('utf-8'))
print("aresParameters.cellfile_name.decode('utf-8')",aresParameters.cellfile_name.decode('utf-8'))
print("ppfile_name",aresParameters.ppfile_name) # !why an array 

print("ecut",aresParameters.ecut)
# print(aresParameters.gridn)
# print(aresParameters.kgrid)
# print(aresParameters.kspacing)
# print(aresParameters.nspin)
# print(aresParameters.etol)
print("ntype",aresParameters.ntype)

ares.Scalapack_Module.init_scala()
ares.Begin_Module.initial_grid_pbc()
ares.Potential_Module.vlpp()
ares.Scf_Module.electronicscf()
hart2ev = ares.Constants().hart2ev
print("Total Energy",ares.Energy_Module().etot*hart2ev)
print("Band Energy",ares.Energy_Module().eband*hart2ev)
print("XC Energy",ares.Energy_Module().exc*hart2ev)
print("Hartree Energy",ares.Energy_Module().eh*hart2ev)
print("I-E Energy",ares.Energy_Module().eext*hart2ev)
print("Ion-Ion Energy",ares.Struct_Module().eionion*hart2ev)
print("Free Energy",ares.Energy_Module().fe*hart2ev)
print("0K Energy",ares.Energy_Module().fe0*hart2ev)
print("natom, naty",ares.Struct_Module().natom,ares.Struct_Module().naty)

# print("forces",struct.forces)

# print("Force",aresstruct_module.struct_type().forces)



# nspin = ares.Parameters().nspin
# ni1 = ares.Grid_Module().global_n1
# ni2 = ares.Grid_Module().global_n2
# ni3 = ares.Grid_Module().global_n3
# nps = ares.Grid_Module().n
# grid = ares.Grid_Module().grid_type()
# print(grid.rhos)
# rho= np.zeros((ni1,ni2,ni3))
# ares.Smpi_Math_Module.smpi_init()
# ares.Read_Module.read_file('ares.in')
# # print(ares.Parameters().system_name.decode('utf8'))
# # ares.Smpi_Math_Module.smpi_init()

# ares.Smpi_Math_Module.start_time('[Total Time]',True)
# ares.Smpi_Math_Module.start_time('[SCF Time]',True)
# print(ares.Struct_Module.struct_type.pos)
# ares.Scalapack_Module.init_scala()
# ares.Begin_Module.initial_grid_pbc()
# ares.Potential_Module.vlpp()
# ares.Scf_Module.electronicscf()
# ares.Smpi_Math_Module.end_time('[SCF Time]',True)
# out = ares.Output_Module()
# out.output()

# ares.Smpi_Math_Module.end_time('[Total Time]',True)
# ares.Smpi_Math_Module.write_time('[Total Time]',True)
# ares.Smpi_Math_Module.smpi_exit()
