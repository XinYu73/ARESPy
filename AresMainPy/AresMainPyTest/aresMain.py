import sys
sys.path.insert(0, '/work/home/xinyu/workplace/PhdProgram/AresMainPy/src4wrap')
import AresMainPy as Ares
import numpy as np
from mpi4py import MPI

print(">>>>>>>>>>>>>>>Test Wrapped AresMain<<<<<<<<<<<<<<<")

Smpi_Math_Module=Ares.Smpi_Math_Module()
Read_Module=Ares.Read_Module()
AresApi = Ares.Aresmainapi()
Scalapack_Module = Ares.Scalapack_Module()
Relax_Module = Ares.Relax_Module()
Parameters = Ares.Parameters()
Forcestress_Module = Ares.Forcestress_Module()
aresOut = AresApi.aresOut()
#!###################
Smpi_Math_Module.smpi_init()
Read_Module.read_file('ares.in')
Scalapack_Module.init_scala()
Relax_Module.initialize_relax()
print(Parameters.nssp)
for i in range(Parameters.nssp-2):
    print(">>>>Current Ionic Step",i)
    Forcestress_Module.cal_force_stress()
    AresApi.init_alloc_arrays(aresOut,30)
print(aresOut.forces)
print(">>>>>>>>>>>>>>>End Wrapped AresMain<<<<<<<<<<<<<<<")