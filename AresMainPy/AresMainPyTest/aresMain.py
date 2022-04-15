import sys
sys.path.insert(0, '/work/home/xinyu/workplace/PhdProgram/AresMainPy/src4wrap')
import AresMainPy as Ares
import numpy as np
from mpi4py import MPI
from matplotlib import pyplot as plt
plt.style.use(['science','ieee','no-latex'])

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
AresApi.init_alloc_arrays(aresOut,30)
Forcestress_Module.cal_force_stress()
AresApi.assignment(aresOut)
print(">>>>",aresOut.comm)
# #!###################
# posStore=np.zeros((3,30,Parameters.nssp))

# #!################
# for i in range(Parameters.nssp):
#     print(">>>>Current Ionic Step",i)
#     Forcestress_Module.cal_force_stress()
#     AresApi.assignment(aresOut)
#     posStore[:,:,i]=aresOut.pos
#     Relax_Module.relaxer(i)
# #!###############
# for i in range(Parameters.nssp):
#     print(">>>>pos at ",i,"\n")
#     print(posStore[:,:,i])
# #!Plot>>>>>>>>>>>>>

# plt.figure()
# plt.plot(posStore[0,:,2]-posStore[0,:,1])
# plt.plot(posStore[1,:,2]-posStore[1,:,1])
# plt.plot(posStore[2,:,2]-posStore[2,:,1])
# plt.legend([r"$dim1$",r"$dim2$",r"$dim3$"])
# plt.savefig('./Figure/ionStep3_2.png')

# plt.figure()
# plt.plot(posStore[0,:,5]-posStore[0,:,0])
# plt.plot(posStore[1,:,5]-posStore[1,:,0])
# plt.plot(posStore[2,:,5]-posStore[2,:,0])
# plt.legend([r"$dim1$",r"$dim2$",r"$dim3$"])
# plt.savefig('./Figure/ionStep2_1.png')

print(">>>>>>>>>>>>>>>End Wrapped AresMain<<<<<<<<<<<<<<<")