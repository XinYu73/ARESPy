import sys

sys.path.insert(0, '/work/home/xinyu/workplace/PhdProgram/AresMainPy/src4wrap')
import AresMainPy as Ares
import ase.io
import numpy as np
atoms = ase.io.read("ares.cell", format='vasp')
print(atoms.positions.T)
Ares.Smpi_Math_Module().smpi_init()
Ares.Read_Module().read_file("ares.in")
Ares.Aresmainapi().updateions(atoms.positions.T, atoms.cell)
Ares.Scalapack_Module().init_scala()
Ares.Relax_Module().initialize_relax()
aresOut = Ares.Aresmainapi().aresOut()
Ares.Aresmainapi().init_alloc_arrays(aresOut, Ares.Struct_Module().natom)
Ares.Forcestress_Module().cal_force_stress()
Ares.Aresmainapi().assignment(aresOut)
print(aresOut.pos)
print(aresOut.apilat_mat)