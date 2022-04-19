import sys

sys.path.insert(0, '/work/home/xinyu/workplace/PhdProgram/AresMainPy/src4wrap')
import AresMainPy as Ares



def aa():
    Ares.Smpi_Math_Module().smpi_init()
    Ares.Read_Module().read_file('/work/home/xinyu/workplace/PhdProgram/AresMainPy/AresMainPyTest/ares.in')
    print("fff")
    Ares.Scalapack_Module().init_scala()
    print("fff")
    Ares.Relax_Module().initialize_relax()
    print("fff")
    aresOut = Ares.Aresmainapi().aresOut()
    print("fff")
    Ares.Aresmainapi().init_alloc_arrays(aresOut, Ares.Struct_Module().natom)
    print("fff")
    Ares.Forcestress_Module().cal_force_stress()
    Ares.Aresmainapi().assignment(aresOut)

aa()