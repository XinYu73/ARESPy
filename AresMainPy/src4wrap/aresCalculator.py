import sys

sys.path.insert(0, '/work/home/xinyu/workplace/PhdProgram/AresMainPy/src4wrap')
import AresMainPy as Ares
from ase.calculators.calculator import Calculator
import ase.io


class aresCalculator(Calculator):
    """Ares calculator for ase"""
    implemented_properties = ['energy', 'forces', 'stress']

    def __init__(self,
                 atoms=None,
                 inputfile=None,
                 input_data=None,
                 structureFile=None,
                 **kwargs):
        Calculator.__init__(self, atoms=atoms, input_data=input_data, **kwargs)
        self.restart()
        self.atoms = atoms
        self.inputfile = inputfile
        self.structureFile = structureFile
        if self.structureFile is not None:
            self.atoms = ase.io.read(self.structureFile, format='vasp')
        self.driver_initialise()

    def driver_initialise(self):
        Ares.Smpi_Math_Module().smpi_init()
        Ares.Read_Module().read_file(self.inputfile)
        Ares.Scalapack_Module().init_scala()
        Ares.Relax_Module().initialize_relax()
        self.aresOut = Ares.Aresmainapi().aresOut()
        Ares.Aresmainapi().init_alloc_arrays(self.aresOut,
                                             Ares.Struct_Module().natom)

    def restart(self):
        self._energy = None
        self._forces = None
        self._stress = None

    def get_potential_energy(self, atoms=None, **kwargs):
        self.atoms = atoms 
        Ares.Aresmainapi().updateions(self.atoms.positions.T,self.atoms.cell)
        Ares.Forcestress_Module().cal_force_stress()
        Ares.Aresmainapi().assignment(self.aresOut)
        self._energy = Ares.Energy_Module().etot * Ares.Constants().hart2ev
        return self._energy

    def get_forces(self, atoms=None):
        self._forces = self.aresOut.forces.T
        return self._forces

    def get_stress(self, atoms=None):
        self._stress = self.aresOut.stress * Ares.Constants().au2gpa
        return self._stress

    def __del__(self):
        Ares.End_Module().destroy_beast()