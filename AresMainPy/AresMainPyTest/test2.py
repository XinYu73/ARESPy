import sys

sys.path.insert(0, '/work/home/xinyu/workplace/PhdProgram/AresMainPy/src4wrap')
import AresMainPy as Ares
from aresCalculator import aresCalculator
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111, add_adsorbate

cal = aresCalculator()
atoms = Atoms('2N', positions=[(0., 0., 0.), (0., 0., 1)])
atoms.calc = cal
print(atoms.get_potential_energy())
print(atoms.get_forces())
print(atoms.get_stress())
