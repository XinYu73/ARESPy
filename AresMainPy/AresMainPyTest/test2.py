import sys

sys.path.insert(0, '/work/home/xinyu/workplace/PhdProgram/AresMainPy/src4wrap')
import AresMainPy as Ares
from aresCalculator import aresCalculator
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111, add_adsorbate
from ase.io.trajectory import Trajectory
from ase.spacegroup.symmetrize import FixSymmetry, check_symmetry
from ase.optimize import BFGS, LBFGS
import ase.io

atoms = ase.io.read("ares.cell", format='vasp')
print(atoms.cell)
cal = aresCalculator(inputfile="ares.in",atoms=atoms)
atoms.calc = cal
atoms.set_constraint(FixSymmetry(atoms))
opt = LBFGS(atoms, trajectory='opt.traj', use_line_search=False)
opt.run(fmax=0.01)