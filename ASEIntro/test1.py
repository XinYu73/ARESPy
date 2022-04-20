from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111, add_adsorbate

h = 1.85
d = 1.10

slab = fcc111('Cu', size=(4, 4, 2), vacuum=10.0) #生成Cu的slab

slab.calc = EMT() # 定slab的calculator
e_slab = slab.get_potential_energy() # slab 对象还有calc方法，神奇

molecule = Atoms('2N', positions=[(0., 0., 0.), (0., 0., d)]) 
molecule.calc = EMT()
e_N2 = molecule.get_potential_energy() 

add_adsorbate(slab, molecule, h, 'ontop') # 建模，铜板上放N
constraint = FixAtoms(mask=[a.symbol != 'N' for a in slab])
slab.set_constraint(constraint)
dyn = QuasiNewton(slab, trajectory='N2Cu.traj')
dyn.run(fmax=0.01)

print('Adsorption energy:', e_slab + e_N2 - slab.get_potential_energy())