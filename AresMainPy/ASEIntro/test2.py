from ase.md.verlet import VelocityVerlet
from ase import units
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111, add_adsorbate

h = 1.85
d = 1.10

slab = fcc111('Cu', size=(4, 4, 2), vacuum=10.0)  #生成Cu的slab

slab.calc = EMT()  # 定slab的calculator
e_slab = slab.get_potential_energy()  # slab 对象还有calc方法，神奇

molecule = Atoms('2N', positions=[(0., 0., 0.), (0., 0., d)])
molecule.calc = EMT()
dyn = VelocityVerlet(molecule, dt=1.0 * units.fs)
for i in range(10):
    pot = molecule.get_potential_energy()
    kin = molecule.get_kinetic_energy()
    print('%2d: %.5f eV, %.5f eV, %.5f eV' % (i, pot + kin, pot, kin))
    dyn.run(steps=20)
