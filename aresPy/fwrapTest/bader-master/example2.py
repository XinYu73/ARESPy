import numpy as np

import ase.units as units

from ase.io import write
from gpaw import restart

from bader import *

atoms, calc = restart('si8.gpw')
rho = calc.get_pseudo_density()

write('si8.cube', atoms, data=rho * units.Bohr**3)

bdr, ions, chg = bader(atoms, rho, full_output=True, refine_edge_itrs=-3)

bader_mod.bader_output(bdr, ions, chg)

w = np.zeros(list(chg.npts) + [bdr.bnum], order='F')
charge_mod.get_weights(chg, w)
