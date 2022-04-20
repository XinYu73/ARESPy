import sys

sys.path.insert(0, '/work/home/xinyu/workplace/PhdProgram/AresMainPy/src4wrap')
import AresMainPy as Ares
from ase.visualize import view
import ase.io
from matplotlib import pyplot as plt
import numpy as np
from aresCalculator import aresCalculator
from ase import Atoms

energylog = np.zeros(85)
atoms = ase.io.read("optMPI.traj", 0)
cal = aresCalculator(inputfile="ares.in", atoms=atoms)
for index in range(85):
    atoms = ase.io.read("optMPI.traj", index)
    energylog[index] = cal.get_potential_energy(atoms=atoms)
np.savetxt("energylog.csv", energylog)
