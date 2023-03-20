import quippy
from quippy import descriptors

import ase
import ase.build
from ase.build import bulk
from ase.visualize import view

import matplotlib.pyplot as plt
import numpy as np

###################
# SOAP DESCRIPTOR #
###################

# n_Z=Number of Different atoms; Z=NUCLEAR charges of the atoms (this case N_z =1 atoms of carbon with Z=6. 
desc = descriptors.Descriptor("soap l_max=6 n_max=12 cutoff=5.0 atom_sigma=0.5")

# These are vectors containing the environmet of the atom up to the cutoff
si = bulk('Si', cubic=True) * 3
D1=desc.calc(si)['data']


si_fcc = bulk('Si', 'fcc', a=5.429) * 3
si_fcc.rattle(0.05)

D2= desc.calc(si_fcc)['data']

D = np.r_[D1, D2]
labels = np.array([1] * len(si) + [2] * len(si_fcc))

#This si the kernel
K = np.power(D @ D.T, 2.0)
