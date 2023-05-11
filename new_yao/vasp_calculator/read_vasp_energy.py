import numpy as np
import re
import os

from ase.io import read, write
from ase.constraints import FixAtoms
from ase import Atoms
from ase.geometry import cell_to_cellpar
from ase.visualize import view

e0_pattern = re.compile(r"energy\(sigma->0\)\s*=\s+([\d\-\.]+)")

textfile = open(filename, 'r')
filetext = textfile.read()
textfile.close()

matches = re.findall(e0_pattern, filetext)

print (matches)
