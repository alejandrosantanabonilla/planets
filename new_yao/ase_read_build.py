from ase.io import read, write
from ase.visualize import view

# How to easily read a POSCAR structure
graph=read('structure.poscar')

# How to get the lattice constants from a poscar file
# as a numpy.array

latt_const=graph.get_cell()

# View the structure you have read.
#view(a)

# Positions of the Graphene sheet
positions_graph=graph.get_positions()

# Read cluster positions
cluster=read('cluster.xyz')

clu_positions=cluster.get_positions()
