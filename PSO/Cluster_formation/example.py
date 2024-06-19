from element_creation import *
from create_spherical_points import *
import numpy
from pyscf import gto, scf, dft

#Get characters without shuffling
combinations = [("In",5), ("P",5)]
split_chars = generate_and_split_characters(combinations,shuffle=False)

#Generate 3D points
center = np.array([0, 0, 0])  # Sphere center
radius = 10.0
min_dist = 3.0  # Minimum distance between points
num_points = 10  # Desired number of points
points = tuple(sample_poisson_disk_sphere(center, radius, min_dist, num_points))

# Convert the points from a list to a numpy array for easier slicing
points = np.array(points)
result = np.column_stack((split_chars, points))

cluster = convert_to_atom_list(result)

#mol = gto.Mole()
#mol.atom = cluster
#mol.basis='def2-svp'
#mol.build()

#mf = dft.ROKS(mol)
#mf.xc='PBE'
#mf.nlc='VV10'
#mf.kernel()

for idx, values in enumerate(cluster):
    print (values[0],values[1][0],values[1][1],values[1][2])
