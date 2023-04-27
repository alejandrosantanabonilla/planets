import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms
from ase import Atoms
from ase.geometry import cell_to_cellpar
from ase.visualize import view
import re

# Calculate the cluster rotation matrix
def cluster_rotation_matrix(angle_x, angle_y, angle_z, clusters):
    """ Function to perform a Geometrical Rotation 
        for an user-provided set of coordinates.

    """

    # Get the centroid of the Atoms object
    center_of_mass = clusters.get_center_of_mass()

    # Define the rotation angle in the x, y, z directions, in radians

    rotation_matrix_x = np.array([[1, 0, 0],
                                  [0, np.cos(angle_x), -np.sin(angle_x)],
                                  [0, np.sin(angle_x), np.cos(angle_x)]])

    rotation_matrix_y = np.array([[np.cos(angle_y), 0, np.sin(angle_y)],
                                  [0, 1, 0],
                                  [-np.sin(angle_y), 0, np.cos(angle_y)]])

    rotation_matrix_z = np.array([[np.cos(angle_z), -np.sin(angle_z), 0],
                                  [np.sin(angle_z), np.cos(angle_z), 0],
                                  [0, 0, 1]])

    rotation_matrix = np.dot(rotation_matrix_z, np.dot(rotation_matrix_y, rotation_matrix_x))

    # Rotate all atoms of the Atoms object around the center of mass and output absolute coordinates
    rotated_positions = np.dot(rotation_matrix, (clusters.positions - center_of_mass).T).T + center_of_mass

    return rotated_positions


# Define the cluster displacement function

def cluster_dis_xyz(x_direct, y_direct, angle_x, angle_y, angle_z, clusters, substrate):
    """ Function to rotate and displace an user-provided cluster
    """
    # Read the lattice parameters of the substrate from the POSCAR file
    

    # Define the direct coordinates of the offset
    direct_coords = np.array([x_direct, y_direct, 0])

    # Convert the direct coordinates to cartesian coordinates
    cart_coords = direct_coords @ substrate.cell

    translation_vector_xy = np.array(cart_coords)

    # Translate the rotated coordinates along the x and y axes
    translated_positions = cluster_rotation_matrix(angle_x, angle_y, angle_z, clusters) + translation_vector_xy

    # Subtract the smallest z-coordinate from the z-coordinates of all atoms,
    #and translate it to a minimum of 2.5 angstrom

    min_z = np.min(translated_positions[:, 2])
    translated_positions[:, 2] -= min_z
    translation_vector_z = np.array([0, 0, z_min])
    translated_positions += translation_vector_z
    
    return translated_positions

def substrate_preprocess(substrate):

    # Selective dynamics information, whether or not it is possible to move the atomic capital owned by the owner
    if not substrate.constraints:
        substrate.set_constraint(True)

    # General Atoms object export to selective dynamics format, Cartesian system
    write('substrate_process.poscar', substrate, vasp5=True, direct=False, sort=False, long_format=True, symbol_count=None)

def combine_substrate_cluster(substrate, rot_cluster):
    
    # Combine the two Atoms objects
    combined_atoms = substrate + rot_cluster

    combined_atoms.wrap()
    # Output the combined Atoms object to a new POSCAR file

    return combined_atoms

def structure_creator(cluster_file, substrate_file, output_name, x_direct, y_direct, angle_x, angle_y, angle_z):
   """ Function to create new structures based on user-defined
        parameters
   """

   # Read the cluster.POSCAR file and convert to Cartesian with selective dynamics style
   clusters = read(str(cluster_file), format='vasp')
   rot_clusters = clusters.copy()
 
   #clusters.set_constraint(FixAtoms(mask=[False]*len(clusters)))
   substrate = read(str(substrate_file), format='vasp')
   new_subs=substrate.copy()
   
   if not new_subs.constraints:
      new_subs.set_constraint(True)
      
   # Update the positions of all atoms in the Atoms object to the rotated coordinates plus the translation vector
   rot_clusters.positions = cluster_dis_xyz(x_direct, y_direct, angle_x, angle_y, angle_z, clusters, substrate)

   final=combine_substrate_cluster(new_subs, rot_clusters)
    
   write(str(output_name), final, vasp5=True, direct=False, sort=False, long_format=True, symbol_count=None)

def read_energy(filename):
   e0_pattern = re.compile(r"energy\(sigma->0\)\s*=\s+([\d\-\.]+)")
   textfile = open(filename, 'r')
   filetext = textfile.read()
   textfile.close()
   matches = re.findall(e0_pattern, filetext)

   return matches

all_energies=read_energy('OUTCAR')

if __name__ == "__main__":
   # Define Variables for rotation
   x_direct = 0.6
   y_direct  = 0.4
   z_min = 2.5
   angle_x = np.pi/2
   angle_y = 0
   angle_z = 0

   structure_creator("cluster.poscar", "substrate.poscar","combined.poscar", x_direct, y_direct, angle_x, angle_y, angle_z)
   final=read("combined.poscar", format='vasp')
   view(final)

   # Reading an OUTCAR and getting the last energy
   all_energies=read_energy('OUTCAR')
   print (all_energies[-1])
