import random
import numpy as np
from pyscf import gto, scf, dft

import numpy as np

def sum_displacements(A, B):
  """
  This function sums matrix A with matrix B, adding 
  A to the second and third columns of B.

  Args:
      A: A numpy array representing matrix A.
      B: A numpy array representing matrix B.

  Returns:
      A numpy array representing the sum of A and B.
  """

  # Check if matrices have compatible shapes
  if A.shape[0] != B.shape[0]:
    raise ValueError("Matrices must have the same number of rows.")

  # Select the numeric columns from B (assuming the first column contains strings)
  B_numeric = B[:, 1:].astype(float)

  # Add A to the second and third columns of B_numeric
  C = np.column_stack((B[:, 0], A + B_numeric))

  return C

def write_xyz_file(data, energy, filename):
  """
  This function writes a list of atomic data to a .xyz file with four columns.

  Args:
      data: A list of tuples, where each tuple represents an atom (element, x, y, z).
      filename: The name of the output .xyz file (string).
  """
  with open(filename, 'w') as f:
    # Write the number of atoms
    f.write(str(len(data)) + '\n')

    # Optional comment line (uncomment if desired)
    f.write('SCF energy: {} \n'.format(energy))

    # Write each atom data with four columns
    for element, x, y, z in data:
      f.write(f"{element} {x:.10f} {y:.10f} {z:.10f}\n")

def run_b3lyp_calculation(atomic_coordinates):
  """
  This function performs a B3LYP calculation using PySCF with the provided coordinates.

  Args:
      atomic_coordinates: A numpy array of shape (num_atoms, 4) containing atomic coordinates and element.

  Prints the coordinates before the calculation.
  """
  # Convert coordinates to PySCF format

  mol = gto.Mole()
  mol.atom=list(transform_data(atomic_coordinates))
  mol.unit = 'B'
  mol.basis = {'C': "sto-3g"}
  mol.build()
  
  # Perform B3LYP calculation (uncomment to run the calculation)
  mf_hf = dft.RKS(mol)
  mf_hf.xc="b3lyp"
  mf_hf = mf_hf.newton()

  return mf_hf.kernel()


def metropolis_monte_carlo(atomic_coordinates, num_steps, temperature, get_energy, dx):
  """
  This function performs a Metropolis Monte Carlo simulation for energy minimization.

  Args:
      atomic_coordinates: A numpy array of shape (num_atoms, 4) containing atomic coordinates and element.
      num_steps: The number of Monte Carlo steps.
      temperature: The temperature parameter for the simulation (controls acceptance rate).
      get_energy: A function that takes atomic coordinates and returns the energy.
      dx: Parameter to define the length of the perturbation. It ranges from (-0.5,0.5)

  Returns:
      A tuple containing the optimized atomic coordinates and the corresponding energy.
  """
  best_coordinates = atomic_coordinates.copy()  # Copy for tracking best state
  best_energy = run_b3lyp_calculation(best_coordinates)  # Initial energy
  print("MC initial energy:", best_energy)

  count=0
  for _ in range(num_steps):
    count=count+1
    print ("Iteration_{}:".format(count))
    #Creating displacements for X and Y. 
    displacements = np.random.uniform(low=-dx, high=dx, size=atomic_coordinates.shape[:2])
    new_coordinates=sum_displacements(displacements[:, :2], best_coordinates.copy())
    
    # Calculate the energy of the trial configuration
    trial_energy = run_b3lyp_calculation(new_coordinates)

    # Metropolis acceptance criterion
    delta_energy = trial_energy - best_energy
    print("MC delta energy:", delta_energy)
    
    acceptance_prob = np.exp(-delta_energy / temperature)

    # Accept or reject the trial configuration
    if np.random.rand() < acceptance_prob:
      best_coordinates = new_coordinates.copy()
      best_energy = trial_energy

    write_xyz_file(transform_data(best_coordinates), best_energy, ''.join(["mol_",str(count),".xyz"]))

  print ("After {} iterations the MC search has finished".format(count))

def add_to_indices(data, value):
  """
  This function adds a value to indices 1, 2, and 3 of each tuple in a list.

  Args:
      data: A list containing tuples.
      value: The value to be added to the specified indices.

  Returns:
      A new list containing modified tuples.
  """
  dx,dy,dz=value
  return [(element[0], element[1] + dx, element[2] + dy, element[3]) for element in data]
      
def generate_points(radius, num_points, min_neighbor_distance):
  """
  Generates a set of uniformly distributed points within a circle, ensuring
  a minimum distance between each point.

  Args:
      radius: The radius of the circle.
      num_points: The desired number of points.
      min_neighbor_distance: The minimum distance between points.

  Returns:
      A list of tuples representing the (x, y) coordinates of the points.
  """
  points = []
  while len(points) < num_points:
    x = random.uniform(-radius, radius)
    y = random.uniform(-radius, radius)
    is_valid = True
    # Check for existing neighbors within minimum distance
    for existing_x, existing_y in points:
      if (x - existing_x)**2 + (y - existing_y)**2 < min_neighbor_distance**2:
        is_valid = False
        break
    if is_valid:
      points.append((x, y))

  return np.array(points)

def transform_data(data):
  """
  This function transforms a list of NumPy arrays into a list of tuples with an added third column.

  Args:
      data: A list of NumPy arrays, where each array represents a row of data.

  Returns:
      A list of tuples, where each tuple has the same elements as the corresponding NumPy array
      with an additional third element of value 0.
  """
  return [(row[0], float(row[1]), float(row[2]), 0) for row in data]

def generate_atomic_coordinates(radius, num_points, min_neighbor_distance, element):
  """
  This function generates atomic coordinates within a circle and adds an element symbol.

  Args:
      center: A tuple representing the center coordinates (x, y) of the circle.
      radius: The radius of the circle.
      num_points: The number of random points to generate.
      element: The element symbol (e.g., 'H', 'C', 'O')

  Returns:
      A numpy array of shape (num_points, 4) containing atomic coordinates and element.
  """
  coordinates = generate_points(radius, num_points, min_neighbor_distance)

  # Fix for dimension mismatch: Repeat element symbol for each coordinate
  element_repeated = np.repeat(element, num_points)

  # Now the dimensions for concatenation match
  atomic_coordinates = np.column_stack((element_repeated, coordinates))

  return atomic_coordinates
  
# Example usage
radius = 5.0
num_points = 5
min_neighbor_distance = 1.2
element = "C"

atomic_coordinates=generate_atomic_coordinates(radius, num_points, min_neighbor_distance, element)
energy=run_b3lyp_calculation(atomic_coordinates)
write_xyz_file(transform_data(atomic_coordinates), energy, "initial.xyz")

#Monte Carlo Optimization
num_steps = 35  # Number of Monte Carlo steps
temperature = 1.25  # Temperature parameter (adjust for acceptance rate)

# Perform Monte Carlo optimization
metropolis_monte_carlo(atomic_coordinates, num_steps, temperature, run_b3lyp_calculation, 0.5)

