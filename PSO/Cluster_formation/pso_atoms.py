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
  mol.basis = {'C': "def2-svp"}
  mol.build()
  
  # Perform B3LYP calculation (uncomment to run the calculation)
  mf_hf = dft.ROKS(mol)
  mf_hf.xc="b3lyp"
  mf_hf = mf_hf.newton()

  return mf_hf.kernel()

def generate_matrix_list(swarm_size, num_atoms, num_coords, dx=0.5):
  """
  Generates a list of  matrices with entries between -dx and dx.

  Args:
      dx: The upper and lower bound for the random matrix entries.

  Returns:
      A list of 20 randomly generated 5x3 NumPy matrices.
  """
  matrices = []
  for _ in range(swarm_size):
    matrices.append(np.random.uniform(low=-dx, high=dx, size=(num_atoms, num_coords)))
  return matrices

def pso_optimize_structure(atomic_coordinates, max_iter=100, swarm_size=15, c1=2.0, c2=2.0,
                           w_min=0.1, w_max=0.9, dx=0.25):
  """
  Optimizes the energy of a structure using the Particle Swarm Optimization (PSO) algorithm.

  Args:
      atomic_coordinates: A numpy array of shape (num_atoms, 4) containing
                          atomic coordinates and element (e.g., [['C', ...], ...]).
      run_b3lyp_calculations: A function that takes the atomic coordinates array
                              and returns the calculated energy.
      max_iter: Maximum number of iterations (default: 100).
      swarm_size: Number of particles in the swarm (default: 20).
      c1: Cognitive learning rate (default: 1.49618).
      c2: Social learning rate (default: 1.49618).
      w_min: Minimum inertia weight (default: 0.4).
      w_max: Maximum inertia weight (default: 0.9).
      dx: Maximum displacement for each coordinate (default: 0.1).

  Returns:
      The optimized atomic coordinates and the corresponding energy.
  """

  num_atoms, num_coords = atomic_coordinates.shape

  # Set up particle swarm with displacements as array of arrays repeated for swarm_size
  displacements = generate_matrix_list(swarm_size, num_atoms, num_coords)
  swarm = np.zeros_like(displacements)  
  velocity = np.zeros_like(swarm)
  pbest = np.zeros_like(displacements)
  pbest_energy = np.zeros(swarm_size)
 
  for i in range(swarm_size):
    displaced_coordinates = sum_displacements(displacements[i][:, -2:], atomic_coordinates.copy())
    pbest_energy[i] = run_b3lyp_calculation(displaced_coordinates)

  gbest = pbest[pbest_energy.argmin()]  # Global best
  gbest_energy = min(pbest_energy)

  print ("Main PSO loop")
  for iter in range(max_iter):
    # Update inertia weight
    print ("Iteration number:", iter)
    
    w = w_max - (iter / (max_iter - 1)) * (w_max - w_min)
    r1=np.random.uniform(0,1,1)[0]
    r2=np.random.uniform(0,1,1)[0]

    # Update velocity
    velocity = w * velocity + c1 * r1 * (pbest - swarm) + c2 * r1 * (gbest - swarm)

    # Update displacements (swarm remains empty)
    
    displacements = displacements + velocity
    
    print ("new swarm")
    # Evaluate new positions using displacements
    for i in range(swarm_size):
      displaced_coordinates = sum_displacements(displacements[i][:, -2:], atomic_coordinates.copy())
      new_energy = run_b3lyp_calculation(displaced_coordinates)
      
      if new_energy < pbest_energy[i]:
        pbest[i] = displacements[i]
        pbest_energy[i] = new_energy
      
    if min(pbest_energy) < gbest_energy:
       gbest = pbest[np.argmin(pbest_energy)]
       gbest_energy = min(pbest_energy)
       print (gbest_energy)
      
  best_coordinates = sum_displacements(gbest[:, -2:], atomic_coordinates.copy())

  print ("Succesfully terminated after: ", iter)
  return best_coordinates, gbest_energy
   
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
radius = 4.0
num_points = 3
min_neighbor_distance = 0.75
element = "C"
atomic_coordinates=generate_atomic_coordinates(radius, num_points, min_neighbor_distance, element)
energy=run_b3lyp_calculation(atomic_coordinates)
write_xyz_file(transform_data(atomic_coordinates), energy, "initial.xyz")

#PSO Optimization Example
max_iter=750
swarm_size=20
c1=1.42
c2=1.42
w_min=0.6
w_max=0.9
#Maximum displacement for the random number
dx=0.25

best_coord, best_energy=pso_optimize_structure(atomic_coordinates, max_iter, swarm_size, c1, c2, w_min, w_max, dx)
write_xyz_file(transform_data(best_coord), best_energy, "final.xyz")
