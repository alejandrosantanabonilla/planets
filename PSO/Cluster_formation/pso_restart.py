import random
import numpy as np
from pyscf import gto, scf, dft
import numpy as np
import json

def save_structure_json(data, energy, iteration, filename="structure_data.json"):
    """Saves structure coordinates and energy to a JSON file."""
    structure = {
        "iteration": iteration,
        "coordinates": data.tolist(),  # Convert to list for JSON serialization
        "energy": energy
    }
    with open(filename, "a") as f:  # Append to file
        json.dump(structure, f)
        f.write("\n")  # Add newline for readability

def load_last_structure(filename="structure_data.json"):
    """Loads the last structure from the JSON file."""
    try:
        with open(filename, "r") as f:
            lines = f.readlines()
        if lines:
            last_structure = json.loads(lines[-1])
            return np.array(last_structure["coordinates"]), last_structure["iteration"], last_structure["energy"]
    except FileNotFoundError:
        return None, 0, None

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

def run_pbe_calculation(atomic_coordinates):
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
   
  mf = dft.ROKS(mol)
  mf.xc='PBE'
  mf.nlc='VV10'
  
  return mf.kernel()

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

def pso_optimize_structure(atomic_coordinates, max_iter=10, swarm_size=15, c1=2.0, c2=2.0,
                           w_min=0.1, w_max=0.9, dx=0.25, restart=False):
    """Optimizes the energy of a structure using the Particle Swarm Optimization (PSO) algorithm.

    Args:
        atomic_coordinates: A numpy array of shape (num_atoms, 4) containing atomic coordinates and element (e.g., [['C', ...], ...]).
        max_iter: Maximum number of iterations (default: 100).
        swarm_size: Number of particles in the swarm (default: 15).
        c1: Cognitive learning rate (default: 2.0).
        c2: Social learning rate (default: 2.0).
        w_min: Minimum inertia weight (default: 0.1).
        w_max: Maximum inertia weight (default: 0.9).
        dx: Maximum displacement for each coordinate (default: 0.25).
        restart: Whether to restart from the last saved structure (default: False).

    Returns:
        The optimized atomic coordinates and the corresponding energy.
    """

    num_atoms, num_coords = atomic_coordinates.shape

    # Load last structure if restarting
    start_iter = 0
    if restart:
        last_coords, start_iter, last_energy = load_last_structure()
        if last_coords is not None:
            try:  # Try to prioritize loading coordinates from file
                if not np.array_equal(last_coords, atomic_coordinates):
                    raise ValueError("New coordinates provided with restart=True. Using coordinates from structure_data.json.")
            except ValueError as e:
                print(e) 
            atomic_coordinates = last_coords

    # Set up particle swarm with displacements as array of arrays repeated for swarm_size
    displacements = generate_matrix_list(swarm_size, num_atoms, num_coords, dx)
    swarm = np.zeros_like(displacements)
    velocity = np.zeros_like(swarm)
    pbest = np.zeros_like(displacements)
    pbest_energy = np.zeros(swarm_size)

    for i in range(swarm_size):
        displaced_coordinates = sum_displacements(displacements[i][:, -2:], atomic_coordinates.copy())
        pbest_energy[i] = run_pbe_calculation(displaced_coordinates)

    gbest = pbest[pbest_energy.argmin()]
    gbest_energy = min(pbest_energy)

    print("Main PSO loop")
    best_coordinates = atomic_coordinates.copy()  # Initialize best_coordinates
    iter = 0  # Initialize iter before the loop
    for iter in range(start_iter, max_iter):  # Start from loaded iteration if restarting
        print("Iteration number:", iter)

        # Update inertia weight
        w = w_max - (iter / (max_iter - 1)) * (w_max - w_min)
        r1 = np.random.uniform(0, 1, 1)[0]
        r2 = np.random.uniform(0, 1, 1)[0]

        # Update velocity
        velocity = w * velocity + c1 * r1 * (pbest - swarm) + c2 * r2 * (gbest - swarm)

        # Update displacements
        displacements = displacements + velocity
        swarm = displacements

        # Evaluate new positions using displacements
        for i in range(swarm_size):
            displaced_coordinates = sum_displacements(displacements[i][:, -2:], atomic_coordinates.copy())
            new_energy = run_pbe_calculation(displaced_coordinates)

            if new_energy < pbest_energy[i]:
                pbest[i] = displacements[i]
                pbest_energy[i] = new_energy

        if min(pbest_energy) < gbest_energy:
            gbest = pbest[np.argmin(pbest_energy)]
            gbest_energy = min(pbest_energy)
            print("Gbest energy:", gbest_energy)
            best_coordinates = sum_displacements(gbest[:, -2:], atomic_coordinates.copy())

        # Save structure after each iteration
        save_structure_json(best_coordinates, gbest_energy, iter + 1)  # Save with iter + 1 

    print("Successfully terminated after:", iter)
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
radius = 4.5
num_points = 10
min_neighbor_distance = 1.25
element = "C"

atomic_coordinates=generate_atomic_coordinates(radius, num_points, min_neighbor_distance, element)
#energy=run_pbe_calculation(atomic_coordinates)
#write_xyz_file(transform_data(atomic_coordinates), energy, "initial.xyz")

#PSO Optimization Example
max_iter=10
swarm_size=15
c1=1.42
c2=1.42
w_min=0.35
w_max=0.95

#Maximum displacement for the random number
dx=0.95

restart = True  # Set to True if you want to restart
best_coord, best_energy = pso_optimize_structure(atomic_coordinates, max_iter, swarm_size, c1, c2, w_min, w_max, dx, restart=restart)
write_xyz_file(transform_data(best_coord), best_energy, "final.xyz")

