import numpy as np
from pyscf import gto, scf, dft

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

def generate_random_points_in_circle(center, radius, num_points, min_distance):
  """
  This function generates a specified number of random points within a circle,
  ensuring a minimum distance between points.

  Args:
      center: A tuple representing the center coordinates (x, y) of the circle.
      radius: The radius of the circle.
      num_points: The number of random points to generate.
      min_distance: The minimum allowed distance between generated points.

  Returns:
      A numpy array of shape (num_points, 2) containing the random points.
  """

  # Initialize an empty list to store valid points
  valid_points = []

  # Loop until the desired number of valid points is generated
  while len(valid_points) < num_points:
    # Generate random angles and distances as before
    theta = np.random.rand(num_points) * 2 * np.pi
    r = np.random.rand(num_points) * radius

    # Convert polar coordinates to cartesian coordinates
    x = r * np.cos(theta) + center[0]
    y = r * np.sin(theta) + center[1]
    candidate_points = np.array([x, y]).T  # Combine coordinates

    # Check if any candidate point is too close to existing valid points
    valid = True
    for valid_point in valid_points:
      # Calculate distances between candidate points and existing points
      distances = np.linalg.norm(candidate_points - valid_point, axis=1)
      if np.min(distances) < min_distance:
        valid = False
        break  # Exit inner loop if a close point is found

    # If all candidate points are valid, add them to the list
    if valid:
      valid_points.extend(candidate_points.tolist())

  # Return the generated points after ensuring minimum distance
  return np.array(valid_points)

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

def generate_atomic_coordinates(center, radius, num_points, element, min_distance):
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
  coordinates = generate_random_points_in_circle(center, radius, num_points, min_distance)

  # Fix for dimension mismatch: Repeat element symbol for each coordinate
  element_repeated = np.repeat(element, num_points)

  # Now the dimensions for concatenation match
  atomic_coordinates = np.column_stack((element_repeated, coordinates))

  return atomic_coordinates

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
  best_energy = get_energy(best_coordinates)  # Initial energy

  for _ in range(num_steps):
    displacements = np.random.uniform(low=-dx, high=dx, size=3)  # Adjust shape
    new_atom=list(transform_data(atomic_coordinates))
    trial_coordinates=add_to_indices(new_atom, displacements)

    # Calculate the energy of the trial configuration
    trial_energy = get_energy(trial_coordinates)

    # Metropolis acceptance criterion
    delta_energy = trial_energy - best_energy
    acceptance_prob = np.exp(-delta_energy / temperature)

    # Accept or reject the trial configuration
    if np.random.rand() < acceptance_prob:
      best_coordinates = trial_coordinates.copy()
      best_energy = trial_energy

  return best_coordinates, best_energy

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

  
# Example usage
center = (0, 0)  # Center coordinates
radius = 2.5        # Circle radius
num_points = 3    # Number of points to generate
element = 'C'     # Element symbol
min_distance = 1.0

atomic_coordinates = generate_atomic_coordinates(center, radius, num_points, element, min_distance)
energy=run_b3lyp_calculation(atomic_coordinates)

#Monte Carlo Optimization
num_steps = 10  # Number of Monte Carlo steps
temperature = 0.75  # Temperature parameter (adjust for acceptance rate)

# Perform Monte Carlo optimization
optimized_coordinates, optimized_energy = metropolis_monte_carlo(atomic_coordinates, num_steps, temperature, run_b3lyp_calculation, 0.15)
print("Optimized Energy:", optimized_energy)
write_xyz_file(optimized_coordinates, optimized_energy, "result.xyz")
