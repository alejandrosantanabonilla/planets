import numpy as np
import matplotlib.pyplot as plt

# Define spin values (H = Hydrophobic, P = Polar)
H, P = -1, 1

def count_minus_ones(arr):
  """Counts the number of -1s in a numpy array.

  Args:
    arr: A numpy array.

  Returns:
    The number of -1s in the array.
  """

  return (arr == -1).sum()

def get_energy(spin_matrix, J):
  """
  Calculates the total energy of the protein configuration.

  Args:
      spin_matrix: A 2D numpy array representing the protein configuration (H or P).
      J: The interaction strength between neighboring residues.

  Returns:
      The total energy of the configuration.
  """
  N = len(spin_matrix)
  energy = 0
  for i in range(N):
    for j in range(N):
      # Consider periodic boundary conditions
      di = (i + 1) % N
      dj = (j + 1) % N
      energy += -J * spin_matrix[i, j] * spin_matrix[di, dj]
  return energy

def monte_carlo_step(spin_matrix, J, beta):
  """
  Performs a single Monte Carlo step using Metropolis acceptance criterion.

  Args:
      spin_matrix: A 2D numpy array representing the protein configuration.
      J: The interaction strength between neighboring residues.
      beta: The inverse temperature (controls acceptance probability).

  Returns:
      A new spin matrix after the Monte Carlo step.
  """
  N = len(spin_matrix)
  # Randomly select a spin to flip
  i, j = np.random.randint(0, N, size=2)
  # Calculate energy difference after flipping
  delta_E = 2 * J * spin_matrix[i, j] * (spin_matrix[(i+1)%N, j] + spin_matrix[i, (j+1)%N])
  # Metropolis acceptance probability
  acceptance_prob = np.exp(-beta * delta_E)

  # Flip the spin based on acceptance probability
  if np.random.random() < acceptance_prob:
    spin_matrix[i, j] *= -1
    
  return spin_matrix

def simulate(spin_matrix, J, beta, n_steps):
  """
  Simulates protein folding using Monte Carlo simulations.

  Args:
      spin_matrix: A 2D numpy array representing the initial protein configuration.
      J: The interaction strength between neighboring residues.
      beta: The inverse temperature.
      n_steps: The number of Monte Carlo steps.

  Returns:
      The final spin matrix after the simulation.
  """
  for _ in range(n_steps):
    spin_matrix = monte_carlo_step(spin_matrix.copy(), J, beta)
  return spin_matrix

def plot_configuration(spin_matrix):
  """
  Plots the protein configuration as a colored grid.

  Args:
      spin_matrix: A 2D numpy array representing the protein configuration.
  """
  plt.matshow(spin_matrix, cmap='bwr')
  plt.colorbar(label='Spin (H: -1, P: 1)')
  plt.title('Protein Configuration')
  plt.show()

# Example: Simple 4x4 protein with random initial configuration
L = 25
# Simulation parameters
J = 1.0  # Interaction strength
beta = 0.5  # Inverse temperature
n_steps = 50000  # Number of Monte Carlo steps
spin_matrix=np.random.choice([H, P], size=(L,L), p=[0.5, 0.5])
print (spin_matrix)
print ("Number of spin downs:", count_minus_ones(spin_matrix))
print ("Initial energy: {}".format(get_energy(spin_matrix,J)))
plot_configuration(spin_matrix)


# Simulate protein folding
final_config = simulate(spin_matrix.copy(), J, beta, n_steps)
# Plot the final configuration
plot_configuration(final_config)
print (final_config)
print ("Spin down:", count_minus_ones(final_config))
print ("Final energy: {}".format(get_energy(final_config,J)))
