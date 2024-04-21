import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import Locator
import random

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
      energy += -J * spin_matrix[i, j] * (spin_matrix[di,j] + spin_matrix[i,dj])
  return energy

def BPSO(spin_matrix, J, max_iter, num_particles):
  # Initialize particles with random positions (0 or 1)
  #particles = np.random.randint(2, size=(num_particles, spin_matrix.shape[0], spin_matrix.shape[1]))
  
  num_zeros = int(num_particles * spin_matrix.shape[0] * spin_matrix.shape[1] / 2)
  particles = np.random.choice([-1, 1], size=(num_particles, spin_matrix.shape[0], spin_matrix.shape[1]), p=[0.5, 0.5])
  velocities = np.zeros_like(particles)
  pbest = particles.copy()  # Personal best positions
  gbest = particles[np.random.randint(num_particles)]  # Pick a random particle as initial gbest
  gbest_energy = get_energy(gbest, J)  # Calculate initial gbest energy

  for _ in range(max_iter):
    # Update velocities
    for i in range(num_particles):
      cognitive = np.random.rand() * (pbest[i] - particles[i])
      social = np.random.rand() * (gbest != particles[i])  # Convert to -1, 1 for difference
      velocities[i] = velocities[i] + cognitive + social

    # Apply velocity clamping (optional)
    velocities = np.clip(velocities, -1, 1)

    # Update positions (using sigmoid for probability)
    particles = particles + velocities
    particles = (1 / (1 + np.exp(-particles))) > 0.5  # Convert to binary based on probability

    # Calculate energy for each particle
    energies = [get_energy(particle, J) for particle in particles]

    # Update personal best
    for i, energy in enumerate(energies):
      if energy < get_energy(pbest[i], J):
        pbest[i] = particles[i].copy()

    # Update global best (check for strictly lower energy)
    best_idx = np.argmin(energies)
    if energies[best_idx] < gbest_energy:
      gbest = particles[best_idx].copy()
      gbest_energy = energies[best_idx]

  # Convert gbest from 0/1 to -1/1
  return gbest * 2 - 1, gbest_energy

# Example usage
# Define spin matrix (replace with your actual protein data)
spin_matrix=np.random.choice([H, P], size=(25,25), p=[0.5, 0.5])
# Parameters for BPSO
J = 1  # Interaction strength
max_iter = 100
num_particles = 20

print ("Initial matrix")
print (spin_matrix)
print ("Number spin down:", count_minus_ones(spin_matrix))
print ("Initial energy", get_energy(spin_matrix,J))

# Run BPSO
best_config, best_energy = BPSO(spin_matrix.copy(), J, max_iter, num_particles)
print("Best spin configuration:")
print(best_config)
# Print best energy
print("Best energy:", best_energy)
print ("Number spins down:", count_minus_ones(best_config))

# Optional: Visualization (modify to use -1 and 1 for colors)
plt.imshow(best_config, cmap='coolwarm')  # Red for -1, blue for 1
plt.colorbar(label='Spin')
plt.title('Best Spin Configuration (-1/1)')
plt.show()





