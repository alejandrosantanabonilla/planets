import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def sample_poisson_disk_sphere(center, radius, min_dist, num_points):
    cell_size = min_dist / np.sqrt(3)  # Cell size for 3D grid
    grid_size = np.ceil(radius * 2 / cell_size).astype(int)  # Grid dimensions
    grid = {}  # Dictionary to store points in the grid
    
    # Helper function to generate a random candidate point within the sphere
    def generate_candidate():
        angle1 = np.random.uniform(0, 2 * np.pi)
        angle2 = np.arccos(np.random.uniform(-1, 1))
        r = np.random.uniform(0, radius)
        x = center[0] + r * np.sin(angle2) * np.cos(angle1)
        y = center[1] + r * np.sin(angle2) * np.sin(angle1)
        z = center[2] + r * np.cos(angle2)
        return np.array([x, y, z])

    # Helper function to check if a point is valid (not too close to existing points)
    def is_valid(p):
        grid_x = int(np.floor((p[0] - center[0] + radius) / cell_size))
        grid_y = int(np.floor((p[1] - center[1] + radius) / cell_size))
        grid_z = int(np.floor((p[2] - center[2] + radius) / cell_size))

        for i in range(max(0, grid_x - 2), min(grid_x + 3, grid_size)):
            for j in range(max(0, grid_y - 2), min(grid_y + 3, grid_size)):
                for k in range(max(0, grid_z - 2), min(grid_z + 3, grid_size)):
                    if (i, j, k) in grid:
                        if np.linalg.norm(p - grid[(i, j, k)]) < min_dist:
                            return False
        return True

    active_list = []  # List of active points
    initial_point = generate_candidate()
    active_list.append(initial_point)
    grid[tuple(np.floor((initial_point - center + radius) / cell_size).astype(int))] = initial_point

    points = [initial_point]  # List to store the generated points

    while active_list and len(points) < num_points:
        index = np.random.randint(len(active_list))
        point = active_list[index]
        
        for _ in range(30):  # Number of attempts to find a valid neighbor
            new_point = point + np.random.uniform(-1, 1, 3) * min_dist
            if np.linalg.norm(new_point - center) <= radius and is_valid(new_point):
                active_list.append(new_point)
                grid[tuple(np.floor((new_point - center + radius) / cell_size).astype(int))] = new_point
                points.append(new_point)
                if len(points) >= num_points:
                    break
        else:
            del active_list[index] 
    
    return points

# Example usage
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

center = np.array([0, 0, 0])  # Sphere center
radius = 2.5
min_dist = 1.42  # Minimum distance between points
num_points = 20  # Desired number of points

points = sample_poisson_disk_sphere(center, radius, min_dist, num_points)

# Convert the points from a list to a numpy array for easier slicing
points = np.array(points)

print(f'Generated {len(points)} points')

# Plot the points
ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='b', s=10)

# Plot the sphere surface (optional)
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = radius * np.outer(np.cos(u), np.sin(v))
y = radius * np.outer(np.sin(u), np.sin(v))
z = radius * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='gray', alpha=0.2)

# Set axis labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title(f'Poisson Disk Sampling on a Sphere (r={radius}, min_dist={min_dist}, num_points={num_points})')

# Show the plot
plt.show()
