import matplotlib.pyplot as plt
import numpy as np

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

def plot_random_points_in_circle(center, radius, num_points, min_distance):
  """
  This function generates random points and plots them within a circle using matplotlib.

  Args:
      center: A tuple representing the center coordinates (x, y) of the circle.
      radius: The radius of the circle.
      num_points: The number of random points to generate.
  """
  points = generate_random_points_in_circle(center, radius, num_points, min_distance)

  # Create the circle patch
  circle = plt.Circle(xy=center, radius=radius, color='blue', alpha=0.2)

  # Plot the points and the circle
  plt.scatter(points[:, 0], points[:, 1], color='red')
  plt.gca().add_patch(circle)

  # Set aspect ratio
  plt.gca().set_aspect('equal')

  # Label axes
  plt.xlabel("X")
  plt.ylabel("Y")

  # Set limits slightly bigger than radius
  plt.xlim(-radius - radius/5 + center[0], radius + radius/5 + center[0])
  plt.ylim(-radius - radius/5 + center[1], radius + radius/5 + center[1])

  plt.title(f"{num_points} Random Points in Circle")
  plt.show()

# Example usage
center = (0, 0)  # Center coordinates
radius = 20.0        # Circle radius
num_points = 10  # Number of points to generate
min_distance=10.0

plot_random_points_in_circle(center, radius, num_points, min_distance)
