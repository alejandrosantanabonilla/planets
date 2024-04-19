import random
import numpy as np
import matplotlib.pyplot as plt

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

def plot_points(points, radius):
  """
  Plots the generated points within a circle.

  Args:
      points: A list of tuples representing the (x, y) coordinates of the points.
      radius: The radius of the circle.
  """
  x, y = zip(*points)
  plt.scatter(x, y)
  circle = plt.Circle((0, 0), radius, color='blue', fill=False)
  plt.gca().add_patch(circle)
  plt.xlabel("X")
  plt.ylabel("Y")
  plt.xlim(-radius - radius/5, radius + radius/5)
  plt.ylim(-radius - radius/5, radius + radius/5)
  plt.title(f"{len(points)} points with minimum distance of {min_neighbor_distance}")
  plt.show()
  
# Example usage
radius = 10
num_points = 20
min_neighbor_distance = 2
generated_points = generate_points(radius, num_points, min_neighbor_distance)
plot_points(generated_points, radius)
