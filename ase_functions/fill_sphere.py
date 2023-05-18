import numpy as np 
  
phi = (1 + np.sqrt(5)) / 2  # the golden ratio
long_incr = 2 * np.pi / phi  # how much to increment the longitude
dz = 2.0 / float(n)  # a unit sphere has diameter 2
bands = np.arange(int(n))  # each band will have one point placed on it
z = bands * dz - 1.0 + (dz / 2.0)  # the height z of each band/point
r = np.sqrt(1.0 - z * z)  # project onto xy-plane
az = bands * long_incr  # azimuthal angle of point modulo 2 pi
x = r * np.cos(az)
y = r * np.sin(az)
points = np.column_stack((x, y, z))
