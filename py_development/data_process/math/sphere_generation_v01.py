#!/usr/bin/env python3
#spherical surface point generation

import math
import numpy as np

def generate_sphere_points(n_points): 
  inc = math.pi * (3.0 - np.sqrt(5.0))
  offset = 2.0 / float(n_points)
  sphere_points = np.empty((n_points,3),dtype=np.float32)

  for i in range(n_points):
    y = i*offset - 1.0 + (offset / 2.0)
    r = np.sqrt(1.0 - y*y)
    phi = i * inc
    sphere_points[i] = np.array([math.cos(phi) * r, y, math.sin(phi) * r], dtype=np.float32)

  return sphere_points

sphere_points=generate_sphere_points(100)

for row in sphere_points:
  print('X {:8.3f} {:8.3f} {:8.3f}'.format(row[0],row[1],row[2]))

