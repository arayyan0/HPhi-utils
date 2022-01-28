from lib_post import TriangularSymmetryPoints
import numpy as np
from math import sqrt

geometry_file   = '../out/geometry.dat'
with open(geometry_file,'r') as f:
    file_info = f.readlines()

A1, A2, A3  = [np.array(
                        str.replace('\n','').split(' ')
                        ).astype(np.int) for str in file_info[4:7]]
A = np.array([A1, A2, A3])

path = ['X', 'K', 'G', 'M2', 'Kp', 'Gp', 'M1', 'G']
sym_points = TriangularSymmetryPoints(A, path)

line_density    = 20 #number of k along the high symmetry points
surface_density = 24 #linear dimension of k-grid for ssf surface
sym_points.append_to_geometry(geometry_file, line_density, surface_density)
