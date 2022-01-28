import numpy as np
from dataclasses import dataclass
from math import sqrt
import os

pi = np.pi

# @dataclass
# class ConventionalUnitCell:
#     a = np.eye(3)      #primitive unit cell is given by rows of a
#     b = np.eye(3)*2*pi #primitive reciprocal unit cell is given by rows of b
#     # a and b satisfy a @ b.T = 2pi np.eye(3)
#
#     R: np.ndarray      #conventional unit cell is given by rows of A = R @ a
#     def find_reciprocal_vectors(self):
#         #returns B such that A @ B.T = 2pi np.eye(3)
#         B = (self.b @ np.linalg.inv(self.R)).T
#         return B

class TriangularSymmetryPoints:
    symmetry_points = { }
    def add_symmetrypoint(self, label, vector):
        self.symmetry_points[label] = vector

    def __init__(self, path):
        self.add_symmetrypoint( 'X', np.array([-1/2, 1/2, 0]))
        self.add_symmetrypoint( 'K', np.array([-1/4, 1/4, 0]))
        self.add_symmetrypoint( 'G', np.array([0, 0, 0]))
        self.add_symmetrypoint('M2', np.array([ 1/2, 1/2, 0]))
        self.add_symmetrypoint('Kp', np.array([2/3, 1/3, 0]))
        self.add_symmetrypoint('Gp', np.array([1, 0, 0]))
        self.add_symmetrypoint('M1', np.array([1/2, 0, 0]))
        self.path = path

    def create_path(self):
        texts = [point + ' ' +
                  ' '.join(self.symmetry_points[point].round(3).astype(str))+'\n'
                 for point in self.path]
        text = ''.join(texts)
        return text #removes last newline character

    def append_to_geometry(self, filename, line_density, surface_density):
        existsQ = os.path.isfile(filename)
        if existsQ == False:
            print('geometry.dat does not exist!')
            raise SystemExit

        first_line = f'{len(path)} {line_density}\n'
        body = self.create_path()
        last_line  = f'{surface_density} {surface_density} 1\n'
        text = first_line + body + last_line

        f = open(filename, 'a')
        f.write(text)
        f.close()

# R  = np.array([[1, -1, 0], [1, 0, 0], [0, 0, 1]])
# A = R @ np.eye(3)
# B = ConventionalUnitCell(R).find_reciprocal_vectors()

path = ['X', 'K', 'G', 'M2', 'Kp', 'Gp', 'M1', 'G']
sym_points = TriangularSymmetryPoints(path)

geometry_file   = '../preHPHI/out/geometry.dat'
line_density    = 20 #number of k along the high symmetry points
surface_density = 24 #linear dimension of k-grid for ssf surface
sym_points.append_to_geometry(geometry_file, line_density, surface_density)
