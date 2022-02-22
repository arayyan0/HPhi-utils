from lib_pre import HoneycombClusterPresets, TwoBodyHamiltonian, StandardInput, \
                    parameterize_multipole_by_angles, \
                    parameterize_special_points_by_angle, pi
from math import sin, cos, atan, sqrt
import numpy as np
import os

#--------------------parameters for type of calculation
model          = 'SpinGC'
method         = 'CG'
lattice        = 'Honeycomb Lattice'
#--------------------lattice parameters
sites          = 6
shape          = 'RH120'
a0, a1         = HoneycombClusterPresets().clusters[f'{sites}-'+shape]
#--------------------conserved quantities
two_sz         = 0              #will be set to None if 'SpinGC' is selected
#--------------------Hamiltonian parameters
#model 1: general jtau, jb, jq, jo
# theta, phi, jb = atan(1/sqrt(2))/pi, 0.75, 0
# jt, jb, jq, jo = parameterize_multipole_by_angles(theta, phi, jb)
# H0, H1, H2     = TwoBodyHamiltonian().make_multipole_hamiltonian(jt, jb, jq, jo)
#model 2: jtau and jq fixed, tuning jo and jb
xi             = 0
jt, jb, jq, jo = parameterize_special_points_by_angle(xi)
H0, H1, H2     = TwoBodyHamiltonian().make_multipole_hamiltonian(jt, jb, jq, jo)
#model 3: j-k-g-gp model
# j, k, g, gp    = 0, 1, 1, 0
# H0, H1, H2     = TwoBodyHamiltonian().make_kitaev_hamiltonian(j, k, g, gp)
#--------------------parameters for the numerical condition
restart        = 'None'
lanczos_max    = 2000           #number of Lanczos/LOBCG steps
exct           = 3              #number of states to converge
lanczos_target = exct+1         #target eigenenergy for convergence
output_mode    = 'none'
eigenvec_io    = 'None'
ham_io         = 'None'         #will be set to None unless 'FullDiag' is selected
#--------------------create instance create file, and output.
simulation     = StandardInput(model, method, lattice,
                               a0, a1,
                               two_sz,
                               H0, H1, H2,
                               restart, lanczos_max, exct, lanczos_target, output_mode,
                               eigenvec_io, ham_io)
simulation.check_availability() #since I only implemented a subset of HPhi features
simulation.ensure_consistency()

output_path = f'out/'
if not os.path.exists(output_path):
    os.makedirs(output_path)
simulation.create_file(output_path+'stan.in')
