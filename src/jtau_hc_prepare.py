from lib_prepare import HoneycombClusterPresets, TwoBodyHamiltonian, StandardInput, pi
from math import sin, cos, atan
import os

def parameterize_multipole_by_angles(theta, phi, jb):
    return [cos(theta*pi), jb, sin(theta*pi)*cos(phi*pi), sin(theta*pi)*sin(phi*pi)]

#--------------------parameters for type of calculation
model          = 'SpinGC'
method         = 'Lanczos'
lattice        = 'Honeycomb Lattice'
#--------------------lattice parameters
sites          = 4
shape          = 'RE'
a0, a1         = HoneycombClusterPresets().clusters[f'{sites}-'+shape]
#--------------------conserved quantities
two_sz         = 0              #will be set to None if 'SpinGC' is selected
#--------------------Hamiltonian parameters
theta, phi, jb = atan(1/2)/pi, 0.5, 0
jt, jb, jq, jo = parameterize_multipole_by_angles(theta, phi, jb)
H0, H1, H2     = TwoBodyHamiltonian().make_multipole_hamiltonian(jt, jb, jq, jo)
#--------------------parameters for the numerical condition
restart        = 'None'
lanczos_max    = 2000           #number of Lanczos/LOBCG steps
exct           = 3              #number of states to converge
lanczos_target = exct+1         #target eigenenergy for convergence
output_mode    = 'correlation'
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
filename = 'stan.in'
simulation.create_file(output_path+filename)
