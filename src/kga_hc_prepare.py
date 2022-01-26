from lib_prepare import HoneycombPresets, TwoBodyHamiltonian, StandardInput, pi
from math import sin, cos, sqrt
import os
#--------------------parameters for type of calculation
model          = 'SpinGC'
method         = 'CG'
lattice        = 'Honeycomb Lattice'
#--------------------lattice parameters
a0, a1         = HoneycombPresets().clusters['24-site']
#--------------------conserved quantities
two_sz         = 0              #will be set to None if 'SpinGC' is selected
#--------------------Hamiltonian parameters
j, k, g, gp    = -1, 0, 0, 0
H0, H1, H2     = TwoBodyHamiltonian().make_kitaev_hamiltonian(j, k, g, gp)
#--------------------parameters for the numerical condition
restart        = 'None'
lanczos_max    = 2000           #number of Lanczos/LOBCG steps
exct           = 2              #number of states to converge
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
filename = 'stan.in'
simulation.create_file(output_path+filename)
