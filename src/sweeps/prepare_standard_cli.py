import sys

path_of_lib_pre = '/'.join(sys.argv[0].split('/')[:-2])
sys.path.append(path_of_lib_pre)
from lib_pre import HoneycombClusterPresets, TwoBodyHamiltonian, StandardInput, \
                    parameterize_multipole_by_angles, parameterize_special_points_by_angle
#--------------------parameters for type of calculation
model          = 'SpinGC'
method         = sys.argv[1]
lattice        = 'Honeycomb Lattice'
#--------------------lattice parameters
sites          = int(sys.argv[2])
shape          = sys.argv[3]
a0, a1         = HoneycombClusterPresets().clusters[f'{sites}-'+shape]
#--------------------conserved quantities
two_sz         = 0              #will be set to None if 'SpinGC' is selected
#--------------------parameters for the numerical condition
restart        = sys.argv[4]      #select 'None', 'Restart_out', 'Restart_in', 'Restart'
lanczos_max    = int(sys.argv[5]) #number of Lanczos/LOBCG steps
exct           = int(sys.argv[6]) #number of states to converge
lanczos_target = exct+1           #target eigenenergy for convergence
output_mode    = sys.argv[7]      #select 'none', 'correlation', 'full'
eigenvec_io    = 'None'           #select 'None', 'Out', 'In'
ham_io         = 'None'
#--------------------Hamiltonian parameters
ham_model = sys.argv[8] #'gen_jtau', 'xi', 'jkggp'

params = sys.argv[9:]
params_float = list(map(float, params))

if ham_model == 'gen_jtau':
    #model 1: general jtau, jb, jq, jo
    [theta, phi, jb] = params_float
    jt, jb, jq, jo = parameterize_multipole_by_angles(theta, phi, jb)
    H0, H1, H2     = TwoBodyHamiltonian().make_multipole_hamiltonian(jt, jb, jq, jo)
elif ham_model == 'xi':
    #model 2: jtau and jq fixed, tuning jo and jb
    [xi]             = params_float
    jt, jb, jq, jo = parameterize_special_points_by_angle(xi)
    H0, H1, H2     = TwoBodyHamiltonian().make_multipole_hamiltonian(jt, jb, jq, jo)
elif ham_model == 'jkggp':
    #model 3: j-k-g-gp model
    [j, k, g, gp]    = params_float
    H0, H1, H2     = TwoBodyHamiltonian().make_kitaev_hamiltonian(j, k, g, gp)
#--------------------create instance create file, and output.
simulation     = StandardInput(model, method, lattice,
                               a0, a1,
                               two_sz,
                               H0, H1, H2,
                               restart, lanczos_max, exct, lanczos_target, output_mode,
                               eigenvec_io, ham_io)
simulation.check_availability() #since I only implemented a subset of HPhi features
simulation.ensure_consistency()

simulation.create_file('stan.in')
