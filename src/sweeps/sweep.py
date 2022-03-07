from lib_sweep import HPhiSweeps
import os
import numpy as np
from math import cos, sin
#---------------------------------------
#---------------------logistical details
#---------------------------------------
run = 1
what_computer = 'niagara'  #can be either 'laptop', 'home', or 'niagara'

Ncorespernode = 4
hyperthreadQ = False
Nmpi          = pow(2,1)               #MUST BE a power of 2 for S=1/2
if (Ncorespernode % Nmpi == 0):
    Nomp      = int(Ncorespernode/Nmpi/2) #SHOULD BE AN INTEGER
else:
    print("Check the cores! N_corespernode/N_mpi/N_omp should be an integer.")
    raise SystemExit

Nnodes        = 1
time          = '00:15:00'
hpc_settings = [Nnodes, hyperthreadQ, Nomp, Nmpi, Ncorespernode, time]
#-----------------------------------------------------------
#-----------------------------general command line arguments
#-----------------------------------------------------------
method         = 'CG'
sites          = 6
shape          = 'RH120'
restart        = 'None'    #select 'None', 'Restart_out', 'Restart_in', 'Restart'
lanczos_max    = 2000      #number of Lanczos/LOBCG steps
exct           = 1         #number of states to converge
output_mode    = 'none'    #select 'none', 'correlation', 'full'
ham_model      = 'eps'          #gen_jtau, eps, jkggp

stan_cli_list = [method, sites, shape, restart, lanczos_max,
                   exct, output_mode, ham_model]
#-------------------------------------------------------
#---------------------parameter entry: min, max, spacing
#-------------------------------------------------------
if ham_model == 'gen_jtau':
    # model 1: general jtau, jb, jq, jo
    pass
elif ham_model == 'eps':
    #model 2: jtau and jq fixed, tuning jo and jb
    eps_val_list, eps_label = [-1.000, 0.000, 1.000], "eps"
    params_list = [eps_val_list]
    params_label_list = [eps_label]
elif ham_model == 'jkggp':
    #model 3: j-k-g-gp model
    pass

#----magnetic field magnitude
h_val_list,  h_label = [0.000, 0.000, 0.20], "h"

params_list += [h_val_list]
params_label_list += [h_label]

#----magnetic field direction
#angles in degrees
htheta = 0
hphi   = 0

X, Y, Z = np.eye(3)
h0 = sin(htheta/180*np.pi)*cos(hphi/180*np.pi) #"x"
h1 = sin(htheta/180*np.pi)*sin(hphi/180*np.pi) #"y"
h2 = cos(htheta/180*np.pi)                     #"z"

hdirection = h1*X + h2*Y + h0*Z #theta = 0 is +Y, theta=90 and phi=0 is +Z
hdir_norm = hdirection/np.linalg.norm(hdirection)
#-------------------------------------------------------------
#---------------------create instance create file, and output.
#-------------------------------------------------------------
cwd = os.getcwd()
sweep = HPhiSweeps(run, what_computer,
                   hpc_settings,
                   stan_cli_list,
                   params_list, params_label_list, hdir_norm,
                   cwd)
