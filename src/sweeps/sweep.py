from lib_sweep import HPhiSweeps
import os
import numpy as np
from math import cos, sin
#---------------------------------------
#---------------------logistical details
#---------------------------------------
run = 1
what_computer = 'laptop'  #can be either 'laptop', 'home', or 'niagara'

#the following are only really important for mode 'niagara'
num_processes = 4         #isn't used unless 'niagara' is selected
ntasks_per_node = 80
nodes = 1
time = '01:45:00'
hpc_settings = [num_processes, ntasks_per_node, nodes, time]
#-----------------------------------------------------------
#-----------------------------general command line arguments
#-----------------------------------------------------------
method         = 'CG'
sites          = 6
shape          = 'RH120'
restart        = 'None'    #select 'None', 'Restart_out', 'Restart_in', 'Restart'
lanczos_max    = 2000      #number of Lanczos/LOBCG steps
exct           = 4         #number of states to converge
output_mode    = 'correlation'    #select 'none', 'correlation', 'full'
ham_model      = 'xi'          #gen_jtau, xi, jkggp

stan_cli_list = [method, sites, shape, restart, lanczos_max,
                   exct, output_mode, ham_model]
#-------------------------------------------------------
#---------------------parameter entry: min, max, spacing
#-------------------------------------------------------
if ham_model == 'gen_jtau':
    # model 1: general jtau, jb, jq, jo
    pass
elif ham_model == 'xi':
    #model 2: jtau and jq fixed, tuning jo and jb
    xi_val_list, xi_label = [0.300, 0.360, 0.003], "xi"
    params_list = [xi_val_list]
    params_label_list = [xi_label]
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
