from lib_sweep import HPhiSweeps, ComputerPresets, SLURMHelper
import os
import numpy as np
from math import cos, sin
#---------------------------------------
#---------------------logistical details
#---------------------------------------
run = 1

what_computer = 'laptop'  #can be either 'laptop', 'home', or 'niagara'
computer_settings = ComputerPresets().computers[what_computer]

Nnodes        = 1
hyperthreadQ  = False
Nmpi          = pow(2,0)               #MUST BE a power of 2 for S=1/2
Nomp          = 4
time          = '00:15:00'
slurm_helper  = SLURMHelper(computer_settings, Nnodes, hyperthreadQ, Nmpi, Nomp, time)

slurm_helper.calculate_relevant_integers()
slurm_helper.create_local_sim_commands()
slurm_helper.create_submit_script_texts()
#-----------------------------------------------------------
#-----------------------------general command line arguments
#-----------------------------------------------------------
method         = 'CG'
sites          = 24
shape          = 'RH120'
restart        = 'None'    #select 'None', 'Restart_out', 'Restart_in', 'Restart'
lanczos_max    = 2000      #number of Lanczos/LOBCG steps
exct           = 1         #number of states to converge
output_mode    = 'none'    #select 'none', 'correlation', 'full'
ham_model      = 'jtaujbjqjo'          #gen_jtau, eps, jkggp

stan_cli_list = [method, sites, shape, restart, lanczos_max,
                   exct, output_mode, ham_model]
#-------------------------------------------------------
#---------------------parameter entry: min, max, spacing
#-------------------------------------------------------
if ham_model == 'jtaujbjqjo':
    # model 1: general jtau, jb, jq, jo
    jtau_val_list, jtau_label = [-1.000, 1.000, 1.000], "jtau"
    jb_val_list, jb_label     = [ 1.000, 1.000, 1.000], "jb"
    jq_val_list, jq_label     = [ 0.000, 0.000, 1.000], "jq"
    jo_val_list, jo_label     = [ 0.000, 0.000, 1.000], "jo"
    params_list = [jtau_val_list, jb_val_list, jq_val_list, jo_val_list]
    params_label_list = [jtau_label, jb_label, jq_label, jo_label]

elif ham_model == 'eps':
    #model 2: jtau and jq fixed, tuning jo and jb
    eps_val_list, eps_label = [0.000, 0.000, 0.000], "eps"
    params_list = [eps_val_list]
    params_label_list = [eps_label]

elif ham_model == 'jkggp':
    #model 3: j-k-g-gp model
    j_val_list, j_label   = [-1.000, 1.000, 1.000], "j"
    k_val_list, k_label   = [ 1.000, 1.000, 1.000], "k"
    g_val_list, g_label   = [ 0.000, 0.000, 1.000], "g"
    gp_val_list, gp_label = [ 0.000, 0.000, 1.000], "gp"
    params_list = [j_val_list, k_val_list, g_val_list, gp_val_list]
    params_label_list = [j_label, k_label, g_label, gp_label]

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
sweep = HPhiSweeps(run, slurm_helper,
                   stan_cli_list,
                   params_list, params_label_list, hdir_norm,
                   cwd)
