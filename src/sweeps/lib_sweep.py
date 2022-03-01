import itertools as it
import numpy as np
import os
import stat
import sys

class HPhiSweeps:
    def __init__(self, run, what_computer,
                       hpc_settings,
                       stan_cli_list,
                       params_list, params_label_list, hdirection,
                       cwd):

        self.PWD = cwd

        self.JobTitle = f"jobrun_{run}"
        print(self.JobTitle)

        [self.NProc, self.NTasks, self.NNodes, self.Time] = hpc_settings

        self.WhatComputer = what_computer
        if self.WhatComputer == 'laptop':
            self.HPhiBuild = '/Users/ahmed/Documents/University/PhD/Research/'+\
                         'General/HPhi/HPhi.build/'
            self.RunPreamble = ''
        elif self.WhatComputer == 'home':
            #to be implemented later
            #self.HPhiBuild =
            #self.RunPreamble =
            pass
        elif self.WhatComputer == 'niagara':
            self.HPhiBuild = '/scratch/h/hykee/arayyan/HPhi.build/'
            self.RunPreamble = f'mpiexec -np {self.NProc}'

        self.Labels = params_label_list
        self.Params = []
        self.create_parameters_arrays(params_list)

        for i in range(len(self.Labels)):
            print(f"{self.Labels[i]}:\n{self.Params[i]}")

        self.create_directory_structure(stan_cli_list, hdirection)

        self.create_command_list()

        self.create_genesis_script()

    def create_parameters_arrays(self, params_list):
        for params in params_list:
            x_min, x_max, dx = params
            if abs(x_max - x_min) < pow(10,-8):
                self.Params.append(np.array([x_min]))
            else:
                N_x = round(1+ (x_max-x_min)/dx)
                self.Params.append(np.linspace(x_min, x_max, N_x))

    def create_directory_structure(self, stan_cli_list, hdirection):
        self.OutputPath = f"out/{self.JobTitle}"
        if not os.path.exists(self.OutputPath):
            os.makedirs(self.OutputPath)

        a = [list(param) for param in self.Params]
        product =   list(it.product(*a))

        self.Folders = []
        self.NPoints = len(product)
        for prod in product:
            param_label_list = ''.join(map(str, [f'{self.Labels[i]}_{prod[i]:.12f}_' for i in range(len(self.Params))]))
            folder_name = self.OutputPath+'/'+param_label_list+'/'
            if not os.path.exists(folder_name):
                os.makedirs(folder_name)
            self.Folders.append(folder_name)

            sim_filename = folder_name+'local_sim.sh'
            self.create_each_script(sim_filename, stan_cli_list, prod, hdirection)



    def create_each_script(self, filename, stan_cli_list, prod, hdirection):
        '''
        creates bash script that
        1) creates standard input file
        2) performs a dry run
        3) append to standard file
        4) executes expert mode
        5) creates geometry.dat and runs greenr2k
        6) create gp_script.gp
        '''

        # preamble of bash scripts
        f = open(filename, 'w')
        #change permissions to execute by user, group, others
        permissions = stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH | stat.S_IXGRP | stat.S_IXOTH
        os.chmod(filename, permissions)
        f.write('#!/bin/bash\n\n')

        # writing function shortcut for HPhi
        f.write('HPhiSC () {\n')
        f.write(f'  command {self.RunPreamble} {self.HPhiBuild}src/HPhi $1 $2\n')
        f.write('}\n\n')

        # prepare cli and pipe to prep standard python script
        stan_cli_str = ' '.join(list(map(str, stan_cli_list)))

        param_str = ' '.join([f'{p:.12f}' for p in prod[:-1]])
        f.write(f'python3 '+ self.PWD+'/src/sweeps/prepare_standard_cli.py ' + stan_cli_str + ' ' + param_str + '\n\n')

        stan_str = '/'.join(filename.split('/')[:-1])+'/stan.in\n'

        f.write(f'HPhiSC -sdry stan.in\n\n')

        append_cli_str = ' '.join([f'{hd:.12f}' for hd in hdirection])
        f.write(f'python3 '+ self.PWD+'/src/sweeps/append_cli.py ' f'{prod[-1]:.12f} ' + append_cli_str + '\n\n')

        f.write('HPhiSC -e namelist.def\n\n')

        output_mode = stan_cli_list[-2]
        if output_mode == 'correlation':
            #-append to geometry.dat
            f.write(f'python3 '+ self.PWD+'/src/ssf_post.py geometry.dat\n')
            #-run greenr2k to calculate reciprocal lattice properties
            f.write(f'{self.HPhiBuild}tool/greenr2k namelist.def geometry.dat\n')
            #-create gp_script.gp
            gp_script_str = '/'.join(filename.split('/')[:-1])+'/gp_script.gp'
            g = open(gp_script_str, 'w')
            g.write('load "kpath.gp"\n')
            g.write('plot "output/zvo_corr_eigen0.dat" u 1:12 w l\n')
            g.write('pause -1 "Hit any key to continue"\n')
            g.close()

        f.close()

    def create_genesis_script(self):

        self.NBatch = self.NNodes*self.NTasks/self.NProc
        if not (self.NNodes*self.NTasks % self.NProc == 0):
            print( 'Watch out. Number of parallel jobs (or batches) is not an integer. '+\
                  f'Do NOT submit {self.JobTitle}.sh!')

        f = open(f'{self.JobTitle}.sh','w+')
        f.write('#!/bin/bash\n\n')

        f.write(f'#SBATCH --nodes={self.NNodes}'+'\n')
        f.write(f'#SBATCH --ntasks-per-node={self.NTasks}'+'\n')
        f.write(f'#SBATCH --time={self.Time}\n')
        f.write(f'#SBATCH --job-name={self.JobTitle}\n\n')
        f.write('cd $SLURM_SUBMIT_DIR\n\n')
        f.write('module purge\n')
        f.write('module load python/3.9.8\n')      #for executing various scripts
        f.write('module load intel/2019u4 intelmpi/2019u4 scalapack\n') #for HPhi
        f.write('module load gnu-parallel\n\n')      #for the following
        f.write(f'parallel -j {int(self.NBatch)} --joblog {self.OutputPath}/{self.JobTitle}.out --wd $PWD'+
                 f' < {self.JobTitle}.lst\n')
        f.close()

    def create_command_list(self):
        f = open(f'{self.JobTitle}.lst','w+')
        for dir in self.Folders:
            f.write(f'cd {dir} && ./local_sim.sh > sim_output.out\n')
        f.close()
