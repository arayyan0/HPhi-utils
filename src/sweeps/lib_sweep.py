from dataclasses import dataclass
import itertools as it
import numpy as np
import os
import stat
import sys

@dataclass
class ComputerInfo:
    hphi_build_loc: str
    mpiQ          : bool
    Ncorespernode : int

class ComputerPresets:
    computers = { }
    def add_computer(self, label, hphi_build_loc, mpiQ, Ncorespernode):
        self.computers[label] = ComputerInfo(hphi_build_loc, mpiQ,
                                             Ncorespernode)

    def __init__(self):
        hphi_laptop  = '/Users/ahmed/Documents/University/PhD/Research/General/HPhi/HPhi.build/'
        hphi_niagara = '/gpfs/fs0/scratch/h/hykee/arayyan/HPhi.build/'

        self.add_computer( 'laptop',  hphi_laptop, False,  4)
        self.add_computer('niagara', hphi_niagara,  True, 40)

@dataclass
class SLURMHelper:
    computer_settings: ComputerInfo
    Nnodes:            int
    hyperthreadQ:      str
    Ntasksperpoint:    int
    Ncpuspertask:      int
    time:              str

    def calculate_relevant_integers(self):
        Nthreadspercore = 2 if self.hyperthreadQ else 1
        #implement later: ensure these two are integers.
        Nthreadspernode = Nthreadspercore*self.computer_settings.Ncorespernode

        self.Ntaskspernode = 1
        self.Nj = 1

        # if (Nthreadspernode % self.Ncpuspertask == 0):
        #     self.Ntaskspernode = Nthreadspernode/self.Ncpuspertask
        # else:
        #     print('Your tasks per node is not an integer. Bye!')
        #     raise SystemExit
        #
        # if (self.Ntaskspernode % self.Ntasksperpoint == 0):
        #     self.Nj            = self.Ntaskspernode/self.Ntasksperpoint
        # else:
        #     print('Your points per node is not an integer. Bye!')
        #     raise SystemExit

    def create_local_sim_commands(self):
        self.preamble = f'mpiexec -np {self.Ntasksperpoint}' if self.computer_settings.mpiQ else ''
        self.postamble = '--wd $PWD ' if self.Nnodes > 1 else ''

    def create_submit_script_texts(self):
        slurm_settings = ''
        slurm_settings += f'#SBATCH --nodes={self.Nnodes}'+'\n'
        slurm_settings += f'#SBATCH --ntasks-per-node={int(self.Ntaskspernode)}'+'\n'
        slurm_settings += f'#SBATCH --cpus-per-task={self.Ncpuspertask}'+'\n'
        slurm_settings += f'#SBATCH --time={self.time}\n'
        self.slurm_settings = slurm_settings
        self.parallel_command = f'parallel --delay 0.5 -j {int(self.Nj)} '

class HPhiSweeps:
    def __init__(self, run, slurm_helper,
                       stan_cli_list,
                       params_list, params_label_list, hdirection,
                       cwd):

        self.PWD = cwd
        print(self.PWD)

        self.JobTitle = f"jobrun_{run}"
        print(self.JobTitle)

        self.HPhiBuild = slurm_helper.computer_settings.hphi_build_loc
        self.Preamble = slurm_helper.preamble
        self.NThreads = slurm_helper.Ntasksperpoint*slurm_helper.Ncpuspertask

        self.Labels = params_label_list
        self.Params = []
        self.create_parameters_arrays(params_list)

        for i in range(len(self.Labels)):
            print(f"{self.Labels[i]}:\n{self.Params[i]}")

        self.create_directory_structure(stan_cli_list, hdirection)
        self.create_genesis_script(slurm_helper)
        self.create_command_list()

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
        f.write(f'  command {self.Preamble} {self.HPhiBuild}src/HPhi $1 $2\n')
        f.write('}\n\n')

        f.write('HPhiDRY () {\n')
        f.write(f'  command {self.HPhiBuild}src/HPhi -sdry $1\n')
        f.write('}\n\n')

        # prepare cli and pipe to prep standard python script
        stan_cli_str = ' '.join(list(map(str, stan_cli_list)))

        param_str = ' '.join([f'{p:.12f}' for p in prod[:-1]])
        f.write(f'python3 '+ self.PWD+'/src/sweeps/prepare_standard_cli.py ' + stan_cli_str + ' ' + param_str + '\n\n')

        stan_str = '/'.join(filename.split('/')[:-1])+'/stan.in\n'

        f.write(f'HPhiDRY stan.in\n\n')

        append_cli_str = ' '.join([f'{hd:.12f}' for hd in hdirection])
        f.write(f'python3 '+ self.PWD+'/src/sweeps/append_cli.py ' f'{prod[-1]:.12f} ' + append_cli_str + '\n\n')

        f.write(f'export OMP_NUM_THREADS={self.NThreads}\n')
        f.write('HPhiSC -e namelist.def\n')
        f.write(f'export OMP_NUM_THREADS=1\n\n')

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

    def create_genesis_script(self, slurm_helper):
        f = open(f'{self.JobTitle}.sh','w+')
        f.write('#!/bin/bash\n\n')
        f.write(slurm_helper.slurm_settings)
        f.write(f'#SBATCH --job-name={self.JobTitle}\n\n')
        f.write('cd $SLURM_SUBMIT_DIR\n\n')
        f.write('module purge\n')
        f.write('module load python/3.9.8\n')      #for executing various scripts
        f.write('module load intel/2019u4 intelmpi/2019u4 scalapack\n') #for HPhi
        f.write('module load gnu-parallel\n\n')      #for the following
        f.write(slurm_helper.parallel_command+\
            f'--joblog {self.OutputPath}/{self.JobTitle}.out '  +\
            slurm_helper.postamble  +\
            f'< {self.JobTitle}.lst\n')
        f.close()

    def create_command_list(self):
        f = open(f'{self.JobTitle}.lst','w+')
        for dir in self.Folders:
            f.write(f"cd {dir} && ./local_sim.sh > sim_output.out && echo '{dir} finished'\n")
        f.close()
