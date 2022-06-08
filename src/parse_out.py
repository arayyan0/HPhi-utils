from lib_post import ParameterSpace, HPhiOutput, OneDParameterSweep, create_plot_filename
import glob
import os
import numpy as np
import matplotlib.pyplot as plt

################################################################################
################################################################################
################################################################################
which = '5.19.2022_24_eps_h_detailed'
number = 1
################################################################################
################################################################################
################################################################################

data_folder = f'archive/{which}/jobrun_{number}/'
print(data_folder)
plot_folder = data_folder + 'plots/'
if not os.path.exists(plot_folder):
    os.makedirs(plot_folder)

file_lst = sorted(glob.glob(data_folder+'*_/'))
paramslist, energieslist = [], []
for i, file in enumerate(file_lst):
    # print(file)

    param_block = file.split('/')[-2]
    split_file = param_block.split('_')[:-1]
    labels = [      split_file[0+2*i]  for i in range(0, int(len(split_file)/2))]
    params = [float(split_file[1+2*i]) for i in range(0, int(len(split_file)/2))]

    hphi_dir = HPhiOutput(file)
    if hphi_dir.EnergyQ == True:
        paramslist.append(params)
        energieslist.append(hphi_dir.Energies)

    ###plotting lanczos step, if you want it
    # fig = hphi_dir.plot_lanczos_step()
    # plt.savefig( create_plot_filename('lanczos_', plot_folder, param_block))
    # plt.close()

p_space = ParameterSpace(paramslist, labels, energieslist, plot_folder, data_folder)
