from lib_post import HPhiOutput, OneDParameterSweep, create_plot_filename
import glob
import os
import numpy as np
import matplotlib.pyplot as plt

################################################################################
################################################################################
################################################################################
which = '3.02.2022_18_xi_spectrum'
number = 1
which_parameter_to_sort = 'xi'
################################################################################
################################################################################
################################################################################

data_folder = f'archive/{which}/jobrun_{number}/'
print(data_folder)
plot_folder = data_folder + f'plots/'
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

sweep = OneDParameterSweep(paramslist,energieslist,labels,which_parameter_to_sort)

# spectrum plot
# if sweep.numstates > 1:
#     ylim = 0.01
#     fig = sweep.plot_spectrum(ylim)
#     plt.savefig(create_plot_filename('spectrum_', plot_folder, ''))
#     plt.show()
#     plt.close()

#gs and derivs plot
# fig = sweep.plot_gs_properties()
# plt.savefig(create_plot_filename('gsenergy_', plot_folder, ''))
# plt.show()
# plt.close()
