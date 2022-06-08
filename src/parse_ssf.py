'''
this code calculates correlation functions for HPhi output

step 1
read in HPhi output
'''
from lib_post import HPhiOutput, TriangularReciprocalSpaceGrid
import glob
import os
################################################################################
################################################################################
################################################################################
which = '6.01.2022_24_eps_h_ssf'
number = 2
################################################################################
################################################################################
################################################################################
data_folder = f'archive/{which}/jobrun_{number}/'
# print(data_folder)
plot_folder = data_folder + 'plots/'
if not os.path.exists(plot_folder):
    os.makedirs(plot_folder)

file_lst = sorted(glob.glob(data_folder+'*_/'))
paramslist = []

density=36
size =2
grid = TriangularReciprocalSpaceGrid(density,size)
#fraction, orientation, colormap
cb_options = [0.1, 'vertical', 'inferno']

for i, file in enumerate(file_lst):
    print(file)

    param_block = file.split('/')[-2]
    split_file = param_block.split('_')[:-1]
    labels = [      split_file[0+2*i]  for i in range(0, int(len(split_file)/2))]
    params = [float(split_file[1+2*i]) for i in range(0, int(len(split_file)/2))]

    hphi_dir = HPhiOutput(file)
    if hphi_dir.CorrOutputQ == True:
        paramslist.append(params)
        hphi_dir.create_corr_funcs()
        hphi_dir.create_recip_corr_funcs(grid)
        hphi_dir.plot_ssf(grid, cb_options,plot_folder,i)
