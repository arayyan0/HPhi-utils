import sys
import os
# code should be run from within HPhi-utils
filename = 'speedrun.sh'
#HPhi_build = '/Users/ahmed/Documents/University/PhD/Research/General/HPhi/HPhi.build/'
HPhi_build = '/scratch/h/hykee/arayyan/HPhi.build/'
output_folder = 'out'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
what_simulation = 'SSF'  #(either 'SSF' or 'simple')
what_computer  = 'home' #(either 'home' or 'cluster')

#-------------------------------beginning of bash script
# preamble of bash scripts
f = open(filename, 'w')
f.write('#!/bin/bash\n')
f.write('\n')

# writing function shortcut for HPhi
f.write('HPhiSC () {\n')
f.write(f'  command {HPhi_build}src/HPhi $1 $2\n')
f.write('}\n')
f.write('\n')

# edit prepare standard file
f.write(f'vim src/prepare_standard.py\n')
f.write(f'python3 src/prepare_standard.py\n')
f.write('\n')

#get into output folder
f.write(f'cd {output_folder}\n')
#run in dry standard mode
f.write(f'HPhiSC -sdry stan.in\n')
#get out back into HPhi-utils
f.write(f'cd ..\n')
f.write('\n')

# append to standard file
f.write(f'vim src/append_to_standard.py\n')
f.write(f'python3 src/append_to_standard.py\n')
f.write('\n')

#get into output folder
f.write(f'cd {output_folder}\n')
#run in expert mode
f.write(f'HPhiSC -e namelist.def\n')
f.write('\n')

if what_simulation == 'SSF':
    #append to geometry.dat
    f.write(f'python3 ../src/ssf_post.py\n')
    #run greenr2k to calculate reciprocal lattice properties
    f.write(f'{HPhi_build}tool/greenr2k namelist.def geometry.dat\n')
    f.write('\n')
    #create gp_script.gp
    g = open(output_folder+'/gp_script.gp', 'w')
    g.write('load "kpath.gp"\n')
    g.write('plot "output/zvo_corr_eigen0.dat" u 1:12 w l\n')
    g.write('pause -1 "Hit any key to continue"\n')
    g.close()

    if what_computer == 'home':
        #plot the kpath
        f.write(f'gnuplot gp_script.gp\n')
        f.write('\n')

#get out back into HPhi-utils
f.write(f'cd ..\n')
f.write('\n')
f.close()
