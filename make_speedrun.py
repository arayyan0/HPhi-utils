import sys
import os
# code should be run from within HPhi-utils
filename = 'speedrun.sh'
what_simulation = 'simple'  #(either 'SSF' or 'simple')
what_computer  = 'niagara' #(either 'laptop', 'home', or 'niagara')

if what_computer == 'laptop':
    HPhi_build = '/Users/ahmed/Documents/University/PhD/Research/General/HPhi/HPhi.build/'
    run_preamble = ''
elif what_computer == 'home':
    #HPhi_build = 
    #run_preamble =     
    pass
elif what_computer == 'niagara':
    HPhi_build = '/scratch/h/hykee/arayyan/HPhi.build/'
    num_processes = 4
    run_preamble = f'mpiexec -np {num_processes}'

output_folder = 'out'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

#-------------------------------beginning of bash script
# preamble of bash scripts
f = open(filename, 'w')
f.write('#!/bin/bash\n')
f.write('\n')

# writing function shortcut for HPhi
f.write('HPhiSC () {\n')
f.write(f'  command {run_preamble} {HPhi_build}src/HPhi $1 $2\n')
f.write('}\n')
f.write('\n')

# edit prepare standard file
f.write('vim src/prepare_standard.py\n')
f.write('python3 src/prepare_standard.py\n')
f.write('\n')

#get into output folder
f.write(f'cd {output_folder}\n')
#run in dry standard mode
f.write('HPhiSC -sdry stan.in\n')
#get out back into HPhi-utils
f.write('cd ..\n')
f.write('\n')

# append to standard file
f.write('vim src/append_to_standard.py\n')
f.write('python3 src/append_to_standard.py\n')
f.write('\n')

#get into output folder
f.write(f'cd {output_folder}\n')
#run in expert mode
f.write('HPhiSC -e namelist.def\n')
f.write('\n')

if what_simulation == 'SSF':
    #append to geometry.dat
    f.write('python3 ../src/ssf_post.py\n')
    #run greenr2k to calculate reciprocal lattice properties
    f.write(f'{HPhi_build}tool/greenr2k namelist.def geometry.dat\n')
    f.write('\n')
    #create gp_script.gp
    g = open(output_folder+'/gp_script.gp', 'w')
    g.write('load "kpath.gp"\n')
    g.write('plot "output/zvo_corr_eigen0.dat" u 1:12 w l\n')
    g.write('pause -1 "Hit any key to continue"\n')
    g.close()

    if what_computer != 'niagara':
        #plot the kpath
        f.write('gnuplot gp_script.gp\n')
        f.write('\n')

#get out back into HPhi-utils
f.write('cd ..\n')
f.write('\n')
f.close()
