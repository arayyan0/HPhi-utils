import sys

# code should be run from within HPhi-utils
filename = 'speedrun.sh'
HPhi_build = '/Users/ahmed/Documents/University/PhD/Research/General/HPhi/HPhi.build/'
project = 'kga'
output_folder = 'out'
input_file = 'stan.in'

what_simulation = sys.argv[1] # can be 'simple' or 'SSFsimple'

# writing bash script
f = open(filename, 'w')
f.write('#!/bin/bash\n')
f.write('\n')
# writing function shortcut for HPhi
f.write('mode=$1\n')
f.write('HPhiSC () {\n')
f.write(f'  command {HPhi_build}src/HPhi $mode $1\n')
f.write('}\n')
f.write('\n')
# writing function shortcut for HPhi
f.write(f'vim src/prepare_standard.py\n')
f.write(f'python3 src/prepare_standard.py\n')
f.write('\n')
#get into output folder, then figure out what to do
f.write(f'cd {output_folder}\n')
if what_simulation == 'simple':
    #run
    f.write(f'HPhiSC ' + input_file + '\n')
    #get out back into HPhi-utils
    f.write(f'cd ..\n')
    f.close()
elif what_simulation == 'SSFsimple':
    #run
    f.write(f'HPhiSC ' + input_file + '\n')
    #append to geometry.dat
    f.write(f'python3 ../src/ssf_post.py\n')
    #run greenr2k to calculate reciprocal lattice properties
    f.write(f'{HPhi_build}tool/greenr2k namelist.def geometry.dat\n')
    #create gp_script.gp
    g = open(f'{output_folder}/gp_script.gp', 'w')
    g.write('load "kpath.gp"\n')
    g.write('plot "output/zvo_corr_eigen0.dat" u 1:12 w l\n')
    g.write('pause -1 "Hit any key to continue"\n')
    g.close()
    #plot the kpath
    f.write(f'gnuplot gp_script.gp\n')
    #get out back into HPhi-utils
    f.write(f'cd ..\n')
    f.close()
