import sys

# code should be run from within HPhi-utils
filename = 'speedrun.sh'
HPhi_build = '/Users/ahmed/Documents/University/PhD/Research/General/HPhi/HPhi.build/'
output_folder = 'out'

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

# edit append to standard file
f.write(f'vim src/append_to_standard.py\n')
f.write(f'python3 src/append_to_standard.py\n')
f.write('\n')

# #get into output folder
# f.write(f'cd {output_folder}\n')
# #run in expert mode
# f.write(f'HPhiSC -e namelist.def\n')
# #get out back into HPhi-utils
# f.write(f'cd ..\n')
# f.write('\n')
f.close()
