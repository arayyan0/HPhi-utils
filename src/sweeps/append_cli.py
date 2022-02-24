import sys
import numpy as np

path_of_lib_pre = '/'.join(sys.argv[0].split('/')[:-2])
sys.path.append(path_of_lib_pre)
from lib_pre import StandardInputAux

loc_file = 'locspn.def'
trans_file = 'trans.def'

hstrength = float(sys.argv[1])
hdirection = np.array(list(map(float,sys.argv[2:5])))

#--------------------read necessary info and modify the trans file
aux_simulation = StandardInputAux(loc_file, trans_file,
                                  hstrength, hdirection)
aux_simulation.extract_num_sites()
aux_simulation.add_to_trans()
