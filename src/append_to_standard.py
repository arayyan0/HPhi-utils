from lib_pre import StandardInputAux, pi
from math import sin, cos
import numpy as np
import os

output_path = f'out/'
loc_file = output_path + 'locspn.def'
trans_file = output_path + 'trans.def'
#--------------------parameters for magnetic field
X, Y, Z = np.eye(3)

hstrength      = 0
htheta         = 0 #(in degrees, pointing from Y direction towards Z direction)
hdirection     = cos(htheta/180*pi)*Y + sin(htheta/180*pi)*Z
#--------------------read necessary info and modify the trans file
aux_simulation = StandardInputAux(loc_file, trans_file,
                                  hstrength, hdirection)
aux_simulation.extract_num_sites()
aux_simulation.add_to_trans()
