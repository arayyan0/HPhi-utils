import numpy as np
import matplotlib.pyplot as plt
################################################################################
################################################################################
################################################################################
which = '3.19.2022_24_eps_h'
number = 1
################################################################################
################################################################################
################################################################################
min_prom=0
data_folder = f'archive/{which}/jobrun_{number}/'
data_file = data_folder+f'prom_{min_prom:.3f}.txt'
data = np.loadtxt(data_file)

#scale magnetic field to the proper limit
data[:, 1] = data[:,1]*2/3/(1/2)

#remove points greater than h_max
h_max = 0.5
data = data[data[:,1]<=h_max+0.05]

#remove points with min_prom less than some value
min_prom_max=0.19
data = data[data[:,2]>min_prom_max]

#test out the data...does it look good?
fig, ax = plt.subplots()
plt.scatter(data[:,0], data[:,1], s=data[:,2])
plt.xlim((-1,1))
plt.ylim((0,h_max))
plt.show()
plt.close()

#filtered data
filtered_data_file = data_folder+f'filtered_prom_{min_prom:.3f}.out'

# np.savetxt(filtered_data_file, data[:, :2],fmt='%.12f')
