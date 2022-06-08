import numpy as np
import matplotlib.pyplot as plt
################################################################################
################################################################################
################################################################################
which = '5.19.2022_24_eps_h_detailed'
number = 1
min_prom=0
data_folder = f'archive/{which}/jobrun_{number}/'
data_file = data_folder+f'prom_{min_prom:.3f}_final.txt'
################################################################################
################################################################################
################################################################################
data = np.loadtxt(data_file)


#scale magnetic field to the proper limit
momentS = 1
scale   = 1
data[:, 1] = data[:, 1]/momentS/scale
# data[:, 1] = data[:, 1]

#remove points greater than h_max
h_max = 0.6
data = data[data[:,1]<=h_max+0.05]

#remove points with min_prom less than some value
min_prom_max=0.19
data = data[data[:,2]>min_prom_max]

#test out the data...does it look good?
fig, ax = plt.subplots()
plt.scatter(data[:,0], data[:,1], s=data[:,2])
plt.xlim((-1.0,1.0))
plt.ylim((0,h_max))
plt.show()
plt.close()

#filtered data
filtered_data_file = data_folder+f'filtered_prom_{min_prom:.3f}.out'

#sort it by x axis values
print(data[np.argsort(data[:,0],axis=0),:])

sorted_data = data[np.argsort(data[:,0],axis=0),:]


np.savetxt(filtered_data_file, sorted_data[:, :2],fmt='%.12f')
