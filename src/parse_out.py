from lib_post import HPhiOutput, is_not_unique, FreeEnergyDerivatives, pi
import glob
import os
import numpy as np
import matplotlib.pyplot as plt

which = '2.24.2022-sweeptests'

number = 1
data_folder = f'archive/{which}/jobrun_{number}/'
print(data_folder)
plot_folder = data_folder + f'plots/'
if not os.path.exists(plot_folder):
    os.makedirs(plot_folder)

file_lst = sorted(glob.glob(data_folder+'*_/'))

paramslist, energieslist = [], []
for file in file_lst:
    split_file = file.replace('/','_').split('_')[4:-2] #should just be labels and values

    labels = [split_file[0+2*i] for i in range(0, int(len(split_file)/2))]
    params_str = [split_file[1+2*i] for i in range(0, int(len(split_file)/2))]

    params = list(map(float, params_str))

    sim_directory = HPhiOutput(file)

    paramslist.append(params)
    energieslist.append(sim_directory.Energies)

paramslist = np.array(paramslist)
which_params = is_not_unique(paramslist)
energieslist = np.array(energieslist)
################################################################################
################################################################################
################################################################################
which_parameter_to_sort = 0
################################################################################
################################################################################
################################################################################
idx = np.argsort(paramslist[:, which_parameter_to_sort])
paramslist = paramslist[idx,:]
energieslist = energieslist[idx]

numstates = energieslist.shape[1]

# spectrum plot
fig, ax = plt.subplots()
for i in range(numstates):
    ax.scatter(paramslist[:,which_parameter_to_sort], energieslist[:,i],
                marker="o",
                # clip_on=False,
                s=20,
                facecolors='none',
                edgecolors='k',
                linewidth=1.5)
    ax.plot(paramslist[:,which_parameter_to_sort], energieslist[:,i], ls='--')#, color=c)
# plt.show()
# plt.savefig()
plt.close()

derivs = FreeEnergyDerivatives(paramslist[:,which_parameter_to_sort], energieslist[:,0], pi)
colors, color_order = derivs.Colors, [0,2,1]
colors = [colors[color_index] for color_index in color_order]

# # figure out the y-labels
# ylabel = [
# r"$\frac{E_0}{N}$",
# r"$-\frac{1}{N}\frac{\mathrm{d}^2E_0}{\mathrm{d}%s^2}\quad$" % (paramsTeXlabel[i]),
# r"$-\frac{1}{N}\frac{\mathrm{d}E_0}{\mathrm{d}%s}$" % (paramsTeXlabel[i]),
# ]

# plot the derivs with proper labels
fig = derivs.PlotSweep()
# for j, [ax, color] in enumerate(zip(fig.axes, colors)):
    # ax.set_ylabel(     ylabel[j], rotation = "horizontal",
    #                  fontsize=12,             labelpad=20,
    #                  color=color)
# fig.axes[1].set_xlabel(r"$%s$ " % xlabel )
fig.tight_layout(rect=[0, 0.03, 1, 0.95])

# save and show the plot after creating filename
# s = ''
# for k, boolean in enumerate(which_params):
    # if ~boolean:
        # s = s + paramslabel[k] + '_' + f'{df[paramslabel[k]][0]:.6f}_'
# print(s)
# plt.savefig(plot_folder + 'energy_' + s + '.pdf')
plt.show()
plt.close()
