import numpy as np
from dataclasses import dataclass
from math import sqrt
import os
import glob
import matplotlib.pyplot as plt
import re

pi = np.pi

class TriangularSymmetryPoints:
    symmetry_points = { }
    def add_symmetrypoint(self, label, vector):
        self.symmetry_points[label] = self.conventional_unit_cell @ vector

    def __init__(self, conventional_unit_cell, path):
        self.path = path
        self.conventional_unit_cell = conventional_unit_cell
        #positions are given in the primitive unit cell enclosing HC z-bond
        self.add_symmetrypoint( 'X', np.array([-1/2, 1/2, 0]))
        self.add_symmetrypoint( 'K', np.array([-1/3, 1/3, 0]))
        self.add_symmetrypoint( 'G', np.array([0, 0, 0]))
        self.add_symmetrypoint('M2', np.array([ 1/2, 1/2, 0]))
        self.add_symmetrypoint('Kp', np.array([2/3, 1/3, 0]))
        self.add_symmetrypoint('Gp', np.array([1, 0, 0]))
        self.add_symmetrypoint('M1', np.array([1/2, 0, 0]))

    def create_path(self):
        texts = [point + ' ' +
                  ' '.join(self.symmetry_points[point].round(3).astype(str))+'\n'
                 for point in self.path]
        text = ''.join(texts)
        return text #removes last newline character

    def append_to_geometry(self, filename, line_density, surface_density):
        existsQ = os.path.isfile(filename)
        if existsQ == False:
            print('geometry.dat does not exist!')
            raise SystemExit

        first_line = f'{len(self.path)} {line_density}\n'
        body = self.create_path()
        last_line  = f'{surface_density} {surface_density} 1\n'
        text = first_line + body + last_line

        f = open(filename, 'a')
        f.write(text)
        f.close()

class HPhiOutput:
    def __init__(self, folder):
        lanczos_step_file = glob.glob(folder + 'output/*Lanczos_Step*')
        energy_file = glob.glob(folder + 'output/*energy.dat*')
        self.plot_lanczos_step(lanczos_step_file[0])
        self.extract_energies(energy_file[0])
        pass

    def extract_lanczos_step(self, file):
        with open(file, 'r') as f:
            content = f.readlines()[1:]

        content_num = [line.replace('\n',' ').split() for line in content]
        content_arr = np.array(content_num)

        steps    = np.array(content_arr[:, 0 ], np.int  )
        energies = np.array(content_arr[:, 3:], np.float)
        return steps, energies

    def plot_lanczos_step(self, file):
        steps, energies = self.extract_lanczos_step(file)

        states_num = energies.shape[1]

        # plot of states over iterations
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        fig, ax = plt.subplots()
        for i, color in zip(range(states_num), colors):
            # c=np.random.rand(3,)
            ax.scatter(steps, energies[:,i], marker="o",
                                        # clip_on=False,
                                        s=20,
                                        facecolors='none',
                                        edgecolors=color,
                                        linewidth=1.5,
                                        label=f'{i}th state')
            ax.plot(steps, energies[:,i], ls='--')#, color=c)
        plt.legend()
        # plt.show()
        # plt.savefig()
        plt.close()

    def extract_energies(self, file):
        with open(file, 'r') as f:
            content = f.readlines()

        line_number_lst = []
        for line_number, file_line in enumerate(content):
            z = re.match('  Energy', file_line)
            if z:
                line_number_lst.append(line_number)

        num_states = len(line_number_lst)
        energies = [np.array(content)[line_number_lst][i].split(' ')[4] for i in range(num_states)]

        self.Energies = np.array(list(map(float, energies)))
        self.GroundState = np.min(self.Energies)

def is_not_unique(params):
    bool_array = []
    num_params = params.shape[1]
    for i in range(num_params):
        tol = 10**-6
        a = params[:,i]
        bool_array.append( ~(np.abs(a[0] - a) < tol).all() )
    return bool_array

class FreeEnergyDerivatives:
    Colors = ["blue", "magenta", "green"] #nondark background
    # Colors = ["turquoise", "limegreen", "orange"] #dark background
    # Colors = ["turquoise", "limegreen", "orange", "red"] #dark background

    def __init__(self, x_list, y_list, factor):
        self.XList = x_list
        self.YList = y_list
        self.Factor = factor

    def PseudoMagnetization(self):
        m = -np.gradient(self.YList, self.XList, edge_order=2) / self.Factor
        return m

    def PseudoSusceptibility(self):
        m = self.PseudoMagnetization()
        chi = np.gradient(m, self.XList, edge_order=2) / self.Factor
        return chi

    def ThirdDerivative(self):
        chi = self.PseudoSusceptibility()
        f = np.gradient(chi, self.XList, edge_order=2) / self.Factor
        return f

    def PlotSweep(self):
        m = self.PseudoMagnetization()
        chi = self.PseudoSusceptibility()
        # f = self.ThirdDerivative()

        functions = [self.YList, m, chi]
        # functions = [self.YList, chi, m, f]
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        # fig.subplots_adjust(top=1.5)
        axes = [ax1, ax1.twinx(), ax2]
        # axes = [ax1, ax2, ax1.twinx(), ax2.twinx()]

        for function, ax, color in zip(functions, axes, self.Colors):
            # print(function, ax, color)
            ax.scatter(
                self.XList,
                function,
                marker=".",
                # clip_on=False,
                s=20,
                facecolors='none',
                edgecolors=color,
                linewidth=1.5)
            ax.tick_params(axis="y", colors=color)
        # axes[2].axhline(c='gray',ls="-.")
        # ax2.axhline(color=self.Colors[1], ls="-.")
        # ax2.set_ylim([-0.25,1.25])
        # axes[2].set_ylim([-0.25,1.25])

        # ax2.set_ylim([-10,10])

        ax1.grid(True, axis='x')
        ax2.grid(True, axis='x')

        plt.xlim(min(self.XList), max(self.XList))
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])

        return fig

    def PseudoSusceptibilityPeaks(self, prom):
        f = self.PseudoSusceptibility()

        x_peak_list, f_peak_list = [], []
        f_peaks, f_prominences = find_peaks(f, prominence=prom)

        for f_peak_index in f_peaks:
            x_peak_list.append(self.XList[f_peak_index])
            f_peak_list.append(f[f_peak_index])

        return x_peak_list, f_peak_list, f_prominences["prominences"]
