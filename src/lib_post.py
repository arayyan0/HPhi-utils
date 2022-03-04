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
        #extract number of sites
        self.LocSpinFile = glob.glob(folder + 'locspn.def')[0]
        self.NSites = self.extract_number_sites()

        #extract lanczos step energies. Warning: full precision is not provided in this file.
        self.LanczosStepFile = glob.glob(folder + 'output/*Lanczos_Step*')[0]
        self.LanczosSteps, lanczos_energies = self.extract_lanczos_step()
        self.LanczosEnergies = lanczos_energies/self.NSites

        self.NStates = self.LanczosEnergies.shape[1]

        #extract accurate energies IF the simulation finished
        self.EnergyQ = True
        energy_file = glob.glob(folder + 'output/*energy.dat*')
        try:
            file = energy_file[0]
        except:
            self.EnergyQ = False
            print("Can't find energy file at " + folder)
        else:
            energies, self.GroundState = self.extract_energies(file)
            self.Energies        =        energies/self.NSites


    def extract_number_sites(self):
        '''
        extract the number of sites from the LocSpin File
        '''
        with open(self.LocSpinFile, 'r') as f:
            content = f.readlines()
        num_sites=int(content[1].replace('\n',' ').split(' ')[-4])
        return num_sites

    def extract_energies(self, file):
        '''
        extract the energies from the Energy file
        '''
        with open(file, 'r') as f:
            content = f.readlines()

        line_number_lst = []
        for line_number, file_line in enumerate(content):
            z = re.match('  Energy', file_line)
            if z:
                line_number_lst.append(line_number)

        energies = [np.array(content)[line_number_lst][i].split(' ')[4] for i in range(self.NStates)]
        energy_arrays = np.array(list(map(float, energies)))

        return energy_arrays, np.min(energy_arrays)

    def extract_lanczos_step(self):
        '''
        extract the energies from the LanczosStep file
        '''
        with open(self.LanczosStepFile, 'r') as f:
            content = f.readlines()[1:]

        content_num = [line.replace('\n',' ').split() for line in content]
        content_arr = np.array(content_num)

        steps    = np.array(content_arr[:, 0 ], np.int  )
        energies = np.array(content_arr[:, 3:], np.float)
        return steps, energies

    def plot_lanczos_step(self):
        subtracted_energies = self.LanczosEnergies - np.tile(self.LanczosEnergies[:,0], (self.NStates, 1)).T

        # plot of states over iterations
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        fig, ax = plt.subplots()
        for i, color in zip(range(self.NStates), colors):
            ax.scatter(self.LanczosSteps, subtracted_energies[:,i], marker=".",clip_on=False,
                       s=20, facecolors='none', edgecolors=color, linewidth=1.5, label=f'state {i}')
            ax.plot(subtracted_energies[:,i], ls='--')
        plt.legend()
        ax.set_xlabel('Lanczos step')
        ax.set_ylabel(r'$\frac{E-E_0}{N}$')
        ax.set_ylim(0,np.max(subtracted_energies))

def is_not_unique(params):
    bool_array = []
    num_params = params.shape[1]
    for i in range(num_params):
        tol = 10**-6
        a = params[:,i]
        bool_array.append( ~(np.abs(a[0] - a) < tol).all() )
    return bool_array

@dataclass
class EnergyDerivatives:
    paramlist:   np.ndarray
    elist:       np.ndarray
    factor:      float

    Colors = ["blue", "magenta", "green"] #nondark background

    def calculate_derivs(self):
        self.m = -np.gradient(self.elist, self.paramlist, edge_order=2) / self.factor
        self.chi = np.gradient(self.m, self.paramlist, edge_order=2) / self.factor

    def plot_energy_derivs(self):
        self.calculate_derivs()
        functions = [self.elist, self.m, self.chi]
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        axes = [ax1, ax2, ax2.twinx()]
        colors=["blue", "magenta", "green"]

        for function, ax, color in zip(functions, axes, colors):
            ax.scatter(
                self.paramlist,
                function,
                marker=".",
                # clip_on=False,
                s=20,
                facecolors='none',
                edgecolors=color,
                linewidth=1.5)
            ax.tick_params(axis="y", colors=color)

        ax1.grid(True, axis='x')
        ax2.grid(True, axis='x')

        plt.xlim(min(self.paramlist), max(self.paramlist))
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        return fig

    def return_labels(self, TeXlabel_for_derivs):
        ylabel = [
        r"$\frac{E_0}{N}$",
        r"$-\frac{1}{N}\frac{\mathrm{d}E_0}{\mathrm{d} %s }$" % (TeXlabel_for_derivs),
        r"$-\frac{1}{N}\frac{\mathrm{d}^2E_0}{\mathrm{d}%s^2}$" % (TeXlabel_for_derivs)
        ]
        return ylabel

def create_plot_filename(plot_type, plot_folder, param_block):
    str = plot_folder + plot_type + ''.join(param_block) + '.pdf'
    return str

@dataclass
class ParameterInfo:
    TeXlabel: str
    angleQ:   bool

    def add_pi_if_angle(self):
        f   =      pi if self.angleQ else 1
        str = r'/\pi' if self.angleQ else ''
        self.factor    = f
        self.TeXlabel_for_derivs  = self.TeXlabel
        self.TeXlabel += str

class ParameterPresets:
    parameters = { }
    def add_parameter(self, label, TeXlabel, angleQ):
        x = ParameterInfo(TeXlabel, angleQ)
        x.add_pi_if_angle()
        self.parameters[label] = x

    def __init__(self):
        #j-k-g-gp
        self.add_parameter(  'j',         'J', False)
        self.add_parameter(  'k',         'K', False)
        self.add_parameter(  'g',   r'\Gamma', False)
        self.add_parameter( 'gp', r"\Gamma'",  False)
        #multipolar_1
        self.add_parameter(  't',   r'\theta',  True)
        self.add_parameter(  'p',     r'\phi',  True)
        self.add_parameter( 'jb',      r'J_B', False)
        #multipolar_2
        self.add_parameter( 'xi',      r'\xi',  True)
        self.add_parameter('eps', r'\epsilon', False)
        #magnetic field
        self.add_parameter(  'h',         'h', False)

class OneDParameterSweep:
    def __init__(self, paramslist, energieslist, paramlabel, which_parameter_to_sort):
        self.params, self.energies= list(map(np.array, [paramslist, energieslist]))
        self.param_labels = paramlabel

        self.sort_parameters(which_parameter_to_sort)

        self.numstates = self.energies.shape[1]
        # decide what to plot based on this value.

    def sort_parameters(self, which_parameter_to_sort):
        self.swept_param_info = ParameterPresets().parameters[which_parameter_to_sort]

        self.swept_param_index = self.param_labels.index(which_parameter_to_sort)

        idx          = np.argsort(self.params[:, self.swept_param_index])
        self.params   = self.params[idx,:]
        self.energies = self.energies[idx]

    def plot_spectrum(self, ylim):
        subtracted_energies = self.energies - np.tile(self.energies[:,0], (self.numstates, 1)).T
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        fig, ax = plt.subplots()
        for i, color in zip(range(self.numstates),colors):
            ax.scatter(self.params[:,self.swept_param_index], subtracted_energies[:,i],
                        marker=".",
                        # clip_on=False,
                        s=20,
                        facecolors='none',
                        edgecolors=color,
                        label=f'state {i}',
                        linewidth=1.5)
            ax.plot(self.params[:,self.swept_param_index], subtracted_energies[:,i], ls='--')#, color=c)
        ax.set_xlabel('$'+self.swept_param_info.TeXlabel+'$')
        ax.set_ylabel(r'$\frac{E-E_0}{N}$')
        ax.set_ylim(0,ylim)
        plt.legend()
        return fig

    def plot_gs_properties(self):
        e_deriv = EnergyDerivatives(self.params[:,self.swept_param_index],
                                    self.energies[:,0],
                                    self.swept_param_info.factor)
        labels = e_deriv.return_labels(self.swept_param_info.TeXlabel_for_derivs)
        fig = e_deriv.plot_energy_derivs()
        for ax, label, color in zip(fig.get_axes(),labels, EnergyDerivatives.Colors):
            ax.set_ylabel(label, rotation = "horizontal",fontsize=12,labelpad=20,color=color)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        fig.get_axes()[1].set_xlabel('$'+self.swept_param_info.TeXlabel+'$')
        return fig
