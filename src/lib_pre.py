from dataclasses import dataclass
from math import sin, cos, sqrt
import numpy as np

pi = np.pi

class UnitCellPresets:
    #unit cell shape presets
    unit_cells = { }
    def add_unit_cell(self, label, A0, A1):
        self.unit_cells[label] = [A0, A1]

    def __init__(self):
        a0, a1 = np.eye(2)
        self.add_unit_cell(   'rect',   a0 - a1,   a0 + a1)
        self.add_unit_cell('rhom120', 2*a0 - a1, 2*a1 - a0)
        self.add_unit_cell( 'rhom60',   a0 - a1,        a0)

@dataclass
class HoneycombInfo:
    unitcell_shape: str
    W:              int
    L:              int

    def create_lattice_vectors(self):
        a0, a1 = UnitCellPresets().unit_cells[self.unitcell_shape]
        return [self.W * a0, self.L * a1]

class HoneycombClusterPresets:
    clusters = { }
    def add_cluster(self, label, unitcell, W, L):
        #exception will be raised if "label" doesn't match any of the cluster presets
        self.clusters[label] = HoneycombInfo(unitcell, W, L).create_lattice_vectors()

    def __init__(self):
        self.add_cluster(   '24-RE',    'rect', 3, 2)
        self.add_cluster('24-RH120', 'rhom120', 2, 2)
        self.add_cluster( '18-RH60',  'rhom60', 3, 3)
        self.add_cluster(   '12-RE',    'rect', 3, 1)
        self.add_cluster( '12-RH60',  'rhom60', 2, 3)
        self.add_cluster(  '8-RH60',  'rhom60', 1, 4)
        self.add_cluster( '6-RH120', 'rhom120', 1, 1)
        self.add_cluster(  '6-RH60',  'rhom60', 1, 3)
        self.add_cluster(    '4-RE',    'rect', 1, 1)
        self.add_cluster(  '4-RH60',  'rhom60', 1, 2)
        self.add_cluster(  '2-RH60',  'rhom60', 1, 1)

#parameters for model 1
def parameterize_multipole_by_angles(theta, phi, jb):
    return [cos(theta*pi), jb, sin(theta*pi)*cos(phi*pi), sin(theta*pi)*sin(phi*pi)]

#parameters for model 2
def parameterize_multipole_by_epsilon(eps):
    return [1, (1+eps)/np.sqrt(2), 0, (1-eps)/2]

@dataclass
class TwoBodyHamiltonian:
    ga: float
    def make_kitaev_hamiltonian(self, j, k, g, gp):
        __, xy, xz, yz = self.some_helpful_matrices()
        Hz = j*np.eye(3) + k*np.diag([0,0,1]) + g*xy + gp*(xz+yz)
        w = np.array([
                      [0,1,0],
                      [0,0,1],
                      [1,0,0]
                    ])
        Hx = w.T @ Hz @ w
        Hy = w.T @ Hx @ w
        return [(1+self.ga)*Hz, Hx, Hy]#z,x,y hamiltonian

    def make_multipole_hamiltonian(self, jt, jb, jq, jo):
        xxzz, xy, xz, yz = self.some_helpful_matrices()
        H_indep = np.diag([jq+jt/2, jo, jq+jt/2])
        hams = [] #z,x,y hamiltonian
        anisotropy = [self.ga, 0, 0]
        for ga, angle in zip(anisotropy, [n*2*pi/3 for n in [0, 1, 2]]):
            s, c = sin(angle), cos(angle)
            H_dep = (jt/2)*(c*xxzz-s*xz) + jb*(s*xy + c*yz)
            hams.append((1+ga)*(H_indep + H_dep))
        return hams

    def some_helpful_matrices(self):
        xxzz = np.diag([-1, 0, 1])
        xy = np.array([
                       [0,1,0],
                       [1,0,0],
                       [0,0,0]
                     ])
        xz = np.array([
                       [0,0,1],
                       [0,0,0],
                       [1,0,0]
                     ])
        yz = np.array([
                       [0,0,0],
                       [0,0,1],
                       [0,1,0]
                     ])
        return xxzz, xy, xz, yz


@dataclass
class StandardInput:
    model:          str
    method:         str
    lattice:        str

    a0:             np.ndarray
    a1:             np.ndarray

    two_sz:         int

    H0:             np.ndarray
    H1:             np.ndarray
    H2:             np.ndarray

    restart:        str
    lanczos_max:    int
    exct:           int
    lanczos_target: int
    lanczos_eps:    int
    output_mode:    str
    eigenvec_io:    str
    ham_io:         str

    def check_availability(self):
        modelQ  = self.model in ['Spin', 'SpinGC']
        methodQ = self.method in ['Lanczos', 'Full Diag', 'CG']
        latticeQ = self.lattice in ['Honeycomb Lattice']
        restartQ = self.restart in ['None', 'Restart_out', 'Restart_in', 'Restart']
        outputQ = self.output_mode in ['none', 'correlation', 'full']
        eigenvecQ = self.eigenvec_io in ['None', 'Out', 'In']
        hamQ = self.ham_io in ['None', 'Out', 'In']

        bool_array = np.array([modelQ, methodQ, latticeQ, restartQ,
                               outputQ, eigenvecQ, hamQ])
        if bool_array.any() == False:
            print('Your standard input goes beyond what I can do so far!')
            raise SystemExit

    def ensure_consistency(self):
        self.two_sz = self.two_sz if (self.model != 'SpinGC') else None
        self.ham_io = self.ham_io if (self.method == 'Full Diag') else None
        convergingenoughQ = self.lanczos_target > self.exct
        if convergingenoughQ == False:
            print('Inconsistent file!')
            raise SystemExit

    def process_hamiltonians(self):
        hams = [self.H0, self.H1, self.H2]
        ids = np.arange(0,len(hams))
        texts = []
        for ham, id in zip(hams, ids):
            labels = np.array([f'J{id}x  = ', f'J{id}xy  = ', f'J{id}xz  = ',
                               f'J{id}yx = ', f'J{id}y   = ', f'J{id}yz  = ',
                               f'J{id}zx = ', f'J{id}zy  = ', f'J{id}z   = '])
            values = ham.round(10).reshape(9)
            texts.append(self.process(labels, values))
        ''.join(texts)
        return ''.join(texts)

    def process_vectors(self):
        labels = np.array(['a0W = ', 'a0L = ', 'a1W = ', 'a1L = '])
        values = np.concatenate((self.a0, self.a1))
        return self.process(labels, values)

    def process(self, labels, values):
        values = values.astype(str)
        nls = np.array(['\n']*len(labels))
        text = np.char.add(np.char.add(labels,values),nls)
        return ''.join(text)

    def create_file(self, filename):
        f = open(filename, 'w')
        f.write('//------type of calculation\n')
        f.write(f'model = "{self.model}"\n')
        f.write(f'method = "{self.method}"\n')
        f.write(f'lattice = "{self.lattice}"\n')
        f.write('//------lattice parameters\n')
        f.write(self.process_vectors())
        if not self.two_sz == None:
            f.write('//------conserved quantities\n')
            f.write(f'2Sz = {self.two_sz}\n')
        f.write('//------Hamiltonians on each bond\n')
        f.write(self.process_hamiltonians())
        f.write('//------Numerical parameters for simulation\n')
        f.write(f'Restart = "{self.restart}"\n')
        f.write(f'Lanczos_max = {self.lanczos_max}\n')
        f.write(f'exct = {self.exct}\n')
        f.write(f'LanczosTarget = {self.lanczos_target}\n')
        f.write(f'LanczosEps = {self.lanczos_eps}\n')
        f.write(f'OutputMode = "{self.output_mode}"\n')
        f.write(f'EigenVecIO = "{self.eigenvec_io}"\n')
        if not self.ham_io == None:
            f.write(f'HamIO = "{self.ham_io}"\n')
        f.close()

@dataclass
class StandardInputAux:
    loc_file:        str
    trans_file:      str

    hstrength:       float
    hdirection:      np.ndarray

    def check_file_existence(self):
        #check if required files exist, and if so, slurp em
        necessary_files = [loc_file, trans_file]
        for file in necessary_files:
            if not os.path.exists(file):
                print(file + ' does not exist! run a *_prepare.py file and' +
                             ' perform a dry standard run first.')
                raise SystemExit

    def extract_num_sites(self):
        #extract number of sites using locspin file format in HPhi documentation (4.2.4)
        with open(self.loc_file, 'r') as f:
            self.sites = len(f.readlines())-5

    def add_to_trans(self):
        string_list = self.create_field_information()

        #modify header; NTransfer should be equal to num_sites * 2 * 3
        num_terms = len(string_list)
        header_line = 'NTransfer' + ' '*7 + f'{num_terms}\n'

        with open(self.trans_file, 'r') as f:
            lines = f.readlines()
        lines[1] = header_line
        with open(self.trans_file, 'w') as f:
            f.writelines(lines)

        #add the field information to trans file
        string = ''.join(string_list)
        with open(self.trans_file, 'a') as f:
            f.write(string)

    def create_field_information(self):
        hx, hy, hz = self.hstrength * self.hdirection/np.linalg.norm(self.hdirection)

        #[(first term), (second term)]
        #(create_spin_index, annihilate_spin_index, real part, imag part)
        hx_spin_terms = [(0,1,+0.5*hx,+0.0*hx),(1,0,+0.5*hx,+0.0*hx)]
        hy_spin_terms = [(0,1,+0.0*hy,-0.5*hy),(1,0,+0.0*hy,+0.5*hy)]
        hz_spin_terms = [(0,0,+0.5*hz,+0.0*hz),(1,1,-0.5*hz,+0.0*hz)]

        string_list = []
        for spin_terms in [hx_spin_terms, hy_spin_terms, hz_spin_terms]:
            for site in range(self.sites):
                for spin_term in spin_terms:
                    string  = ' '*4 \
                            + f'{site}' + ' '*5 + f'{spin_term[0]}' + ' '*5 \
                            + f'{site}' + ' '*5 + f'{spin_term[1]}' + ' '*9 \
                            + f'{spin_term[2]:.15f}' + ' '*8 \
                            + f'{spin_term[3]:.15f}' \
                            + f'\n'
                    string_list.append(string)
        return string_list
