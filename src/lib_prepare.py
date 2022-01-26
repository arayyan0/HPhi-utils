from dataclasses import dataclass
from math import sin, cos
import numpy as np

pi = np.pi

@dataclass
class HoneycombInfo:
    edgetype: str
    W:        int
    L:        int

    def create_lattice_vectors(self):
        a0 = np.array([1, 0])
        a1 = np.array([0, 1])
        if self.edgetype == 'armchair':
            a0_t = 2*a0 - a1
            a1_t = 2*a1 - a0
            a0, a1 = a0_t, a1_t
        return [self.W*a0, self.L*a1]

class HoneycombPresets:
    clusters = { }

    def add_cluster(self, label, edgetype, W, L):
        self.clusters[label] = HoneycombInfo(edgetype, W, L).create_lattice_vectors()

    def __init__(self):
        self.add_cluster('24-site', 'armchair', 2, 2)
        self.add_cluster('18-site',   'zigzag', 3, 3)
        self.add_cluster('12-site',   'zigzag', 2, 3)
        self.add_cluster( '6-site', 'armchair', 1, 1)
        self.add_cluster( '4-site',   'zigzag', 1, 2)
        self.add_cluster( '2-site',   'zigzag', 1, 1)
	
class TwoBodyHamiltonian:
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
        return [Hz, Hx, Hy]

    def make_multipole_hamiltonian(self, jt, jb, jq, jo):
        xxzz, xy, xz, yz = self.some_helpful_matrices()
        H_indep = np.diag([jq+jt/2, jo, jq+jt/2])
        hams = []
        for angle in [n*2*pi/3 for n in [0, 1, 2]]:
            s, c = sin(angle), cos(angle)
            H_dep = (jt/2)*(c*xxzz-s*xz) + jb*(s*xy + c*yz)
            hams.append(H_indep + H_dep)
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
        f.write(f'OutputMode = "{self.output_mode}"\n')
        f.write(f'EigenVecIO = "{self.eigenvec_io}"\n')
        if not self.ham_io == None:
            f.write(f'HamIO = "{self.ham_io}"\n')
        f.close()
