import phonopy
import os
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.interface.calculator import read_crystal_structure
from phonopy.interface.phonopy_yaml import PhonopyYaml
from phonopy.phonon.band_structure import BandStructure, get_band_qpoints_by_seekpath
from phonopy.harmonic.dynamical_matrix import DynamicalMatrix
from phonopy.structure.cells import Primitive
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from phonopy.units import VaspToTHz
import numpy as np
from phonopy.structure.cells import *


def cal_band_by_phonopy_method(data_dir):
    phpy = PhonopyYaml(settings={"force_constants": True})
    phpy.read(os.path.join(data_dir, "phonopy.yaml"))
    qpoints, labels, path_connections = get_band_qpoints_by_seekpath(phpy.primitive, npoints=101)

    dm = DynamicalMatrix(phpy.supercell, Primitive(phpy.supercell, [[1 / 2, 0, 0], [0, 1 / 2, 0], [0, 0, 1 / 2]]),
                         phpy.force_constants)
    band = BandStructure(qpoints, dm, path_connections=path_connections, labels=labels,
                         is_legacy_plot=True)

    # n = len([x for x in path_connections if not x])
    # fig = plt.figure()
    # axs = ImageGrid(
    #     fig,
    #     111,  # similar to subplot(111)
    #     nrows_ncols=(1, n),
    #     axes_pad=0.11,
    #     label_mode="L")

    fig, axs = plt.subplots(1, 1)
    band.plot(axs)
    plt.show()


def check_phono_band(data_dir):
    try:
        phonon = phonopy.load(os.path.join(data_dir, "phonopy_disp.yaml"),
                              force_sets_filename=os.path.join(data_dir, "FORCE_SETS"),
                              force_constants_filename=os.path.join(data_dir, "FORCE_CONSTANTS"))
    except FileNotFoundError:
        return False
    qpoints, labels, path_connections = get_band_qpoints_by_seekpath(phonon.primitive, npoints=101)
    prim_mat = np.linalg.inv(phonon.supercell_matrix)
    dm = DynamicalMatrix(phonon.supercell, Primitive(phonon.supercell, prim_mat), phonon.force_constants)
    band = BandStructure(qpoints, dm, path_connections=path_connections, labels=labels,
                         is_legacy_plot=True)
    is_right_band = True
    for freq_per_path in band.frequencies:
        if freq_per_path.min() < -1e-2:
            is_right_band = False
            break
    fig, axs = plt.subplots(1, 1)
    band.plot(axs)
    plt.title(is_right_band)
    plt.show()
    return is_right_band


def cal_band(data_dir):
    phpy = PhonopyYaml(settings={"force_constants": True})
    phpy.read(os.path.join(data_dir, "phonopy.yaml"))
    qpoints, labels, path_connections = get_band_qpoints_by_seekpath(phpy.primitive, npoints=101)
    qpath = []
    for path in qpoints:
        qpath.extend(path)
    vec_info = []
    supercell_matrix = phpy.supercell_matrix
    # primcell = phpy.primitive
    primitive_matrix = np.linalg.inv(supercell_matrix)
    primcell = Primitive(phpy.supercell, primitive_matrix)
    force_constants = phpy.force_constants
    supercell = phpy.supercell
    relative_distances = []
    for i, position in enumerate(supercell.scaled_positions):
        position = np.dot(supercell_matrix, position)
        p_i = primcell.p2p_map[primcell.s2p_map[i]]
        relative_distances.append({"p_i": p_i,
                                   "distance": list(np.array(np.subtract(position, primcell.scaled_positions[p_i]),
                                                             dtype=int))})
    print(relative_distances)
    for i in range(-supercell_matrix[0][0]+1, supercell_matrix[0][0]):
        for j in range(-supercell_matrix[1][1]+1, supercell_matrix[1][1]):
            for k in range(-supercell_matrix[2][2]+1, supercell_matrix[2][2]):
                for p_i, p_position in enumerate(primcell.scaled_positions):
                    for r_i, r_position in enumerate(primcell.scaled_positions):
                        vec = np.add(np.subtract(r_position, p_position), [i, j, k])
                        a = [it if it > 0 else 0 for it in [i, j, k]]
                        b = [abs(it) if it < 0 else 0 for it in [i, j, k]]
                        a_i = relative_distances.index({"p_i": p_i, "distance": a})
                        b_i = relative_distances.index({"p_i": r_i, "distance": b})
                        fc = force_constants[a_i][b_i]
                        vec_info.append({"vector": vec, "l": [i, j, k], "i-j": [p_i, r_i], "fc": fc, "fc-ab": [a_i, b_i]})
    print(vec_info)
    num_atom = len(primcell)

    frequencies = []
    for qpoint in qpath:
        dm = np.zeros((3 * num_atom, 3 * num_atom), dtype="c%d" % (np.dtype("double").itemsize * 2))
        for info in vec_info:
            i, j = info["i-j"]
            sqrt_mm = np.sqrt(primcell.masses[i] * primcell.masses[j])
            vec = info["vector"]
            fc = info["fc"]
            phase = np.vdot(vec, qpoint) * 2j * np.pi
            phase_factor = np.exp(phase)
            dm[(i * 3): (i * 3 + 3), (j * 3): (j * 3 + 3)] += fc * phase_factor / sqrt_mm
        dm = (dm + dm.conj().transpose()) / 2
        eigvalues = np.linalg.eigvalsh(dm).real
        frequencies.append(np.sqrt(abs(eigvalues)) * np.sign(eigvalues) * VaspToTHz)
    plt.plot(frequencies, marker="o", c="red", linestyle=None, markersize=1, linewidth=0)
    plt.show()


def cal_band_ori(data_dir):
    phpy = PhonopyYaml(settings={"force_constants": True})
    phpy.read(os.path.join(data_dir, "phonopy.yaml"))
    qpoints, labels, path_connections = get_band_qpoints_by_seekpath(phpy.primitive, npoints=101)
    num_atom = len(phpy.primitive)
    frequencies = []

    qpath = []
    for path in qpoints:
        qpath.extend(path)

    for q in qpath:
        dm = np.zeros((3 * num_atom, 3 * num_atom), dtype="c%d" % (np.dtype("double").itemsize * 2))
        primitive = phpy.primitive
        for i in range(len(phpy.supercell)):
            p_i = i // 8
            for j in range(len(phpy.supercell)):
                p_j = j // 8
                sqrt_mm = np.sqrt(phpy.supercell.masses[i] * phpy.supercell.masses[j])
                vec = phpy.supercell.scaled_positions[j]-phpy.supercell.scaled_positions[i]
                vec = np.dot(phpy.supercell_matrix, vec)
                phase = np.vdot(vec, q) * 2j * np.pi
                phase_factor = np.exp(phase)
                dm[(p_i * 3): (p_i * 3 + 3), (p_j * 3): (p_j * 3 + 3)] += phpy.force_constants[i][j] * phase_factor / sqrt_mm
        dm = (dm + dm.conj().transpose()) / 2
        eigvalues = np.linalg.eigvals(dm)
        frequencies.append(np.sqrt(abs(eigvalues)) * np.sign(eigvalues) / VaspToTHz)
    plt.plot(frequencies, marker="o", c="red", linestyle=None, markersize=1, linewidth=0)
    plt.show()


class PhonoParser:

    def __init__(self, phonopy_disp_yaml, force_constants_file, force_sets_file):
        try:
            phonon = phonopy.load(phonopy_yaml=phonopy_disp_yaml,
                                  force_constants_filename=force_constants_file,
                                  force_sets_filename=force_sets_file)
        except FileNotFoundError:
            return
        self.force_constants = phonon.force_constants
        self.supercell = phonon.supercell
        self.supercell_mat = phonon.supercell_matrix
        self.primcell_mat = np.linalg.inv(self.supercell_mat)
        self.primcell = phonon.primitive

        self.s2p_map = []
        self.supercell_relative_lattice = []

        self.dynamical_matrix = None
        self.set_dynamical_matrix_zero()
        self.qpoints = self.get_qpoints()

        self.parse_supercell()
        self.supercell_info = [{"IndexInPrimcell": self.s2p_map[i], "Distance": self.supercell_relative_lattice[i]}
                               for i in range(len(self.supercell))]

    def set_dynamical_matrix_zero(self):
        num_atom = len(self.primcell)
        self.dynamical_matrix = np.zeros((3 * num_atom, 3 * num_atom), dtype="c%d" % (np.dtype("double").itemsize * 2))

    def parse_supercell(self):
        for s_i, s_p in enumerate(self.supercell.scaled_positions):
            s_p = np.dot(self.supercell_mat, s_p)
            for p_i, p_p in enumerate(self.primcell.scaled_positions):
                distance = np.subtract(s_p, p_p)
                j = np.mod(distance, 1) < 1e-5
                if np.sum(j) == len(j):
                    self.s2p_map.append(p_i)
                    self.supercell_relative_lattice.append(list(np.array(distance, dtype=int)))
                    continue

    def get_lattice(self):
        lattice = []
        for i in range(-self.supercell_mat[0][0]+1, self.supercell_mat[0][0]):
            for j in range(-self.supercell_mat[1][1]+1, self.supercell_mat[1][1]):
                for k in range(-self.supercell_mat[2][2]+1, self.supercell_mat[2][2]):
                    lattice.append([i, j, k])
        return lattice

    def get_qpoints(self):
        qpoints = []
        qpaths, labels, path_connections = get_band_qpoints_by_seekpath(self.primcell, npoints=101)
        for qpath in qpaths:
            qpoints.extend(qpath)
        return qpoints

    def get_fc_with_pair_atoms(self, a_i, b_i, r_l):
        a_l = [abs(d) if d < 0 else 0 for d in r_l]
        b_l = [d if d > 0 else 0 for d in r_l]
        a_info = {"IndexInPrimcell": a_i, "Distance": a_l}
        b_info = {"IndexInPrimcell": b_i, "Distance": b_l}
        i = self.supercell_info.index(a_info)
        j = self.supercell_info.index(b_info)
        fc = self.force_constants[i][j]
        return fc

    def cal_dynamical_matrix_with_qpoint(self, qpoint, lang="C"):
        if lang == "C":
            dm = DynamicalMatrix(self.supercell,
                                 Primitive(self.supercell, self.primcell_mat),
                                 self.force_constants)
            dm.run(qpoint, lang)
            self.dynamical_matrix = dm.dynamical_matrix
        else:
            lattice = self.get_lattice()

            for r_l in lattice:
                for a_i, a_p in enumerate(self.primcell.scaled_positions):
                    for b_i, b_p in enumerate(self.primcell.scaled_positions):
                        vec = np.subtract(np.add(b_p, r_l), a_p)
                        fc = self.get_fc_with_pair_atoms(a_i, b_i, r_l)
                        phase = np.vdot(vec, qpoint) * 2j * np.pi
                        phase_factor = np.exp(phase)
                        sqrt_mm = np.sqrt(self.primcell.masses[a_i] * self.primcell.masses[b_i])
                        self.dynamical_matrix[a_i*3:a_i*3+3, b_i*3:b_i*3+3] += fc * phase_factor / sqrt_mm

    def cal_band(self, qpoints=None):
        if qpoints is None:
            qpoints = self.qpoints
        band = []
        for q in qpoints:
            self.cal_dynamical_matrix_with_qpoint(q, lang="l")
            eigval = np.linalg.eigvalsh(self.dynamical_matrix).real
            freq = np.sqrt(abs(eigval)) * np.sign(eigval) * VaspToTHz
            band.append(freq)
            self.set_dynamical_matrix_zero()
        return band

    def fourier_transform_of_dynamical_matrix(self, q_mesh=None):
        if q_mesh is None:
            q_mesh = [8, 8, 8]
        dm = self.dynamical_matrix




if __name__ == "__main__":
    # data_dir = "C:/Users/11159/Desktop/TMP/GdPt2B/nosoc/phono/"
    data_dir = "C:/Users/11159/Desktop/TMP/Si/phono/"
    phonopy_disp = os.path.join(data_dir, "phonopy_disp.yaml")
    force_sets = os.path.join(data_dir, "FORCE_SETS")
    force_constants = os.path.join(data_dir, "FORCE_CONSTANTS")
    phono_parser = PhonoParser(phonopy_disp, force_constants, force_sets)
    print(phono_parser.s2p_map)
    print(phono_parser.supercell_relative_lattice)
    print(phono_parser.supercell_info)
    print(phono_parser.supercell_info.index({'IndexInPrimcell': 0, 'Distance': [1, 1, 1]}))
    band = phono_parser.cal_band()
    plt.plot(band, marker="o", c="red", linestyle=None, markersize=1, linewidth=0)
    plt.show()
    cal_band_by_phonopy_method(data_dir)
    # cal_band(data_dir)
    # cal_band_ori(data_dir)
    # home_path = "C:/Users/11159/Desktop/TMP/HP/test"
    # for dir in os.listdir(home_path):
    #     if os.path.isdir(os.path.join(home_path, dir)):
    #         print(check_phono_band(os.path.join(home_path, dir)))