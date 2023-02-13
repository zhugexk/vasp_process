import phonopy
from phonopy.harmonic.dynamical_matrix import DynamicalMatrix
from phonopy.phonon.band_structure import BandStructure, get_band_qpoints_by_seekpath
from phonopy.structure.cells import Primitive
import os
from matplotlib import pyplot as plt
import numpy as np
from phonopy.interface.phonopy_yaml import PhonopyYaml
import shutil


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
    plt.savefig(os.path.join(data_dir, "band.jpg"))
    return is_right_band


def check_phono_band_directly(data_dir):
    try:
        phonon = phonopy.load(os.path.join(data_dir, "phonopy_disp.yaml"),
                              force_sets_filename=os.path.join(data_dir, "FORCE_SETS"),
                              force_constants_filename=os.path.join(data_dir, "FORCE_CONSTANTS"))
    except FileNotFoundError:
        return False
    band = PhonopyYaml()
    band.read(os.path.join(data_dir, "band.yaml"))
    print(band)


if __name__ == "__main__":
    home_path = "C:/Users/11159/Desktop/TMP/HP/6001-8000/res"
    useful_res_path = "C:/Users/11159/Desktop/TMP/HP/6001-8000/useful_res"
    if not os.path.exists(useful_res_path):
        os.mkdir(useful_res_path)
    print(len(os.listdir(home_path)))
    for _dir in os.listdir(home_path):
        if os.path.isdir(os.path.join(home_path, _dir)):
            if check_phono_band(os.path.join(home_path, _dir)):
                shutil.copytree(os.path.join(home_path, _dir), os.path.join(useful_res_path, _dir))

