import pyprocar
import matplotlib.pyplot as plt
import matplotlib as mpl
import copy
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import BSDOSPlotter, \
    BSPlotter, BSPlotterProjected, DosPlotter
import os
import pandas as pd
import numpy as np

# mpl.use('Agg')


def plot_band_dos(data_dir):
    # read vasprun.xml，get band and dos information
    bs_vasprun = Vasprun(data_dir + "band/" + "vasprun.xml", parse_projected_eigen=True)
    bs_data = bs_vasprun.get_band_structure(line_mode=True, kpoints_filename=data_dir + "band/" + "KPOINTS")

    dos_vasprun = Vasprun(data_dir + "scf/" + "vasprun.xml")
    dos_data = dos_vasprun.complete_dos

    # set figure parameters, draw figure
    banddos_fig = BSDOSPlotter(bs_projection='elements', dos_projection='elements',
                               vb_energy_range=10, cb_energy_range=10, fig_size=(32, 18))
    banddos_fig.get_plot(bs=bs_data, dos=dos_data)
    plt.savefig(data_dir + "profig/" + "banddos_fig.png")


def plot_fatbands(data_dir, element_dic, orbital_dic):
    if not os.path.exists(data_dir + "profig/fatband/"):
        os.mkdir(data_dir + "profig/fatband/")

    # orbital_dic = {"s": [0], "p": [1, 2, 3], "d": [4, 5, 6, 7, 8]}
    for element, ele_idx in element_dic.items():
        for orbital, orb_idx in orbital_dic.items():
            pyprocar.bandsplot(data_dir + "band/" + 'PROCAR',
                               outcar=data_dir + "band/" + "OUTCAR",
                               elimit=[-10, 10],
                               cmap='jet',
                               mode='parametric',
                               atoms=copy.deepcopy(ele_idx),
                               orbitals=copy.deepcopy(orb_idx),
                               kpointsfile=data_dir + "band/" + "KPOINTS",
                               savefig=data_dir + "profig/fatband/" + element + "_" + orbital)


def sumo_plot(data_dir):
    if not os.path.exists(os.path.join(data_dir, "profig/sumo")):
        os.mkdir(os.path.join(data_dir, "profig/sumo"))
    # from sumo.cli.kgen import kgen
    # kgen(filename=data_dir + "band/" + "POSCAR",
    #      directory=data_dir + "profig/" + "sumo/")
    from sumo.cli.bandplot import bandplot
    bandplot(filenames=data_dir + "band/vasprun.xml",
             directory=data_dir + "profig/sumo/",
             dos_file=data_dir + "scf/vasprun.xml",
             height=18, width=32, ymin=-10, ymax=10)

    from sumo.electronic_structure import dos

    # bandplot(filenames=data_dir + "band/vasprun.xml",
    #          directory=data_dir + "profig/sumo/fatband/",
    #          dos_file=data_dir + "scf/vasprun.xml",
    #          projection_selection=[("Co", "p"), ("Co", "d"), ("S", "p"), ("Sn", "p")],
    #          mode="stacked",
    #          height=18, width=32, ymin=-5, ymax=5)

    from sumo.cli.dosplot import dosplot
    dosplot(filename=data_dir + "scf/vasprun.xml",
            directory=data_dir + "profig/sumo/",
            height=18, width=32, xmin=-10, xmax=10)


def plot_wannier(data_dir):
    with open(data_dir + "wannier/wannier90_band.dat", "r") as f:
        lines = f.readlines()
    wannier_bs_data = []
    for line in lines:
        if len(line.split()) == 2:
            wannier_bs_data.append([float(d) for d in line.split()])
    wannier_bs_data = np.array(wannier_bs_data).transpose()
    print(wannier_bs_data)
    plt.plot(wannier_bs_data[0], wannier_bs_data[1], marker="o", c="red", linestyle=None, markersize=1, linewidth=0)
    plt.show()


def analyse(data_dir, element_dic):
    # data_dir = "C:/Users/11159/Desktop/TMP/GdPt2B/nosoc/"
    # data_dir = "C:/Users/11159/Desktop/TMP/Cu2Hg1Se4Sn1/nosoc/"
    if not os.path.exists(os.path.join(data_dir, "profig/")):
        os.mkdir(os.path.join(data_dir, "profig"))
    plot_band_dos(data_dir)
    sumo_plot(data_dir)
    # element_dic = {"Cu": [0, 1], "Hg": [2], "Se": [3, 4, 5, 6], "Sn": [7]}
    # element_dic = {"Co": [0, 1, 2], "Sn": [3, 4], "S": [5, 6]}
    # element_dic = {"B": [0, 1, 2], "Gd": [3, 4, 5], "Pt": [6, 7, 8, 9, 10, 11]}
    orbital_dic = {"s": [0], "p": [1, 2, 3], "d": [4, 5, 6, 7, 8]}
    plot_fatbands(data_dir, element_dic, orbital_dic)


if __name__ == "__main__":
    analyse(data_dir="C:/Users/11159/Desktop/TMP/Co3Sn2S2/soc/",
            element_dic={"Co": [0, 1, 2], "Sn": [3, 4], "S": [5, 6]})
    # plot_wannier(data_dir="C:/Users/11159/Desktop/TMP/Cu2Hg1Se4Sn1/nosoc/")