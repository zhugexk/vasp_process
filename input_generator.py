from pymatgen.io.vasp import Poscar, Kpoints
from pymatgen.symmetry.kpath import KPathLatimerMunro, KPathSeek
from pyxtal import pyxtal
from pyxtal.lattice import Lattice


def generate_kpath(poscar_file, kpath_file):
    poscar = Poscar.from_file(poscar_file)
    structure = poscar.structure
    kpath = KPathSeek(structure)
    print(kpath.kpath)
    print(kpath.structure)
    Kpoints.automatic_linemode(20, kpath).write_file(kpath_file)


def generate_poscar(poscar_file):
    struc = pyxtal()
    l1 = Lattice(ltype="hexagonal", volume=195.49,
                 matrix=[[2.672965, -4.629711, 0.0], [2.672965, 4.629711, 0.0], [0.0, 0.0, 7.89876]])
    l2 = Lattice(ltype="hexagonal", volume=1,
                 matrix=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    struc.build(group=181, species=['Gd', 'Pt', 'B'], numIons=[3, 6, 3], lattice=l1,
                sites=[{"3c": [1/2, 0, 0]}, {"6i": [0.151669, 0.303329, 0]}, {"3d": [1/2, 0, 1/2]}])
    print(struc)
    struc.to_file(poscar_file, fmt="poscar")


generate_poscar("C:/Users/11159/Desktop/TMP/GdPt2B/POSCAR")
