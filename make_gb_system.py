import openmm as mm
from openmm import app as app
import argparse
from pathlib import Path

p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
p.add_argument('receptor_pdb', type=Path,
               help='path to receptor pdb to build openmm system with.')
p.add_argument('ligand_mol2', type=Path,
               help='path to ligand mol2 file to build openmm system with.')
p.add_argument('output_system', type=Path,
               help='Filename to output xml system to.')

