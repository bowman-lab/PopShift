
#########################################################################################################
## This script prepares both ligand and protein files for docking by using ADFR suite scripts.         ## 
## It uses multiprocessing to speed up the prep time.                                                  ##
## It also aligns protein structures if needed, using user provided atoms in -a                        ##
#########################################################################################################

#* Fixed BUG by not using map: seems to convert all the ligands to a random ligand with different charges in each mol2 file - antechamber

import os
import subprocess as sp
import sys
import argparse
import glob
import multiprocessing as mp
from pathlib import Path
from functools import partial

home = os.getcwd()

#executables = '/project/bowmanlab/bnovak/ADFRsuite_x86_64Linux_1.0/bin' # Path to ADFR suite where prep scripts are #*Put scripts in your path instead

# Adds vina charges to a ligand
def preplig_vina(ligand, a_flag=None):
    if a_flag:
        sp.run(['prepare_ligand', '-l', ligand, '-A', a_flag])
    else:
        sp.run(['prepare_ligand', '-l', ligand])
    return print("Ligand has been converted to pdbqt file with vina charges")

# Adds antechamber charges to the ligand pdb
def add_charges(ligand):
    name = ligand.split('.')[0]
    sp.run(f'antechamber -i {ligand} -fi pdb -o {name}.mol2 -fo mol2 -c bcc -at gaff2',shell=True)
    return f'{name}.mol2'

# Converts to pdbqt using user provided charges
def preplig(ligand, a_flag=None):
    os.chdir(path_lig)
    if a_flag:
        sp.run(['prepare_ligand', '-l', ligand, '-C', '-A', a_flag])
    else:
        sp.run(['prepare_ligand', '-l', ligand, '-C'])
    return print("Ligand has been converted to pdbqt file using provided charges")

# Converts protein pdb to pdbqt file
def prep_receptor(receptor, out, name, a_flag=None):
    if a_flag:
        sp.run(['prepare_receptor', '-r', receptor, '-o', f'{out}/{name}.pdbqt', '-A', a_flag])
    else:
        sp.run(['prepare_receptor', '-r', receptor, '-o', f'{out}/{name}.pdbqt'])
    return f'{out}/{name}.pdbqt'

charge_methods={
    'vina': preplig_vina,
    'antechamber': add_charges
}

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p', '--protein_dir',
                    help='Path to protein directory')
parser.add_argument('-l', '--ligand_dir',
                    help='Path to ligand directory')
parser.add_argument('-c', '--charge', default='vina',
                    choices=['antechamber', 'vina'],
                    help='What charges to add the ligand.')
parser.add_argument('-A', '--a-flag', default=None, type=str,
                    help='If provided, pass this flag (types of repairs to make) to calls to prep receptor/ligand. '
                         'Same options as the "A" flag for prepare_receptor and prepare_ligand. Bond-finding in this '
                         'way is extremely slow for larger receptors.')


args = parser.parse_args()
if len(sys.argv) < 2:
    parser.print_help()
    exit(0)
path_lig = args.ligand_dir
path_prot = args.protein_dir

n_procs = mp.cpu_count()
pool = mp.Pool(processes=n_procs)

# Checking whether ligands need to be prepared as well, by checking if --ligands_dir flag is set
if args.ligand_dir is not None:
    os.chdir(path_lig) #* There is probably a more elegant way than changing directories
    print('Adding charges to the ligand')
    ligands = sorted(glob.glob('*pdb')) #TODO change so that either mol2 or a pdb can be input; might not be necessary
    for ligand in ligands:
        charge_methods[args.charge](ligand, a_flag=args.a_flag)
    os.chdir(home)
    #charged_ligands = list(pool.map(charge_methods[args.charge], ligands)) #! Currently not using map, because it makes all the ligands have the same atoms with different charges for yet unknown reason
    if args.charge == 'antechamber':
        charged_ligands = sorted(glob.glob(f'{path_lig}/*mol2'))
        list(pool.map(preplig,charged_ligands))
else:
    print('Assuming ligands have been already converted to pdbqts with partial charges')

# Converting protein pdbs to pdbqts
if args.protein_dir is not None: #Checking whether proteins should be converted too
    frames = sorted(glob.glob(f'{path_prot}/receptor/*/*pdb'))
    prot_name = [Path(frame).stem for frame in frames]
    directory = [os.path.dirname(frame) for frame in frames]
    prep = partial(prep_receptor, a_flag=args.a_flag)
    arguments = zip(frames, directory, prot_name)
    list(pool.starmap(prep, arguments))

pool.close()

