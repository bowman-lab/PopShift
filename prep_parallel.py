
#########################################################################################################
## This script prepares both ligand and protein files for docking by using ADFR suite scripts.         ## 
## It uses multiprocessing to speed up the prep time.                                                  ##
## It also aligns protein structures if needed, using user provided atoms in -a                        ##
#########################################################################################################


#TODO might be interesting to add the function of adding amber charges instead of only vina charges to the protein
#* Fixed BUG by not using map: seems to convert all the ligands to a random ligand with different charges in each mol2 file - antechamber

import os
import subprocess as sp
import sys
import mdtraj as md
import argparse
import glob
import multiprocessing as mp
from pathlib import Path

home = os.getcwd()

#executables = '/project/bowmanlab/bnovak/ADFRsuite_x86_64Linux_1.0/bin' # Path to ADFR suite where prep scripts are #*Put scripts in your path instead

# Adds vina charges to a ligand
def preplig_vina(ligand):
    sp.run(['prepare_ligand', '-l', ligand,])
    return print("Ligand has been converted to pdbqt file with vina charges")

# Adds antechamber charges to the ligand pdb
def add_charges(ligand):
    name = ligand.split('.')[0]
    sp.run('antechamber -i %s -fi pdb -o %s.mol2 -fo mol2 -c bcc -at gaff2' % (ligand, name),shell=True)
    return '%s.mol2' % name

# Converts to pdbqt using user provided charges
def preplig(ligand):
    os.chdir(path_lig)
    sp.run(['prepare_ligand', '-l', ligand, '-C'])
    return print("Ligand has been converted to pdbqt file using provided charges")

# Converts protein pdb to pdbqt file
def prep_receptor(receptor, out, name):
    sp.run(['prepare_receptor', '-r', receptor, '-o', '%s/%s.pdbqt' % (out, name)])
    return '%s/%s.pdbqt' % (out, name)

# Aligns protein files using user provided atoms
# def align(protein_file,reference_file,output, atoms):
#     pdb = md.load(protein_file)
#     atoms = pdb.topology.select(atoms)
#     aligned = pdb.superpose(reference_file, atom_indices=atoms)
#     prot_name = protein_file.split('/')[-1].split('.')[0]
#     aligned.save_pdb('%s/%s_aligned.pdb' % (output,prot_name))
#     return '%s/%s_aligned.pdb' % (output, prot_name)

charge_methods={
    'vina': preplig_vina,
    'antechamber': add_charges
}

parser = argparse.ArgumentParser()
parser.add_argument('-p','--protein_dir',
                    help='Path to protein directory')
parser.add_argument('-l','--ligand_dir',
                    help='Path to ligand directory')
parser.add_argument('-c','--charge', default='vina',
                    choices=['antechamber', 'vina'],
                    help='What charges to add the ligand, Default: vina')

args = parser.parse_args()
path_lig = args.ligand_dir
path_prot = args.protein_dir

n_procs = mp.cpu_count()
pool = mp.Pool(processes=n_procs)

# try:
#     os.makedirs('%s' % path_output)
# except FileExistsError:
#     pass

# Checking whether ligands need to be prepared as well, by checking if --ligands_dir flag is set
if args.ligand_dir is not None:
    os.chdir(path_lig) #* There is probably a more elegant way than changing directories
    print('Adding charges to the ligand')
    ligands = sorted(glob.glob('*pdb')) #TODO change so that either mol2 or a pdb can be input; might not be necessary
    for ligand in ligands:
        charge_methods[args.charge](ligand)
    os.chdir(home)
    #charged_ligands = list(pool.map(charge_methods[args.charge], ligands)) #! Currently not using map, because it makes all the ligands have the same atoms with different charges for yet unknown reason
    if args.charge == 'antechamber':
        charged_ligands = sorted(glob.glob('%s/*mol2' % path_lig))
        list(pool.map(preplig,charged_ligands))
else:
    print('Assuming ligands have been already converted to pdbqts with partial charges')

# Aligning and converting protein pdbs to pdbqts
if args.protein_dir is not None: #Checking whether proteins should be converted too
    frames = sorted(glob.glob('%s/receptor/*/*pdb' % path_prot))
    prot_name = [Path(frame).stem for frame in frames]
    directory = [os.path.dirname(frame) for frame in frames]
    #output_list = ['%s/%s.pdbqt' % (directory[i],prot_name[i]) for i in range(len(directory))]
    arguments = zip(frames, directory,prot_name)
    list(pool.starmap(prep_receptor, arguments))

pool.close()

