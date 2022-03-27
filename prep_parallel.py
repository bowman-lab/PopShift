
#########################################################################################################
## This script prepares both ligand and protein files for docking by using ADFR suite scripts.         ## 
## It uses multiprocessing to speed up the prep time.                                                  ##
## It also aligns protein structures if needed, currently using backbone atoms                         ##
#########################################################################################################

#python prep_parallel.py -p /full/path/to/protein/files -l /full/path/to/ligands/files -o /full/path/to/output -a "no" if not aligning, /full/path/to/reference/structure/otherwise

#TODO might be interesting to add the function of adding amber charges instead of only vina charges to the protein
#! BUG: seems to convert all the ligands to a random ligand with different charges in each mol2 file - antechamber

import os
import subprocess as sp
import sys
import mdtraj as md
import argparse
import glob
import multiprocessing as mp

executables = '/project/bowmanlab/bnovak/ADFRsuite_x86_64Linux_1.0/bin' # Path to ADFR suite where prep scripts are

# Adds vina charges to a ligand
def preplig_vina(ligand):
    os.chdir(path_lig)
    sp.run(['%s/prepare_ligand' % executables, '-l', ligand,])
    return print("Ligand has been converted to pdbqt file with vina charges")

# Adds antechamber charges to the ligand pdb
def add_charges(ligand):
    #os.chdir(path_lig)
    #name = (ligand.split('/')[-1]).split('.')[0]
    name = ligand.split('.')[0]
    sp.run('antechamber -i %s -fi pdb -o %s.mol2 -fo mol2 -c bcc -at gaff2' % (ligand, name),shell=True)
    return '%s.mol2' % name

# Converts to pdbqt using user provided charges
def preplig(ligand):
    os.chdir(path_lig)
    sp.run(['%s/prepare_ligand' % executables, '-l', ligand, '-C'])
    return print("Ligand has been converted to pdbqt file using provided charges")

# Converts protein pdb to pdbqt file
def prep_receptor(receptor, out, name):
    sp.run(['%s/prepare_receptor' % executables, '-r', receptor, '-o', '%s/%sqt' % (out, name)])
    return '%s/%sqt' % (out, name)

# Aligns protein files using user provided atoms
def align(protein_file,reference_file,output, atoms):
    pdb = md.load(protein_file)
    atoms = pdb.topology.select(atoms)
    aligned = pdb.superpose(reference_file, atom_indices=atoms)
    prot_name = protein_file.split('/')[-1].split('.')[0]
    aligned.save_pdb('%s/%s_aligned.pdb' % (output,prot_name))
    return '%s/%s_aligned.pdb' % (output, prot_name)

charge_methods={
    'vina': preplig_vina,
    'antechamber': add_charges
}

parser = argparse.ArgumentParser()
parser.add_argument('protein_dir',
                    help='Path to protein directory')
parser.add_argument('output',
                    help='Path to the output')
parser.add_argument('-l','--ligand_dir',
                    help='Path to ligand directory')
parser.add_argument('-a','--align', default='no',
                    choices=['yes','no'],
                    help='Whether to align the structures. Default: no')
parser.add_argument('-r','--reference', default='None',
                    help='Path to reference structure to use for aligning')
parser.add_argument('-x','--atoms', default='backbone',
                    help='What atoms to use for aligning. Default: backbone')
parser.add_argument('-c','--charge', default='vina',
                    choices=['antechamber', 'vina'],
                    help='What charges to add the ligand, Default: vina')

args = parser.parse_args()
path_lig = args.ligand_dir
path_prot = args.protein_dir
path_output = args.output
to_align = args.align
reference = args.reference

protein_name = path_output.split('/')[-1]
n_procs = mp.cpu_count()
pool = mp.Pool(processes=n_procs)

try:
    os.makedirs('%s' % path_output)
except FileExistsError:
    pass

# Checking whether ligands need to be prepared as well, by checking if --ligands_dir flag is set
#TODO change so that either mol2 or a pdb can be input
if args.ligand_dir is not None:
    print('Adding charges to the ligand')
    ligands = sorted(glob.glob('%s/*pdb' % path_lig))  #! Only takes in pdbs for now
    charged_ligands = list(pool.map(charge_methods[args.charge], ligands))
    if args.charge == 'antechamber':
        list(pool.map(preplig,charged_ligands))
else:
    print('Assuming ligand charges are already added')

if to_align == 'yes':
    reference_protein = md.load(reference)
    frames = sorted(glob.glob('%s/*pdb' % path_prot))
    prot_name = [frame.split('/')[-1] for frame in frames]
    output_list = [path_output for l in range(len(frames))]
    reference_protein = [i for i in reference_protein for l in range(len(frames))]
    atoms = [args.atoms for l in range(len(frames))]
    arguments = zip(frames, reference_protein, output_list,atoms)
    aligned = list(pool.starmap(align, arguments))
    arguments = zip(aligned, output_list, prot_name)
    list(pool.starmap(prep_receptor, arguments))
else:
    frames = sorted(glob.glob('%s/*pdb' % path_prot))
    prot_name = [frame.split('/')[-1] for frame in frames]
    output_list = [i for i in [path_output] for l in range(len(frames))]
    arguments = zip(frames, output_list,prot_name)
    list(pool.starmap(prep_receptor, arguments))

pool.close()

