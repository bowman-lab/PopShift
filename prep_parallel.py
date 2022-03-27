
#########################################################################################################
## This script prepares both ligand and protein files for docking by using ADFR suite scripts.         ## 
## It uses multiprocessing to speed up the prep time.                                                  ##
## It also aligns protein structures if needed, currently using backbone atoms                         ##
#########################################################################################################

#python prep_parallel.py -p /full/path/to/protein/files -l /full/path/to/ligands/files -o /full/path/to/output -a "no" if not aligning, /full/path/to/reference/structure/otherwise

#TODO might be interesting to add the function of adding amber charges instead of only vina charges to the protein

import os
import subprocess as sp
import sys
import mdtraj as md
#from optparse import OptionParser
import argparse
import glob
import multiprocessing as mp

executables = '/project/bowmanlab/bnovak/ADFRsuite_x86_64Linux_1.0/bin' # Path to ADFR suite where prep scripts are

def preplig(ligand):
    sp.run(['%s/prepare_ligand' % executables, '-l', ligand, '-C']) # This uses the charges prepared using antechamber
    return print("Ligand has been converted to pdbqt file with vina charges")

def prep_receptor(receptor, out, name):
    sp.run(['%s/prepare_receptor' % executables, '-r', receptor, '-o', '%s/%sqt' % (out, name)])
    return '%s/%sqt' % (out, name)

def align(protein_file,reference_file,output, atoms):
    pdb = md.load(protein_file)
    atoms = pdb.topology.select(atoms)
    print('Aligning using %s' % atoms)
    aligned = pdb.superpose(reference_file, atom_indices=atoms)
    prot_name = protein_file.split('/')[-1].split('.')[0]
    aligned.save_pdb('%s/%s_aligned.pdb' % (output,prot_name))
    return '%s/%s_aligned.pdb' % (output, prot_name)

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

if args.ligand_dir is not None:
    print('ligand_dir is set') #! This is where setting up ligands should go
else:
    print('not set')

#TODO prep_ligand and prep_protein should be combines; prep_ligand can only be done if a certain flag exists --ligand yes
#for ligand in sorted(glob.glob('%s/*mol2' % path_lig)):  #Having pdb here should only be a temporary solution, make sure to change this and also take mol2
#    os.chdir(path_lig)
#    preplig(ligand)
#    os.chdir('..')

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
    pool.close()
else:
    frames = sorted(glob.glob('%s/*pdb' % path_prot))
    prot_name = [frame.split('/')[-1] for frame in frames]
    output_list = [i for i in [path_output] for l in range(len(frames))]
    arguments = zip(frames, output_list,prot_name)
    list(pool.starmap(prep_receptor, arguments))
    pool.close()

