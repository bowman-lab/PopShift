
#########################################################################################################
## This script prepares both ligand and protein files for docking by using ADFR suite scripts.         ## 
## It uses multiprocessing to speed up the prep time.                                                  ##
## It also aligns protein structures if needed, currently using backbone atoms                         ##
#########################################################################################################

#python prep_parallel.py -p /full/path/to/protein/files -l /full/path/to/ligands/files -o /full/path/to/output -a "no" if not aligning, /full/path/to/reference/structure/otherwise

import numpy as np
import os
import subprocess as sp
import sys
import mdtraj as md
from optparse import OptionParser
import glob
import multiprocessing as mp

executables = '/project/bowmanlab/bnovak/ADFRsuite_x86_64Linux_1.0/bin' # Path to ADFR suite where prep scripts are

def preplig(ligand):
    sp.run(['%s/prepare_ligand' % executables, '-l', ligand, '-C']) # This uses the charges prepared using antechamber
    return print("Ligand has been converted to pdbqt file with vina charges")

def prep_receptor(receptor, out, name):
    sp.run(['%s/prepare_receptor' % executables, '-r', receptor, '-o', '%s/%sqt' % (out, name)])
    return '%s/%sqt' % (out, name)

def align(protein_file,reference_file,output):
    pdb = md.load(protein_file)
    #atoms = pdb.topology.select(' or '.join(f'(backbone and resSeq {r})' for r in pocket_residues['myh2-5n6a']))
    atoms = pdb.topology.select('backbone')
    aligned = pdb.superpose(reference_file, atom_indices=atoms)
    prot_name = protein_file.split('/')[-1].split('.')[0]
    aligned.save_pdb('%s/%s_aligned.pdb' % (output,prot_name))
    return '%s/%s_aligned.pdb' % (output, prot_name)

parser = OptionParser()
parser.add_option('-l','--ligand_dir', dest='ligand_dir', help='Path to ligand directory')
parser.add_option('-p','--protein_dir', dest='protein_dir', help='Path to protein directory')
parser.add_option('-a','--align', dest='aligning', help='Do you want to align?', default='no')
parser.add_option('-r','--reference', dest='reference', help='Reference structure to align to', default='None')
parser.add_option('-o','--out', dest='output', help='Path to the output')

(options, args) = parser.parse_args()
path_lig = options.ligand_dir
path_prot = options.protein_dir
path_output = options.output
to_align = options.aligning
reference = options.reference

protein_name = path_output.split('/')[-1]
n_procs = mp.cpu_count()
pool = mp.Pool(processes=n_procs)

try:
    os.makedirs('%s' % path_output)
except FileExistsError:
    pass

#for ligand in sorted(glob.glob('%s/*mol2' % path_lig)):  #Having pdb here should only be a temporary solution, make sure to change this and also take mol2
#    os.chdir(path_lig)
#    preplig(ligand)
#    os.chdir('..')

if to_align == 'yes':
    reference_protein = md.load(reference)
    frames = sorted(glob.glob('%s/*pdb' % path_prot))
    prot_name = [frame.split('/')[-1] for frame in frames]
    output_list = [i for i in [path_output] for l in range(len(frames))]
    print(output_list)
    reference_protein = [i for i in reference_protein for l in range(len(frames))]
    arguments = zip(frames, reference_protein, output_list)
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

