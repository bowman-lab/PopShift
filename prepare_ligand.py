"""
#########################################################################################################
##          This script prepares only ligand files for docking by using ADFR suite scripts.            ## 
#########################################################################################################
    Copyright (C) 2023 Borna Novak

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""
import os
import subprocess as sp
from optparse import OptionParser
import glob

executables = '/project/bowmanlab/bnovak/ADFRsuite_x86_64Linux_1.0/bin' # Path to ADFR suite where prep scripts are

def preplig(ligand):
    sp.run(['%s/prepare_ligand' % executables, '-l', ligand]) # This uses vina charges
    return print("Ligand has been converted to pdbqt file with vina charges")

parser = OptionParser()
parser.add_option('-l','--ligand_dir', dest='ligand_dir', help='Path to ligand directory')
parser.add_option('-o','--out', dest='output', help='Path to the output')

(options, args) = parser.parse_args()
path_lig = options.ligand_dir
path_output = options.output

os.chdir(path_lig)

for ligand in sorted(glob.glob('%s/*mol2' % path_lig)):  #!Having mol2 here should only be a temporary solution, make sure to change this and also take pdb
    print(ligand)
    preplig(ligand)
