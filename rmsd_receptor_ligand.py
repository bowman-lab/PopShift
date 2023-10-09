"""
Computes the RMSD between the receptor and ligand across a 
sampled ensemble and a reference structure. Aligns sampled
receptor and ligand conformations by reference receptor conformation.
    Copyright (C) 2023 Louis G. Smith

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
import argparse
from pathlib import Path
import json
import loos
from loos import pyloos
from sys import argv
import re
from spyrmsd import rmsd
import numpy as np


def atomic_group_to_adjacency(ag: loos.AtomicGroup):
    n = len(ag)
    adjacency = np.zeros((n, n), dtype=int)
    # maybe there's a more efficient way to do this, but this seemed easier
    for i, at in enumerate(ag):
        # assuming that the atom IDs don't have a gap and start at 1
        bi = np.array(at.getBonds(), dtype=int) - 1
        adjacency[i][bi] = 1
    return adjacency


def get_atomic_nums(ag: loos.AtomicGroup, model_path:Path, informat='pdb'):
    from openbabel import openbabel
    mol = openbabel.OBMol()
    conv = openbabel.OBConversion()
    conv.SetInFormat(informat)
    conv.ReadFile(mol, str(model_path))
    atomic_nums = []
    # note that for openbabel the atom id is (apparently) the zeros based
    # ordering of the group of atoms within a molecule, whereas for LOOS (and the PDB)
    # this would be an atom index, and the 'id' is a field in the data file.
    for lat in ag:
        obat = mol.GetAtomById(lat.index())
        atomic_nums.append(obat.GetAtomicNum())
    return np.array(atomic_nums, dtype=int)


def get_meta_from_path(ligp):
    bin_index = ligp.parent.name
    sample = ligp.name
    return bin_index, sample


# right now this only works for SMINA
def get_ligand_affinity(ligp, search_re=re.compile(r'^REMARK minimizedAffinity')):
    with ligp.open() as f:
        for line in f:
            if search_re.match(line):
                return float(line.split(maxsplit=4)[2])
    return None


def purify_expt_multiconf(atomic_group, preferred_loc='A'):
    for i in atomic_group:
        alt = i.altLoc()
        # residues that have no alternative confs will have altLoc() == ''
        if alt != '':
            if alt != preferred_loc:
                atomic_group.remove(i)


# this expects that the expt_ag is SMALLER than, but has the same number of
# residues as, the sample_ag. It makes a new AtomicGroup that, per residue,
# is composed of pAtoms from sample_ag that have the same names as those in
# expt_ag. This can fix some problems with comparing to experimental structures,
# but its behaviour should be checked. Adult supervision required.
def reductively_unify_ags(sample_ag, expt_ag):
    sample_split = sample_ag.splitByResidue()
    expt_split = expt_ag.splitByResidue()
    retgp = loos.AtomicGroup()
    if len(sample_split) != len(expt_split):
        print('Groups do not have the same number of residues.')
        raise IndexError
    for sample_res, expt_res in zip(sample_split, expt_split):
        for at in expt_res:
            retgp += loos.selectAtoms(sample_res,
                                      'name == "{}" '.format(at.name()))
    return retgp


if __name__ == '__main__':
    p = argparse.ArgumentParser("Read ligand and a receptor pdbs (as from docking);"
                                " align and compute RMSD to reference.",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--all-subset', type=str, default=" && ! hydrogen",
                   help='Subset to all selection strings. To disable, provide empty'
                   ' string as argument.')
    p.add_argument('--alt-loc-pref', type=str, default='A',
                   help='For experimental structures with alternative conformations'
                   ' modeled in, remove all but the provided conformer letter.'
                        ' Usually, if there are multiple confs, they are lettered '
                        'A and B.')
    p.add_argument('--debug', action=argparse.BooleanOptionalAction,
                   help='Throw to get debugging messages printed (to stdout).')
    p.add_argument('--extract-score', action=argparse.BooleanOptionalAction,
                   help='Throw to add scores to emitted json.')
    p.add_argument('--ligand-sel', type=str, default='all',
                   help='Selection to apply to ligand PDBs.')
    p.add_argument('--translate-sel', type=str, default=None,
                   help='Use this selection to subset which atoms are translated '
                   'if files are saved. The selection given to "--all-subset" '
                        'is not appended here.')
    p.add_argument('--write-complex-prefix', '-w', type=str, default=None,
                   help='If provided, write complex coordinates to separate files;'
                   ' paths will match ligand path, but with prefix specified.')
    p.add_argument('--ref-receptor-sel', '-p', default=None, type=str,
                   help='If none provided, use main selection string.')
    p.add_argument('--samples', '-S', type=Path, nargs='+', required=True,
                   help='Receptor conformations in order matched to ligand poses.')
    p.add_argument('--pocket-sel', type=str, required=True,
                   help='Selection to align the reference complex onto samples. '
                   'Should operate on samples, but produce atoms matching '
                        'those of receptor-sel.')
    p.add_argument('--poses', '-P', type=Path, nargs='+', required=True,
                   help='Ligand poses in order matched to receptor conformations.')
    p.add_argument('--spyrmsd', action=argparse.BooleanOptionalAction, default=True,
                   help='If thrown, compute spyrmsd using software from Meli and Biggin, DOI:10.1186/s13321-020-00455-2')
    p.add_argument('sample_model', type=Path,
                   help='A model file that can be used to read the samples.')
    p.add_argument('ref_complex', type=Path,
                   help='Reference system that will be aligned onto receptor '
                   'structures and from which the rmsd is computed.')
    p.add_argument('ref_ligand_sel', type=str,
                   help='Use this selection to find the ligand in the receptor '
                   'ligand complex.')
    p.add_argument('rmsd_db', type=Path,
                   help='Path to write the list of RMSDs, and optionally paths to '
                   'aligned complexes, to. Will be formatted as JSON.')
    try:
        args = p.parse_args()
    except:
        for index, argval in enumerate(argv):
            print(index, argval)
        raise

    subset_str = args.all_subset
    pocket_sel = args.pocket_sel + subset_str
    ligand_sel = args.ligand_sel + subset_str
    if args.ref_receptor_sel:
        ref_receptor_sel = args.ref_receptor_sel + subset_str
    ref_ligand_sel = args.ref_ligand_sel + subset_str
    
    # make selections before trajectory loops.
    ligpose_full = loos.createSystem(str(args.poses[0]))
    ligpose = loos.selectAtoms(ligpose_full, ligand_sel)
    if not ligpose.hasBonds():
        ligpose.findBonds(1.85)
    pose_adjacency = atomic_group_to_adjacency(ligpose)
    ligpose_anums = get_atomic_nums(ligpose, args.poses[0])
    samplec_full = loos.createSystem(str(args.samples[0]))
    samplec = loos.selectAtoms(samplec_full, pocket_sel)
    samplec_ca = loos.selectAtoms(samplec, 'name == "CA"')
    # reference sels
    ref_complex_full = loos.createSystem(str(args.ref_complex))
    purify_expt_multiconf(ref_complex_full, preferred_loc=args.alt_loc_pref)
    if args.ref_receptor_sel:
        ref_complex = loos.selectAtoms(ref_complex_full, ref_receptor_sel)
    else:
        ref_complex = loos.selectAtoms(ref_complex_full, pocket_sel)

    samplec = reductively_unify_ags(samplec, ref_complex)
    ref_complex_ca = loos.selectAtoms(ref_complex, 'name == "CA"')
    ref_ligand = loos.selectAtoms(ref_complex_full, ref_ligand_sel)
    if not ref_ligand.hasBonds():
        ref_ligand.findBonds(1.85)
    ref_ligand_path = args.ref_complex.parent / ('ligand-' + args.ref_complex.name)
    ref_ligand_path.write_text(str(loos.PDB.fromAtomicGroup(ref_ligand)))
    ref_anums = get_atomic_nums(ref_ligand, ref_ligand_path)
    ref_adjacency = atomic_group_to_adjacency(ref_ligand)

    outlist = [{} for i in range(len(set(x.parent.name for x in args.poses)))]
    lig_traj = pyloos.VirtualTrajectory(*map(
        lambda fn: pyloos.Trajectory(str(fn), ligpose_full), args.poses))
    sample_traj = pyloos.VirtualTrajectory(*map(
        lambda fn: pyloos.Trajectory(str(fn), samplec_full), args.samples))

    if args.debug:
        print('ref_complex', args.ref_complex)
        print('ligname', ref_ligand[0].resname())
        if len(ref_complex) != len(samplec):
            print('PROBLEM: ref_complex size:', len(
                ref_complex), 'samplec size:', len(samplec))
        if len(ref_complex_ca) != len(samplec_ca):
            print('PROBLEM: ref_complex_ca size:', len(
                ref_complex_ca), 'samplec_ca size:', len(samplec_ca))
        if len(ligpose) != len(ref_ligand):
            print('PROBLEM: ligpose size:', len(ligpose),
                  'ref_ligand size:', len(ref_ligand))
        print('ref_ligand size', len(ref_ligand))
        print('pose size', len(ligpose))
    else:
        for posep, samplep, _, _ in zip(args.poses, args.samples, lig_traj, sample_traj):
            # superposition computes the transform needed to superimpose atoms from
            # ref_complex onto atoms from samplec, but does not perform transform
            xform = loos.XForm(ref_complex.superposition(samplec))
            # applies to ref_ligand and ref_complex because they stem from this AG
            ref_complex_full.applyTransform(xform)
            rmsd_rec = samplec.rmsd(ref_complex)
            rmsd_ca = samplec_ca.rmsd(ref_complex_ca)
            rmsd_lig = ligpose.rmsd(ref_ligand)
            if args.spyrmsd:
                coords_ref = ref_ligand.getCoords()
                coords = ligpose.getCoords()
                lig_symm_rmsd = rmsd.symmrmsd(
                    coords_ref,
                    coords,
                    ref_anums,
                    ligpose_anums,
                    ref_adjacency,
                    pose_adjacency,
                    minimize=False
                )
            else:
                lig_symm_rmsd = None
            bin_index, sample = get_meta_from_path(posep)
            if args.extract_score:
                score = get_ligand_affinity(posep)
            else:
                score = None
            outlist[int(bin_index)][sample] = {
                'receptor': rmsd_rec,
                'receptor CA': rmsd_ca,
                'ligand': rmsd_lig,
                'score': score,
                'ligand-symmrmsd': lig_symm_rmsd
            }
            if args.write_complex_prefix:
                out_ag = ligpose_full + samplec_full + ref_complex_full
                if args.translate_sel:
                    out_ag = loos.selectAtoms(out_ag, args.translate_sel)
                pdb = loos.PDB.fromAtomicGroup(out_ag)
                outpath = posep.parent + args.write_complex_prefix + posep.name
                with outpath.open('w') as f:
                    f.write(str(pdb))
                outlist[int(bin_index)][sample]['aligned'] = str(outpath)

        with args.rmsd_db.open('w') as f:
            json.dump(outlist, f, indent=4)
