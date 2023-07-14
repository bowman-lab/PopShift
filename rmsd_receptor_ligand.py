import argparse
from pathlib import Path
import json
import loos
from loos import pyloos
from sys import argv
import re


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
            bin_index, sample = get_meta_from_path(posep)
            if args.extract_score:
                score = get_ligand_affinity(posep)
            else:
                score = None
            outlist[int(bin_index)][sample] = {
                'receptor': rmsd_rec,
                'receptor CA': rmsd_ca,
                'ligand': rmsd_lig,
                'score': score
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
