import argparse
from pathlib import Path
import json
import loos
from loos import pyloos
from loos.pyloos import options

def get_meta_from_path(ligp):
    bin_index = ligp.parent.name
    sample = ligp.name
    return bin_index, sample


p = options.LoosOptions("Read a ligand and a receptor pdb (as from docking); "
                        "align and compute RMSD to reference.")
p.modelSelectionOptions()
p.add_argument('--all-subset', type=str, default=" && ! hydrogen",
               help='Subset to all selection strings. To disable, provide empty'
                    ' string as argument.')
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
p.add_argument('--poses', '-P', type=Path, nargs='+', required=True,
                help='Ligand poses in order matched to receptor conformations.')
p.add_argument('ref_complex', type=Path,
               help='Reference system that will be aligned onto receptor '
                    'structures and from which the rmsd is computed.')
p.add_argument('ref_ligand_sel', type=str,
               help='Use this selection to find the ligand in the receptor '
                    'ligand complex.')
p.add_argument('rmsd_db', type=Path,
               help='Path to write the list of RMSDs, and optionally paths to '
                    'aligned complexes, to. Will be formatted as JSON.')
args = p.parse_args()

subset_str = args.all_subset
pocket_sel = args.selection + subset_str
ligand_sel = args.ligand_sel + subset_str
if args.ref_receptor_sel:
    ref_receptor_sel = args.ref_receptor_sel + subset_str
ref_ligand_sel = args.ref_liigand_sel + subset_str

# simulated ligand and receptor sels
system = loos.createSystem(args.model)
pocket = loos.selectAtoms(system, pocket_sel)

# reference sels
ref_complex_full = loos.createSystem(args.ref_complex)
if args.ref_receptor_sel:
    ref_complex = loos.selectAtoms(ref_complex_full, ref_receptor_sel)
else:
    ref_complex = loos.selectAtoms(ref_complex_full, pocket_sel)

ref_ligand = loos.selectAtoms(args.ref_ligand_sel,  ref_complex_full)
# make selections before trajectory loops.
ligpose_full = loos.createSystem(str(args.poses[0]))
ligpose = loos.selectAtoms(ligpose_full, ligand_sel)
samplec_full = loos.createSystem(str(args.samples[0]))
samplec = loos.selectAtoms(samplec_full, pocket_sel)
outlist = [{} for i in range(len(set(x.parent.name for x in args.poses)))]
lig_traj = pyloos.VirtualTrajectory(*map(
    lambda fn: pyloos.Trajectory(str(fn), ligpose_full), args.poses))
sample_traj = pyloos.VirtualTrajectory(*map(
    lambda fn: pyloos.Trajectory(str(fn), samplec_full), args.samples))

for posep, samplep, _, _ in zip(args.poses, args.samples, lig_traj, sample_traj):
    xform = ref_complex.alignOnto(samplec)
    # applies to ref_ligand and ref_complex because they stem from this AG
    ref_complex_full.applyTransform(xform)
    rmsd_rec = samplec.rmsd(ref_complex)
    rmsd_lig = ligpose.rmsd(ref_ligand)
    bin_index, sample = get_meta_from_path(posep)
    outlist[int(bin_index)][sample] = {
        'receptor': rmsd_rec,
        'ligand': rmsd_lig,
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
    json.dump(outlist, f)
