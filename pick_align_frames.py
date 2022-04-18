import argparse
from pathlib import Path
import numpy as np
from enspara import ra
import loos
from loos import pyloos


def get_assign_inds(assignments, nstates):
    index_list_by_bin = [[] for i in range(nstates)]
    # handle the case where there is only 1 trajectory,
    # meaning the ragged array doesn't have singleton elements at the tightest level.
    if assignments.shape[0] == 1:
        for i, t in enumerate(assignments):
            for j, f in enumerate(t):
                bin_index = f
                index_list_by_bin[bin_index].append((i, j))
    else:
        for i, t in enumerate(assignments):
            for j, f in enumerate(t):
                # iterating over RA assign contents produces singleton arrays of dtype=int32
                # that need to be unpacked for use
                bin_index = f
                index_list_by_bin[bin_index].append((i, j))
    return [np.array(bin_indices, dtype=[('traj', np.int32), ('frame', np.int32)]) for bin_indices in index_list_by_bin]


def get_random_per_bin(assignments, n_states, number_desired, replace=False, gen=np.random.default_rng()):
    # assumes RA's traj axis is the outer one.S
    all_state_indices = get_assign_inds(assignments, n_states)
    chosen_inds = np.array([gen.choice(state_ixs, number_desired, axis=0, replace=replace)
                            for state_ixs in all_state_indices])
    return chosen_inds


def get_specified_number_per_bin_random(assignments, nstates, numbers_desired, replace=False, gen=np.random.default_rng()):
    all_state_indices = get_assign_inds(assignments, nstates)
    chosen_inds = (gen.choice(state_ixs, number_desired, axis=0, replace=replace)
                   for state_ixs, number_desired in zip(all_state_indices, numbers_desired))
    return list(chosen_inds)


def rip_conformations(chosen_inds, model, subset_selection, align_selection, traj_paths):
    trajectories = [pyloos.Trajectory(traj_path, model, subset=subset_selection) for traj_path in traj_paths]
    full_inds = []
    subset_vec = loos.AtomicGroupVector()
    align_vec = loos.AtomicGroupVector()

    # read selected frames into memory
    for bin_ix, samples in enumerate(chosen_inds):

        for trj_ix, fra_ix in samples:
            frame = trajectories[trj_ix].readFrame(fra_ix)
            align_vec.push_back(loos.selectAtoms(frame, align_selection).copy())
            subset_vec.push_back(frame.copy())
            full_inds.append((bin_ix, trj_ix, fra_ix))

    return full_inds, subset_vec, align_vec


def align_samples(subset_vec, align_vec):
    # iteratively aligns AGs in group. Alignment result contains transforms to
    # rotate/translate the aligned group into alignment with the target.
    alignment_result = loos.iterativeAlignmentPy(align_vec)
    # apply transforms here. This is done in-place
    loos.applyTransforms(subset_vec, alignment_result.transforms)
    print('Iterative alignment final RMSD', alignment_result.rmsd, 'after',
          alignment_result.iterations, 'iterations')


def write_sampled_frames(subset_list, full_inds, out_path, write_bin_trajs):
    # If out_path is already a path, then this returns it.
    # If out_path is a string, turn it into a pathlib Path object here
    op = Path(out_path)
    prev_ix = -1
    for (bin_ix, trj_ix, fra_ix), frame in zip(full_inds, subset_list):
        bp = op / str(bin_ix)
        # This will operate like mkdir -pf; it will overwrite.
        bp.mkdir(parents=True, exist_ok=True)
        if write_bin_trajs:
            if prev_ix != bin_ix:
                outtraj = loos.DCDWriter(str(bp / 'samples.dcd'))
            outtraj.writeFrame(frame)
            prev_ix = bin_ix
        pdb = loos.PDB.fromAtomicGroup(frame)
        pdb_path = bp / '{}-{}.pdb'.format(trj_ix, fra_ix)
        pdb_path.write_text(str(pdb))  # will close file handle after writing


frame_selectors = {
    'random': get_random_per_bin,
    'specified_totals': get_specified_number_per_bin_random
}

parser = argparse.ArgumentParser()
parser.add_argument('receptor_name', type=str,
                    help="Name to use to use as top-level directory for the docking run tree.")
parser.add_argument('model', type=str,
                    help='A loos-interpretable model file that will permit reading the trajectories in "traj_paths".')
parser.add_argument('assignments', type=str,
                    help='h5 file with assignments associated to the MSM used.')
parser.add_argument(metavar='eq_probs|pickled_msm', type=str, dest='eq_probs',
                    help='.npy file with equilibrium probabilities from MSM, or pickled MSM object.')
parser.add_argument('frame_selector', type=str,
                    choices=frame_selectors.keys(),
                    help='Strategy for selecting frames to represent each MSM bin.')
parser.add_argument('align_selection', type=str,
                    help='A file containing a loos selection string to align the selected frames with. Alternatively '
                         'provide a selection string on the command line.')
parser.add_argument('traj_paths', type=str, nargs='+',
                    help='A file containing a list of trajectories corresponding (in matching order) to the supplied '
                         'assignments file. Alternatively, paths to the trajectory files.')
# Optional args below here
parser.add_argument('--subset-selection', type=str, default='all',
                    help='A loos selection string to subset all frames by (for example, if water is present it could be '
                         'stripped here).')

parser.add_argument('--number-frames', default=10, type=int,
                    help='Number of frames to select per-bin. If a bin has fewer total assignments than this value, '
                         'an error is thrown.')
parser.add_argument('--total-per-bin', default=None, type=str,
                    help='Text file with totals to draw from each MSM bin. Should be one column of totals, '
                         'where the row index corresponds to the bin and the entry is the number of frames to draw.')
parser.add_argument('--align-resid-list', type=str, default=None,
                    help='If provided, use numbers in file as a list of residue IDs. Concatenate align selection string '
                         'with one selecting these resids.')
parser.add_argument('--make-receptor-sel-chain-A', action=argparse.BooleanOptionalAction,
                    help='If thrown, make all atoms a member of chain "A" when writing PDBs.')
parser.add_argument('--write-bin-trajs', action=argparse.BooleanOptionalAction,
                    help='If thrown, write a DCD with the selected frames in each bin directory.')

if __name__ == '__main__':
    args = parser.parse_args()
    assignments = ra.load(args.assignments)
    try:
        eq_probs = np.load(args.eq_probs)
    except ValueError:
        eq_probs = np.load(args.eq_probs, allow_pickle=True).item().eq_probs_

    model = loos.createSystem(args.model)
    if args.make_receptor_sel_chain_A:
        for atom in model:
            atom.chainId('A')
    align_sel = ''
    if args.align_resid_list:
        alresids = np.genfromtxt(args.align_resid_list)
        align_sel += ' || '.join(['resid == {}'.format(i) for i in alresids])

    # if filename is provided read its contents to obtain align string
    try:
        align_sel += Path(args.align_selection).read_text()
    except FileNotFoundError:  # if the string is not a file name, interpret it as a loos selection string.
        align_sel += args.align_selection
    if args.total_per_bin:
        total_per_bin = np.loadtxt(args.total_per_bin)
        chosen_frames = frame_selectors[args.frame_selector](
            assignments,
            eq_probs.shape[0],
            total_per_bin
        )
    else:
        low_inds = np.where(
            np.bincount(assignments.flatten()) < args.number_frames)
        if len(low_inds[0]) > 0:
            print(
                'The(se) bin(s) have fewer than your requested samples in them:')
            print(low_inds[0])
            print('You requested this many frames per bin:')
            print(args.number_frames)
            print('Exiting.')
            exit(1)
        chosen_frames = frame_selectors[args.frame_selector](
            assignments,
            eq_probs.shape[0],
            args.number_frames)

    if len(args.traj_paths) == 1:
        p_trjs = Path(args.traj_paths[0])
        if p_trjs.suffix == '.txt':
            traj_paths = p_trjs.read_text().split()
        else:
            traj_path = p_trjs
    else:
        traj_paths = args.traj_paths
    print('aligning with the following selection string:')
    print(align_sel)
    out_path = Path(args.receptor_name) / 'receptor'
    full_inds, subset_vec, align_vec = rip_conformations(
        chosen_frames, model, args.subset_selection, align_sel, traj_paths
    )
    align_samples(subset_vec, align_vec)
    write_sampled_frames(subset_vec, full_inds, out_path, args.write_bin_trajs)


