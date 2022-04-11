import argparse
from pathlib import Path
import numpy as np
from enspara import ra
import loos
from loos import pyloos


def get_assign_inds(assignments, nstates):
    index_list_by_bin = [[] for i in range(nstates)]
    for i, t in enumerate(assignments):
        for j, f in enumerate(t):
            # iterating over RA assign contents produces singleton arrays of dtype=int32
            # that need to be unpacked for use
            bin_index = f[0]
            index_list_by_bin[bin_index].append((i, j))
    return [np.array(bin_indices, dtype=[('traj', np.int32), ('frame', np.int32)]) for bin_indices in index_list_by_bin]


def get_random_per_bin(assignments, n_states, number_desired, replace=False, gen=np.random.default_rng()):
    # assumes RA's traj axis is the outer one.S
    all_state_indices = get_assign_inds(assignments, n_states)
    chosen_inds = np.array([gen.choice(state_ixs, number_desired, axis=0, replace=replace)
                            for state_ixs in all_state_indices])
    return chosen_inds


def get_specific_number_per_bin(assignments, n_states, number_desired, replace=False, gen=np.random.default_rng()):
    # assumes RA's traj axis is the outer one.S
    all_state_indices = get_assign_inds(assignments, n_states)
    chosen_inds = np.array([gen.choice(state_ixs, number_desired[i], axis=0, replace=replace)
                            for i, state_ixs in enumerate(all_state_indices)])
    return chosen_inds

def rip_conformations(chosen_inds, model, subset_selection, align_selection, out_path, traj_paths):
    trajectories = [pyloos.Trajectory(traj_path.rstrip(), model, subset=subset_selection) for traj_path in traj_paths]
    full_inds = []
    subset_list = loos.AtomicGroupVector()
    align_list = loos.AtomicGroupVector()
    bin_trajs = []

    # read selected frames into memory
    for bin_ix, samples in enumerate(chosen_inds):
        p = Path(out_path + '/' + str(bin_ix))
        p.mkdir(parents=True, exist_ok=True)  # This will operate like mkdir -pf; it will overwrite.
        bin_trajs.append(loos.DCDWriter(out_path + '/' + str(bin_ix) +
                                        '/samples.dcd'))
        for trj_ix, fra_ix in samples:
            print(trj_ix, fra_ix)
            frame = trajectories[trj_ix].readFrame(fra_ix)
            align_list.push_back(loos.selectAtoms(frame, align_selection).copy())
            subset_list.push_back(frame.copy())
            full_inds.append((bin_ix, trj_ix, fra_ix))

    alignment_result = loos.iterativeAlignmentPy(align_list)
    # apply transforms here
    loos.applyTransforms(subset_list, alignment_result.transforms)
    print('Iterative alignment final RMSD', alignment_result.rmsd, 'after', alignment_result.iterations, 'iterations')

    # now write everything out
    for (bin_ix, trj_ix, fra_ix), frame in zip(full_inds, subset_list):
        bin_trajs[bin_ix].writeFrame(frame)
        pdb = loos.PDB.fromAtomicGroup(frame)
        pdbfname = out_path + '/{}/{}-{}.pdb'.format(bin_ix, trj_ix, fra_ix)
        with open(pdbfname, "w") as f:
            f.write(str(pdb))


def read_traj_paths(traj_list_filename):
    with open(traj_list_filename) as f:
        return f.readlines()


frame_selectors = {
    'random': get_random_per_bin,
}

parser = argparse.ArgumentParser()
parser.add_argument('receptor_name', type=str,
                    help="Name to use to use as top-level directory for the docking run tree.")
parser.add_argument('model', type=str,
                    help='A loos-interpretable model file that will permit reading the trajectories in "traj_paths".')
parser.add_argument('assignments', type=str,
                    help='h5 file with assignments associated to the MSM used.')
parser.add_argument(metavar='eq_probs | pickled_msm', type=str, dest='eq_probs',
                    help='.npy file with equilibrium probabilities from MSM, or pickled MSM object.')
parser.add_argument('traj_paths', type=str,
                    help='A file containing a list of trajectories corresponding (in matching order) to the supplied '
                         'assignments file.')
parser.add_argument('align_selection', type=str,
                    help='A file containing a loos selection string to align the selected frames with. Alternatively provide a selection string on the command line.')
parser.add_argument('--subset-selection', type=str, default='all',
                    help='A loos selection string to subset all frames by (for example, if water is present it could be stripped here).')
parser.add_argument('frame_selector', type=str,
                    choices=frame_selectors.keys())
parser.add_argument('--frames-per-bin', type=int, default=10,
                    help='Number of frames to select per-bin. If a bin has fewer total assignments than this value, '
                         'an error is thrown. Overridden by "--proportional".')
parser.add_argument('--align-resid-list', type=str, default=None,
                    help='If provided, use numbers in file as a list of residue IDs. Concatenate align selection string '
                         'with one selecting these resids.')
parser.add_argument('--make-receptor-sel-chain-A', action=argparse.BooleanOptionalAction,
                    help='If thrown, make all atoms a member of chain "A" when writing PDBs.')
# parser.add_argument('--proportional', action=argparse.BooleanOptionalAction,
#                   help='Compute number of frames to use for each state as a function of some heuristic scoring frame '
#                          'diversity. Value provided is minumum number of frames to take from least diverse state.')

if __name__ == '__main__':
    argv = [
        "--subset-selection",
        "resid < 793",
        "--make-receptor-sel-chain-A",
        "--frames-per-bin",
        "2",
        "tmh2",
        "/home/louis/testtrajs/myh2-5n6a-holo-prot-masses.pdb",
        "/home/louis/myosin/myh2-5n6a-backbone-chi1-dihedrals-bleb-pocket-correlated-subset-charmm36-sims-lag-500-tica-reduced-k-200-cluster-dtrajs.h5",
        "/home/louis/myosin/myh2-5n6a_input_features_backbone-chi1-dihedrals-bleb-pocket-correlated-subset-charmm36-sims_tica_lag_500_k_200.npy",
        '/home/louis/myosin/traj_files.txt',
        "/home/louis/myosin/myh2-5n6a-pocket-sel.txt",
        "random"
    ]
    args = parser.parse_args(argv)
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
    chosen_frames = frame_selectors[args.frame_selector](
        assignments,
        eq_probs.shape[0],
        args.frames_per_bin)
    traj_paths = read_traj_paths(args.traj_paths)
    print('aligning with the following selection string:')
    print(align_sel)
    rip_conformations(chosen_frames, model, args.subset_selection, align_sel, args.receptor_name+'/receptor',
                      traj_paths)

