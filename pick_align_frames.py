import argparse
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


def rip_conformations(chosen_inds, model, subset_selection, align_selection, out_path, traj_paths):
    trajectories = [pyloos.Trajectory(traj_path, model, subset=subset_selection) for traj_path in traj_paths[1:]]
    full_inds = []
    subset_list = []
    align_list = []
    bin_trj_file_names = [loos.DCDWriter(out_path+'/'+str(i)+'/samples.dcd') for i in range(chosen_inds.shape[0])]
    # cflat = chosen_inds.flatten()
    # chosen_sort_inds = np.argsort(cflat, order=['traj', 'frame'])
    # inverse_sort_inds = np.zeros_like(chosen_sort_inds, dtype=np.int32)
    # inverse_sort_inds[chosen_sort_inds] = np.arange(inverse_sort_inds.shape[0], dtype=np.int32)

    # read selected frames into memory
    for bin_ix, samples in enumerate(chosen_inds):
        for trj_ix, fra_ix in samples:
            frame = trajectories[trj_ix].readFrame(fra_ix)
            align_list.append(loos.selectAtoms(frame, align_selection))
            subset_list.append(frame)
            full_inds.append((bin_ix, trj_ix, fra_ix))

    xforms, rmsd, iters = pyloos.iterativeAlignment(align_list)
    # apply transforms here
    loos.applyTransforms(subset_list, xforms)
    print('Iterative alignment final RMSD', rmsd, 'after', iters, 'iterations')

    # now write everything out
    for (bin_ix, trj_ix, fra_ix), frame in zip(full_inds, subset_list):
        bin_trj_file_names[bin_ix].writeFrame(frame)
        pdb = loos.PDB.fromAtomicGroup(frame)
        pdbfname =out_path + '/{}/{}-{}.pdb'.format(bin_ix, trj_ix, fra_ix)
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
parser.add_argument('eq_probs', type=str,
                    help='.npy file with equilibrium probabilities from MSM.')
parser.add_argument('traj_paths', type=str,
                    help='A file containing a list of trajectories corresponding (in matching order) to the supplied '
                         'assignments file.')
parser.add_argument('align_selection', type=str,
                    help='A loos selection string to align the selected frames with.')
parser.add_argument('--subset-selection', type=str, default='all',
                    help='A loos selection string to subset all frames by (for example, if water is present it could be stripped here).')
parser.add_argument('frame_selector', type=str,
                    choices=frame_selectors.keys())
parser.add_argument('--frames-per-bin', type=int, default=10,
                    help='Number of frames to select per-bin. If a bin has fewer total assignments than this value, '
                         'an error is thrown. Overridden by "--proportional".')
# parser.add_argument('--proportional', action=argparse.BooleanOptionalAction,
#                   help='Compute number of frames to use for each state as a function of some heuristic scoring frame '
#                          'diversity. Value provided is minumum number of frames to take from least diverse state.')

if __name__ == '__main__':
    argv = [
        "tmh2",
        "/home/louis/bowmanlab/j.lotthammer/Simulations/myosin/specificity/pps-bleb-isoforms/myh2/input_files3/myh2-5n6a-holo-prot-masses.pdb",
        '~/myosin/traj_files.txt'
    ]
    args = parser.parse_args()
    assignments = ra.load(args.assignments)
    eq_probs = np.load(args.eq_probs)
    model = loos.createSystem(args.model)
    chosen_frames = frame_selectors[args.frame_selector](
        assignments,
        eq_probs.shape[0],
        args.frames_per_bin)
    traj_paths = read_traj_paths(args.traj_paths)
    rip_conformations(chosen_frames, model, args.subset_selection, args.align_selection, args.receptor_name+'/receptor',
                      traj_paths)

