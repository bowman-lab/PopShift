import argparse
import numpy as np
from enspara import ra


def get_assign_inds(assignments, nstates):
    traj_ind_list = [[] for i in range(nstates)]
    frame_ind_list = [[] for i in range(nstates)]
    for i, t in enumerate(assignments):
        for j, f in enumerate(t):
            # iterating over RA assign contents produces singleton arrays of dtype=int32
            # that need to be unpacked for use
            ind = f[0]
            traj_ind_list[ind].append(i)
            frame_ind_list[ind].append(j)
    return [np.array([t_inds, f_inds], dtype=np.int32)
            for t_inds, f_inds in zip(traj_ind_list, frame_ind_list)]


def get_random_per_bin(assignments, n_states, number_desired, replace=False, gen=np.random.default_rng()):
    # assumes RA's traj axis is the outer one.S
    all_state_indices = get_assign_inds(assignments, n_states)
    selected_inds = np.array([gen.choice(state_ixs, number_desired, axis=1, replace=replace)
                              for state_ixs in all_state_indices], dtype=[('bin', np.int32), ('traj', np.int32), ('frame', np.int32)])
    return selected_inds
chosen = get_random_per_bin(d, 200, 3, gen=gen)

def read_traj_paths(traj_list_filename):
    with open(traj_list_filename) as f:
        return f.readlines()

frame_selectors = {
    'random': get_random_per_bin,
}

parser = argparse.ArgumentParser()
parser.add_argument('assignments', type=str,
                    help='h5 file with assignments associated to the MSM used.')
parser.add_argument('eq_probs', type=str,
                    help='.npy file with equilibrium probabilities from MSM.')
parser.add_argument('traj_paths', type=str,
                    help='A file containing a list of trajectories corresponding (in matching order) to the supplied '
                         'assignments file.')
parser.add_argument('model', type=str,
                    help='A loos-interpretable model file that will permit reading the trajectories in "traj_paths".')
parser.add_argument('selection', type=str,
                    help='A loos selection string to align the selected frames with.')
parser.add_argument('frame_selector', type=str,
                    choices=frame_selectors.keys())
parser.add_argument('--frames-per-bin', type=int, default=10,
                    help='Number of frames to select per-bin. If a bin has fewer total assignments than this value, '
                         'an error is thrown. Overridden by "--proportional".')
# parser.add_argument('--proportional', action=argparse.BooleanOptionalAction,
#                   help='Compute number of frames to use for each state as a function of some heuristic scoring frame '
#                          'diversity. Value provided is minumum number of frames to take from least diverse state.')

if __name__ == '__main__':
    args = parser.parse_args()
    assignments = ra.load(args.assignments)
    eq_probs = np.load(args.eq_probs)
    chosen_frames = frame_selectors[args.frame_selector](
        assignments,
        eq_probs.shape[0],
        args.frames_per_bin)
    traj_paths = read_traj_paths(args.traj_paths)
    for