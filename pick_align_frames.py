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
    return [np.array([np.array(tinds, dtype=np.int32), np.array(finds, dtype=np.int32)])
            for tinds, finds in zip(traj_ind_list, frame_ind_list)]


def get_random_per_bin(assignments, n_states, number_desired, replace=False, gen=np.random.default_rng()):
    # assumes RA's traj axis is the outer one.S
    all_state_indices = get_assign_inds(assignments, n_states)
    selected_inds = [np.choice(state_ixs, number_desired, axis=1) for state_ixs in all_state_indices]
    return selected_inds

frame_selectors = {
    'random': get_random_per_bin,
}

parser = argparse.ArgumentParser()
parser.add_argument('assignments', type=str,
                    help='h5 file with assignments associated to the MSM used.')
parser.add_argument('eq_probs', type=str,
                    help='npy file with equilibrium probabilities for MSM')
parser.add_argument('frame_selector', type=str,
                    choices=frame_selectors.keys())
parser.add_argument('--frames-per-bin', type=int, default=10,
                    help='Number of frames to select per-bin. If a bin has fewer total assignments than this value, '
                         'and error is thrown. Overridden by "--proportional".')
parser.add_argument('--proportional', action=argparse.BooleanOptionalAction,
                    help='Compute number of frames to use for each state as a function of some heuristic scoring frame '
                         'diversity. Value provided is minumum number of frames to take from least diverse state.')