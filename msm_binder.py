import numpy as np
from enspara import ra
from enspara import msm
from pathlib import Path
import json
# import stdout to dump json, import argv to save command line
from sys import stdout, argv
from itertools import repeat
# import pyemma.coordinates as coor  # note this is done below only in the case where it's needed.

import argparse


# helper to condense list appending during input file reading.
def parse_indexed_fe_line(traj_index_list, frame_index_list, fes_list, line):
    split_line = line.split()
    traj_index_list.append(int(split_line[0]))
    frame_index_list.append(int(split_line[1]))
    fes_list.append(float(split_line[2]))


# parse a three-column file giving the trajectory number, frame number, and
# free energy for each binding attempt in a replica
def parse_indexed_fe_file(filename, stride=1):
    traj_indices = []
    frame_indices = []
    scores = []
    with open(filename, 'r') as f:
        line = f.readline()
        # check for a comment
        if "#" in line:
            pass
        else:
            parse_indexed_fe_line(traj_indices, frame_indices, scores, line)
        for line in f:
            parse_indexed_fe_line(traj_indices, frame_indices, scores, line)
    # bundle the coordinates together in a tuple for indexing of RAs.
    # include stride at this point. This could lead to errors if operator
    # 'strides' both when recording indexed scores and also when calling this tool.
    return (np.array(traj_indices, dtype=np.int32), np.array(frame_indices, dtype=np.int32) * stride), \
           np.array(scores)


def process_indexed_fe_file(indexed_fe_filename, assignments, active_states, eq_probs, mapping, stride=1):
    inds, scores = parse_indexed_fe_file(indexed_fe_filename, stride=stride)

    mask = np.isin(assignments[inds], active_states)
    trimmed_fes = scores[mask]
    trimmed_filtered_assigns = assignments[inds][mask]
    state_counts = np.bincount(trimmed_filtered_assigns)
    # scale eqch frame weight by the number of counts in each state
    frame_weights = (eq_probs[mapping] / state_counts)[trimmed_filtered_assigns]
    return frame_weights, trimmed_fes



# Determine weighting for each frame
# weight = (eq prob of state / number of frames in state)
def filter_frame_weights(msm_obj, stride, state_counts, ra_assigns, active_states):
    return np.array(
        [msm_obj.eq_probs_[msm_obj.mapping_.to_mapped[s]] / state_counts[msm_obj.mapping_.to_mapped[s]]
         for traj in ra_assigns
         for s in traj[::stride]
         if s in active_states
         ]
    )

# determine weights for each sample-set taken from a given bin
# takes an array of bin weights and an array of lengths (as from RaggedArray.lengths)
# returns RaggedArray of eq_prob per bin divided by number of samples drawn from that bin.
def expand_bin_weights(eq_probs, lengths):
    return ra.RaggedArray(
        [np.array([p for p in repeat(eq_probs[i]/length)]) for i, length in lengths]
    )


# takes a ragged array of scores, a ragged array of
def filter_trim_binding_fes(binding_fes, active_states, stride, ra_assigns):
    trimmed_binding_fes = np.array(
        [fe for i, binding_fe in enumerate(binding_fes)
         for j, fe in enumerate(binding_fe)
         if ra_assigns[i][j * stride] in active_states]
    )
    return trimmed_binding_fes


# returns an array of active states
def get_active_states(msm_obj):
    return np.array(list(msm_obj.mapping_.to_original.values()))


# returns state counts given striding, and ergodic trimming
def count_strided_states(ra_assigns, stride, active_states):
    return np.bincount([s for traj in ra_assigns for s in traj[::stride]])[active_states]


# compute reweightings.
def reweighted_frames(frame_weights, kas):
    unnormed_weights = frame_weights * (1 + kas)
    return unnormed_weights / np.sum(unnormed_weights)


def free_energy_per_state(frame_weights, reweights, rt):
    state_eq = reweights / frame_weights
    return rt * np.log(state_eq)


def kd_from_kcal_mol(fe, rt):
    return np.exp(fe / rt)

def kcal_mol_from_kd(kd, rt):
    return np.log(kd * rt)


# reductions below here
def msm_binding_dG(frame_weights, trimmed_binding_fes, rt):
    rtn = -rt
    return rtn * np.log(np.sum(frame_weights * np.exp(trimmed_binding_fes / rtn)))


def weighted_avg(frame_weights, trimmed_binding_fes):
    return np.average(trimmed_binding_fes, weights=frame_weights)


def simple_avg(trimmed_binding_fes):
    return np.mean(trimmed_binding_fes)


def calx_store(binding_output, trimmed_fes, frame_weights, rt, tag, kd_scale, reweighted_eq_prefix):
    msm_binding = msm_binding_dG(frame_weights, trimmed_fes, rt)
    kd = kd_from_kcal_mol(msm_binding, rt) * kd_scale  # kd, scaled by user-supplied conversion.
    if reweighted_eq_prefix:
        # convert trimmmed Free energies to association constants
        kas = kd_from_kcal_mol(trimmed_fes, rt)**(-1)
        reweights = reweighted_frames(frame_weights, kas)
        fe_per_state = free_energy_per_state(frame_weights, reweights, rt)
        np.save(reweighted_eq_prefix+tag+'-fe.npy', fe_per_state)
        np.save(reweighted_eq_prefix+tag+'-eq_probs.npy', reweights)
    weighted = weighted_avg(frame_weights, trimmed_fes)
    simple = simple_avg(trimmed_fes)
    binding_output[tag] = {
        'msm dG': msm_binding,
        'msm K_D': kd,
        'weighted avg': weighted,
        'simple avg': simple
    }


def interp_trj_samples(args, rt, binding_output):
    if args.emma_dtraj:  # note this import here allows the tool to not strictly depend on PyEMMA.
        from pyemma import coordinates as coor
        assignments = coor.load(args.assignments)
        for i, trj in enumerate(assignments):
            assignments[i] = trj.astype(int).flatten()

        assignments = ra.RaggedArray(assignments)
    else:
        assignments = ra.load(args.assignments)

    # this may need to be redone to be more compatible with PyEMMA or other builders.
    # in principle all that's needed is an eq-probs array and a mapping array that works like enspara's.
    # check to see if they've fed us some weird pickled MSM
    if '.npy' in args.msm:
        msm_obj = np.load(args.msm, allow_pickle=True).item()
    else:  # assume they used enspara's save method.
        msm_obj = msm.MSM.load(args.msm)
    # subset of states that were not trimmed
    active_states = get_active_states(msm_obj)
    if not args.index_from_file:
        state_counts = count_strided_states(assignments, args.stride, active_states)
        frame_weights = filter_frame_weights(msm_obj, args.stride, state_counts, assignments, active_states)
    else:
        mapping = np.array(
            [msm_obj.mapping_.to_mapped[k] for k in msm_obj.mapping_.to_mapped.keys()],
            dtype=np.int32
        )[active_states]
    for binding_run in args.binding_fes:
        tag = Path(binding_run).stem  # turn base filename no ext into tag for saving later.
        if args.index_from_file:
            frame_weights, trimmed_fes = process_indexed_fe_file(binding_run, assignments,
                                                                 active_states, msm_obj.eq_probs_, mapping, stride=args.stride)
        else:
            fes = ra.RaggedArray(np.load(binding_run, allow_pickle=True))
            trimmed_fes = filter_trim_binding_fes(fes, active_states, args.stride, assignments)

        calx_store(binding_output, trimmed_fes, frame_weights, rt, tag, args.K_D_scale, args.reweighted_eq_prefix)




def interp_bin_samples(args, rt, binding_output):
    eq_probs = args.eq_probs
    for binding_run in args.bindingFE_h5s:
        fes = ra.load(binding_run)
        tag = Path(binding_run).stem
        sample_weights = expand_bin_weights(eq_probs, fes.lengths)
        calx_store(binding_output, fes, sample_weights, rt, tag, args.K_D_scale, args.reweighted_eq_prefix)




def run_cli(raw_args=None):
    # relevant constants
    R = 1.98720425864083  # cal*k^-1mol^-1, wikipedia table
    T = 310.0  # Kelvin
    unit_scale = 0.001  # by default convert free energy to kilo-energy units.

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='Select mode of input. Either samples from trajectories, '
                                       'or samples from each bin in an MSM.')
    parser.add_argument('msm', type=str,
                        help='Name of MSM object (mapping expected) to use in calculation.')

    # supra parser optional args
    parser.add_argument('--append-to', type=str, default=None,
                        help='If the name of a JSON is provided, will append results to that object. '
                             'By default, writes a new one from scratch.')
    parser.add_argument('--out', '-o', type=str, default=None,
                        help='If name is provided, overwrite to it. Otherwise write to stdout. '
                             'Superceded by "--apend-to".')
    parser.add_argument('--gas-constant', '-R', type=float, default=R,
                        help="The ideal gas constant. Choose to match units of T.")
    parser.add_argument('--temperature', '-T', type=float, default=T,
                        help='The temperature to estimate free energy at. Choose to match units of R.')
    parser.add_argument('--unit-scale', type=float, default=unit_scale,
                        help='Scale output free energy by this value. '
                             'Default scales to kilo-energy units (4ex kcals/mol).')
    parser.add_argument('--K_D-scale', type=float, default=10 ** 6,
                        help='Scale dissoc. constant by this value to put entry in more customary range. '
                             'Default converts to micromolar.')
    parser.add_argument('--reweighted-eq-prefix', '-w', type=str, default=None,
                        help='Prefix for filename to write reweighted equilibrium probabilities '
                             'for each state to. Output will be in numpy binary format')

    trj_parser = subparsers.add_parser('traj-samples', help='Use a sequence of docked frames. '
                                                            'Use assigns and MSM to map these scores to MSM bins '
                                                            '(and therefore equilibrium probabilities).')
    trj_parser.set_defaults(func=interp_trj_samples)
    # required positional arg for traj reading
    trj_parser.add_argument('assignments', type=str,
                            help='Name of discretized trajectory.')
    # optional args to control traj reading
    trj_parser.add_argument('--emma-dtraj', '-e', action=argparse.BooleanOptionalAction,
                            help='If thrown, read discretized trajectory as a PyEMMA dtraj.')
    trj_parser.add_argument('binding_fes', type=str, nargs="+",
                            help='Name of binding FE file for some compound. Expects .npy by default.')
    trj_parser.add_argument('--stride', '-s', type=int, default=1,
                            help='Stride-rate through dtrajs. Applied uniformly to frames. If using indexed score '
                                 'input, and indices match those of dtrajs (already strided), set to 1.')
    trj_parser.add_argument('--index-from-file', action=argparse.BooleanOptionalAction,
                            help='If thrown, interprets binding_fes file(s) as list of RA indexes, '
                            'with the last element containing the binding FE/score. Stride is applied to frames still.')

    bin_parser = subparsers.add_parser('bin-samples')
    bin_parser.set_defaults(func=interp_bin_samples)
    bin_parser.add_argument('eq_probs', type=lambda x: np.load(x),
                            help='Path to numpy array of equilibrium probabilities in 1-1 correspondence to the states '
                                 'in each FE ragged array.')
    bin_parser.add_argument('bindingFE_h5s', nargs='+',
                            help='File name(s) of extracted binding scores. Should be ragged arrays of lengths '
                                 'msm_bins, samples_from_bin.')

    args = parser.parse_args(raw_args)

    rt = args.gas_constant * args.temperature * args.unit_scale
    outfn = args.out  # Note, this defaults to none.
    if args.append_to:
        with open(args.append_to) as f:
            binding_output = json.load(f)
        outfn = args.append_to
    else:
        binding_output = {}

    try:
        binding_output['command line'].append(argv)
    except KeyError:
        binding_output['command line'] = []
        binding_output['command line'].append(argv)

    args.func(args, rt, binding_output)

    binding_output['log'] = {}
    binding_output['log']['rt'] = rt

    if outfn:
        with open(outfn, 'w') as f:
            json.dump(binding_output, f, indent=4)
    else:
        json.dump(binding_output, stdout, indent=4)


if __name__ == '__main__':
    run_cli(raw_args=argv[1:])
