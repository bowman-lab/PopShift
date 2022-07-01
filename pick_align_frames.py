import argparse
from pathlib import Path
import numpy as np
from enspara import ra
from enspara.util.load import concatenate_trjs
import loos
from loos import pyloos
import json
from sys import argv
from scipy.spatial.distance import euclidean
from sklearn.cluster import KMeans
import pickle


def floatpair(s, delim=','):
    try:
        start_ind, end_ind = map(float, s.split(delim))
        return start_ind, end_ind
    except:
        raise argparse.ArgumentTypeError('Float pairs must be specified as "decimal.one,decimal.two"')


def get_assigns_no_map_no_unbox(assignments, nstates):
    print('Getting assignments for each bin without any mapping or unboxing.')
    index_list_by_bin = [[] for i in range(nstates)]
    # handle the case where there is only 1 trajectory,
    # meaning the ragged array doesn't have singleton elements at the tightest level.
    for traj_index, dtraj in enumerate(assignments):
        for frame_index, cluster_index in enumerate(dtraj):
            index_list_by_bin[cluster_index].append((traj_index, frame_index))
    return index_list_by_bin


def get_assigns_no_map_unbox(assignments, nstates):
    print('Getting assignments for each bin without any mapping. Unbox singletons at deepest layer of RA.')
    index_list_by_bin = [[] for i in range(nstates)]
    # iterating over RA assign contents produces singleton arrays of dtype=int32
    # that need to be unpacked for use
    for traj_index, dtraj in enumerate(assignments):
        for frame_index, cluster_index in enumerate(dtraj):
            index_list_by_bin[cluster_index].append((traj_index, frame_index))
    return index_list_by_bin


def get_assigns_map_no_unbox(assignments, nstates, mapping):
    print('Getting assignments for each bin from a map, no unboxing')
    # use try statement to filter states that have been trimmed
    index_list_by_bin = [[] for i in range(nstates)]
    for traj_index, dtraj in enumerate(assignments):
        for frame_index, cluster_index in enumerate(dtraj):
            try:
                bin_index = mapping[cluster_index]
                index_list_by_bin[bin_index].append((traj_index, frame_index))
            except KeyError:
                pass
    return index_list_by_bin


def get_assigns_map_unbox(assignments, nstates, mapping):
    # use try statement to filter states that have been trimmed
    # unbox the singleton cluster index arrays
    print('Getting assignments for each bin from a map, unboxing singletons at deepest layer of RA.')
    index_list_by_bin = [[] for i in range(nstates)]
    for traj_index, dtraj in enumerate(assignments):
        for frame_index, cluster_index in enumerate(dtraj):
            try:
                bin_index = mapping[cluster_index[0]]
                index_list_by_bin[bin_index].append((traj_index, frame_index))
            except KeyError:
                pass
    return index_list_by_bin


def get_assign_inds(assignments, nstates, mapping=None):
    # Handle the case where there is only 1 trajectory,
    # meaning the ragged array doesn't have singleton elements at the tightest level.
    if len(assignments.shape) == 1:
        if mapping:
            index_list_by_bin = get_assigns_map_no_unbox(assignments, nstates, mapping)
        else:
            index_list_by_bin = get_assigns_no_map_no_unbox(assignments, nstates)
    # Handle the case where the last element in the shape of the RA is 1 (when it should be 'None').
    # This is interpreted as having singleton elements at the tightest level, and thus unboxing them.
    elif assignments.shape[-1] == 1:
        if mapping:
            index_list_by_bin = get_assigns_map_unbox(assignments, nstates, mapping)
        else:
            index_list_by_bin = get_assigns_no_map_unbox(assignments, nstates)
    else:  # operate on a normal 2-ragged array.
        if mapping:
            index_list_by_bin = get_assigns_map_no_unbox(assignments, nstates, mapping)
        else:
            index_list_by_bin = get_assigns_no_map_no_unbox(assignments, nstates)
    return [np.array(bin_indices, dtype=[('traj', np.int32), ('frame', np.int32)]) for bin_indices in index_list_by_bin]


def get_random_per_bin(assignments, n_states, number_desired, mapping, replace=False, gen=np.random.default_rng()):
    # assumes RA's traj axis is the outer one.S
    all_state_indices = get_assign_inds(assignments, n_states, mapping)
    chosen_inds = np.array([gen.choice(state_ixs, number_desired, axis=0, replace=replace)
                            for state_ixs in all_state_indices])
    return chosen_inds


def get_specified_number_per_bin_random(assignments, nstates, numbers_desired, mapping, replace=False,
                                        gen=np.random.default_rng()):
    all_state_indices = get_assign_inds(assignments, nstates, mapping)
    chosen_inds = (gen.choice(state_ixs, number_desired, axis=0, replace=replace)
                   for state_ixs, number_desired in zip(all_state_indices, numbers_desired))
    return list(chosen_inds)


def calc_euclidean_distance(reference, other_frames_in_the_bin):
    distances = [euclidean(reference, other_frames_in_the_bin[num]) for num in range(len(other_frames_in_the_bin))]
    return distances


def map_features(assignments, nstates, mapping, features):
    mapped_features = []
    all_state_indices = get_assign_inds(assignments, nstates, mapping)
    for bin in all_state_indices:
        bin_features = []
        for frame in bin:
            bin_features.append(features[frame[0]][frame[1]])
        mapped_features.append(bin_features)
    return (ra.RaggedArray(mapped_features), all_state_indices)


def kmeans_cluster(features, n_clusters):
    kmeans = KMeans(n_clusters=n_clusters).fit(features)
    return kmeans.cluster_centers_


def find_closest_frame(cluster_centers, features, indices):
    all_indices = []
    for num, bin in enumerate(cluster_centers):
        distances = [calc_euclidean_distance(frame, features[num]) for frame in bin]
        min_distance = [np.min(dist) for dist in
                        distances]  # ! Not sure if i want/need this, might be interesting to print out
        min_distance_loc = [np.where(dist == np.min(dist)) for dist in distances]
        if len(min_distance_loc[0]) > 1:
            min_distance_loc = np.random.default_rng().choice(min_distance_loc[0], 1, axis=0, replace=False)[0]
        frames_to_extract = [indices[num][min_distance_loc[i][0]] for i in range(len(min_distance_loc))]
        all_indices.append(frames_to_extract)
    return all_indices


def get_frames_using_kmeans(assignments, nstates, features, n_clusters, mapping):
    print('Mapping features to bins')
    mapped_features, all_state_indices = map_features(assignments, nstates, mapping, features)
    print('KMeans clustering of the frames within bins using features provided')
    # ! This line seems slow for some reason (not slow in my jupyter nb)
    cluster_centers = [kmeans_cluster(bin_features, number) for bin_features, number in zip(mapped_features,
                                                                                            n_clusters)]
    print('Looking for frames closest to each kmeans center')
    indices_to_extract = find_closest_frame(cluster_centers, mapped_features,
                                            all_state_indices)  # ! Use multiprocessing here
    indices_to_extract = [np.concatenate(i) for i in indices_to_extract]
    print('Done')
    return indices_to_extract


def get_centers(assigs, nstates, nd, mapping):
    return np.array([[0, i] for i in range(nstates)]).reshape(nstates, 1, 2)


def unpickle_resave_centers(centersfn):
    cens = concatenate_trjs(np.load(centersfn, allow_pickle=True))
    fstem = Path(centersfn).stem
    newname = fstem + '.xtc'
    cens.save(newname)
    return newname


def rip_conformations(chosen_inds, model, subset_selection, align_selection, traj_paths):
    # This'll be a flat array that'll include everything needed to index a frame.
    full_inds = np.array([(bin_ix, trj_ix, fra_ix) for bin_ix, samples in enumerate(chosen_inds)
                          for trj_ix, fra_ix in samples],
                         dtype=[('bin', int), ('traj', int), ('frame', int)])
    sorted_chosen = np.argsort(full_inds, order=['traj', 'frame', 'bin'])

    # instantiate loos vectors to do iterative alignment.
    subset_vec = loos.AtomicGroupVector(len(full_inds))
    align_vec = loos.AtomicGroupVector(len(full_inds))

    # set up the atomic groups that readFrame will update (they'll be references)
    subset = loos.selectAtoms(model, subset_selection)
    align_subset = loos.selectAtoms(model, align_selection)
    # track previously read-from traj to see if we need to change which one is open.
    prev_trj = None
    # loop over the selected frame indices.
    for sort_ix in sorted_chosen:
        _, trj_ix, fra_ix = full_inds[sort_ix]
        trj_name = traj_paths[trj_ix]
        if trj_name != prev_trj:
            trj = pyloos.Trajectory(str(trj_name), model,
                                    subset=subset_selection)
        prev_trj = trj_name
        trj.readFrame(fra_ix)
        # For now, np.int64 is not an acceptable index type for AGVecs.
        retyped_sort_ix = int(sort_ix)
        # readFrame updates ref to the traj's internal AG,
        # so both product AGs need to be copied.
        align_vec[retyped_sort_ix] = align_subset.copy()
        subset_vec[retyped_sort_ix] = subset.copy()

    return full_inds, subset_vec, align_vec


def align_samples(subset_vec, align_vec):
    # iteratively aligns AGs in group. Alignment result contains transforms to
    # rotate/translate the aligned group into alignment with the target.
    alignment_result = loos.iterativeAlignmentPy(align_vec)
    # apply transforms here. This is done in-place
    loos.applyTransforms(subset_vec, alignment_result.transforms)
    print('Iterative alignment final RMSD', alignment_result.rmsd, 'after',
          alignment_result.iterations, 'iterations')


def write_sampled_frames(subset_list, full_inds, out_path, write_bin_trajs: bool):
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


def add_bonds_two_cuts(model: loos.AtomicGroup, heavy_cutoff: float, hydrogen_cutoff: float):
    # select just the heavy atoms from the model
    heavies = loos.selectAtoms(model, '!hydrogen')
    # populate the bonds in the heavy-atom only file using the heavy cutoff
    heavies.findBonds(heavy_cutoff)
    # now populate all the bonds to hydrogen with the hydrogen cutoff.
    model.findBonds(hydrogen_cutoff)
    # return a new atomic group that merges the two separate groups.
    return model.merge(heavies)


def pyemma_mapping(msm_obj):
    mapping = {}
    for i, j in zip(range(msm_obj.nstates_full), msm_obj.active_set):
        if i == j:
            mapping[j] = i
        else:
            mapping[j] = range(msm_obj.nstates_full)[i]
    return mapping


frame_selectors = {
    'random': get_random_per_bin,
    'specified_totals': get_specified_number_per_bin_random,
    'centers': get_centers,
    'kmeans': get_frames_using_kmeans
}

parser = argparse.ArgumentParser()
parser.add_argument('receptor_name', type=str,
                    help="Name to use to use as top-level directory for the docking run tree.")
parser.add_argument('model', type=str,
                    help='A loos-interpretable model file that will permit reading the trajectories in "traj_paths".')
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
parser.add_argument('--align-resid-list', type=str, default=None,
                    help='If provided, use numbers in file as a list of residue IDs. Concatenate align selection string '
                         'with one selecting these resids.')
parser.add_argument('--assignments', type=str, default=None,
                    help='h5 file with assignments from which MSM was built. Obligatory unless using "centers" selector.')
parser.add_argument('--clear-bonds', action=argparse.BooleanOptionalAction,
                    help='If thrown, remove connectivity information from model. If incorrect bonds present, but bonds'
                         ' are needed, use this in conjunction with "--find-bonds" to first clear old bonds then assign'
                         ' new ones.')
parser.add_argument('--find-bonds', type=floatpair, default=None,
                    help='If no bonds (CONECT records in PDB) are provided in model, find bonds using these cutoffs.')
parser.add_argument('--make-receptor-sel-chain-A', action=argparse.BooleanOptionalAction, default=True,
                    help='If thrown, make all atoms a member of chain "A" when writing PDBs.')
parser.add_argument('--mapping', '-m', type=Path, default=None,
                    help='Use a supplied path to a mapping or MSM to handle eq_probs that have been reduced '
                         'relative to the number of clusters in assignments (4x, with ergodic trimming). If not'
                         ' supplying an enspara MSM, mapping should be a .json of a dict providing the "to_mapped" '
                         'trim mapping.')
parser.add_argument('--number-frames', default=10, type=int,
                    help='Number of frames to select per-bin. If a bin has fewer total assignments than this value, '
                         'an error is thrown.')
parser.add_argument('--subset-selection', type=str, default='all',
                    help='A loos selection string to subset all frames by (for example, if water is present it could be'
                         ' stripped here).')
parser.add_argument('--total-per-bin', default=None, type=str,
                    help='Text file with totals to draw from each MSM bin. Should be one column of totals, '
                         'where the row index corresponds to the bin and the entry is the number of frames to draw.')
parser.add_argument('--write-bin-trajs', action=argparse.BooleanOptionalAction,
                    help='If thrown, write a DCD with the selected frames in each bin directory.')
parser.add_argument('--write-bin-dtraj', type=Path, default=None,
                    help='Write an enspara RaggedArray with each frame index selected to provided path.')
parser.add_argument('--features', type=Path, default=None,
                    help='Supply a path to features that were used for clustering. This will then be used to pick sufficiently different frames from the bin by using kmeans clustering within a bin')

if __name__ == '__main__':
    for i, arg in enumerate(argv):
        print(i, arg)
    args = parser.parse_args()
    try:
        eq_probs = np.load(args.eq_probs)
    except ValueError:
        try:
            eq_probs = np.load(args.eq_probs, allow_pickle=True).item().eq_probs_
        except AttributeError:
            with open(args.eq_probs, 'rb') as f:
                eq_probs = pickle.load(f).pi

    if args.frame_selector == 'centers':
        assignments = np.arange(eq_probs.shape[0])
    else:
        assignments = ra.load(args.assignments)
    model = loos.createSystem(args.model)
    # check to see if it was requested that we remove bond information
    # (if it's foouling up parsing, by prep script, for example)
    if args.clear_bonds:
        model.clearBonds()
    # Assign bonds using cutoffs here.
    if not model.hasBonds():
        print('Bond specifiers not found in model file.')
        if args.find_bonds:
            print('Defining bonds using distance cutoffs; ', args.find_bonds[0], 'angstrom for heavy atoms, and ',
                  args.find_bonds[1], 'angstroms for hydrogens.')
            model = add_bonds_two_cuts(model, *args.find_bonds)
    # AutoDock usually needs this to produce working parameterizations.
    if args.make_receptor_sel_chain_A:
        for atom in model:
            atom.chainId('A')
    align_sel = ''
    if args.align_resid_list:
        alresids = np.genfromtxt(args.align_resid_list).astype(int)
        align_sel += ' || '.join(['resid == {}'.format(i) for i in alresids]) + " &&"

    # if filename is provided read its contents to obtain align string
    try:
        align_sel += " (" + Path(args.align_selection).read_text() + ")"
    except FileNotFoundError:  # if the string is not a file name, interpret it as a loos selection string.
        align_sel += " (" + args.align_selection + ")"

    # get frame counts, however it's prescribed

    if args.total_per_bin:
        frame_counts = np.loadtxt(args.total_per_bin, dtype=int)

    elif args.assignments:
        low_inds = np.where(
            np.bincount(assignments.flatten()) < args.number_frames)
        if len(low_inds[0]) > 0:
            print('The(se) bin(s) have fewer than your requested samples in them:')
            print(low_inds[0])
            print('You requested this many frames per bin:')
            print(args.number_frames)
            print('Exiting.')
            exit(1)

        frame_counts = args.number_frames
    else:  # this is in the event that we have no assignments, and are just pulling centers
        frame_counts = 1

    mapping = None
    if args.mapping:
        if args.mapping.suffix == '.json':
            mapping = json.load(args.mapping.open())
        elif args.mapping.suffix == '.npy':
            mapping = np.load(args.mapping, allow_pickle=True).item().mapping_.to_mapped
        elif args.mapping.suffix == '.pickle':
            msm = np.load(args.mapping, allow_pickle=True)
            mapping = pyemma_mapping(msm)

        else:
            print(args.mapping, 'does not have an extension that implies it is either a pickled msm or a mapping '
                                'object (.json, .npy or .pickle). Unsupported format. Exiting.')
            exit(2)

    if args.features:
        features = ra.load(args.features)
        chosen_frames = frame_selectors[args.frame_selector](
            assignments,
            eq_probs.shape[0],
            features,
            frame_counts,
            mapping
        )
    else:
        chosen_frames = frame_selectors[args.frame_selector](
            assignments,
            eq_probs.shape[0],
            frame_counts,
            mapping
        )

    if len(args.traj_paths) == 1:
        traj_paths = []
        p_trjs = Path(args.traj_paths[0])
        if p_trjs.suffix == '.txt':
            traj_paths = p_trjs.read_text().split()
        elif p_trjs.suffix == '.pickle':
            traj_paths.append(unpickle_resave_centers(p_trjs))
        elif p_trjs.suffix == '.npy':
           traj_paths = np.load(str(p_trjs))
        else:
            traj_paths.append(p_trjs)
    else:
        traj_paths = args.traj_paths
        print(traj_paths)
    print('aligning with the following selection string:')
    print(align_sel)
    out_path = Path(args.receptor_name) / 'receptor'
    full_inds, subset_vec, align_vec = rip_conformations(
        chosen_frames, model, args.subset_selection, align_sel, traj_paths
    )
    align_samples(subset_vec, align_vec)
    write_sampled_frames(subset_vec, full_inds, out_path, args.write_bin_trajs)
    if args.write_bin_dtraj:
        ra.save(args.write_bin_dtraj, chosen_frames)
