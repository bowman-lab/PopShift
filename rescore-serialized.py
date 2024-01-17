from openmm.app import Simulation, PDBFile
from openff.toolkit import Topology
from openmm import unit as u
import openmm as mm
from pathlib import Path
import loos
import pickle

import argparse as ap


def float_to_kcal_mol_angstrom(number):
    return float(number) * u.kilocalorie_per_mole/u.angstrom


def get_setup(serialize_dir: Path, fn: str,
              restraint_range=None,
              restraint_constant=1000*u.kilocalorie_per_mole/u.angstrom):
    system_xml_p = (serialize_dir/(fn+'-sys')).with_suffix('.xml')
    system = mm.XmlSerializer.deserialize(system_xml_p.read_text())
    top_json_p = (serialize_dir/(fn+'-top')).with_suffix('.json')
    top = Topology.from_json(top_json_p.read_text()).to_openmm()
    if restraint_range:
        restraint = CustomExternalForce(
            'k*periodicdistance(x, y, z, x0, y0, z0)^2')
        restraint_ix = complex_sys.addForce(restraint)
        restraint.addGlobalParameter(
            'k', restraint_constant * u.kilocalories_per_mole / u.angstrom)
        restraint.addPerParticleParameter('x0')
        restraint.addPerParticleParameter('y0')
        restraint.addPerParticleParameter('z0')

    # apply positional restraints to all heavy atoms
    for atom in rl_complex_ommt.atoms():
        if atom.element != 'H' and atom.index < receptor_count:
            restraint.addParticle(atom.index)
    integrator = mm.VerletIntegrator(0.001*mm.unit.picosecond)
    simulation = Simulation(top, system, integrator)
    return simulation, top


def get_energy_from_coords(simulation: Simulation,
                           traj_ag: loos.AtomicGroup, minimize=True):
    simulation.context.setPositions(traj_ag.getCoords()*u.angstroms)
    if minimize:
        simulation.minimizeEnergy(tolerance=0.001)
    state = simulation.context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(u.kilocalories_per_mole) * u.kilocalories_per_mole

# sort order changed frim pick_align_frames.py


def get_chosen_inds_filepaths(receptor_dir: Path):
    # we sort this because we are assuming the bin dirs are numbered,
    # and iterdir has no guaranteed order
    bin_p_list = list(sorted(receptor_dir.iterdir(),
                      key=lambda x: int(x.stem)))
    # the 'stem' (file name no extension) should be a string following the pattern
    # '{traj_index}-{frame_index}', so splitting on '-' and turning both of the substrings into ints
    # gives us a list ordered by bin where each element in the list is a traj_index, frame_index tuple
    # corresponding to a particular drawn frame.
    chosen_inds = [
        tuple(map(int, frame_p.stem.split('-'))) for bin_p in bin_p_list
        for frame_p in bin_p.iterdir()
    ]

    full_inds = np.array([(bin_ix, trj_ix, fra_ix) for bin_ix, samples in enumerate(chosen_inds)
                          for trj_ix, fra_ix in samples],
                         dtype=[('bin', int), ('traj', int), ('frame', int)])
    sorted_chosen = np.argsort(full_inds, order=['bin', 'traj', 'frame'])
    return sorted_chosen


def save_conf_pdb(omt: mm.app.Topology, simulation: Simulation, outpre: Path, suffix: str):
    # get positions out of postmin state, get convert from nanometers to angstroms.
    positions = simulation.context.getState(getPositions=True).getPositions()
    outfile = outpre.parent / (outpre.stem + suffix)
    with outfile.open('w') as f:
        PDBFile.writeFile(omt, positions, f)


def sort_by_filename(path_iterator):
    return sorted(path_iterator, key=lambda x: x.stem)


def vtraj_by_filename(traj_path_iterator, atomic_group):
    return pyloos.VirtualTrajectory(
        *map(lambda fp: pyloos.Trajectory(str(fp), atomic_group),
             sort_by_filename(traj_path_iterator))
    )


p = ap.ArgumentParser(formatter_class=ap.ArgumentDefaultsHelpFormatter)
p.add_argument('serialized_params', type=Path,
               help='Path to directory holding parameterized topology and systems. '
               'Expects dir to contain six files with names: {complex,receptor,ligand}-{sys.xml,top.json}.')
p.add_argument('receptor_conf', type=Path,
               help='Path to coordinates of receptor. Ordering of atoms must match tops and systems.')
p.add_argument('receptor_dir', type=Path,
               help='Path to directory containing the conformations to use as receptor conformations (expected to be PDBs).')
p.add_argument('ligand_conf', type=Path,
               help='Path to coordinates of ligand. Order of atoms must match top and sys files.')
p.add_argument('ligand_paths', type=Path,
               help='Pickle file containing coordinates of ligand poses (from extract_scores.py).')

p.add_argument('--minimize', action=ap.BooleanOptionalAction, default=True,
               help='If switched to "no-minimize", will not minimize structure before GB calculation.')
p.add_argument('--outconf-prefix', type=Path, default=None,
               help='If provided, write pdb of each finished structure with this prefix.')
p.add_argument('--restraint-k', type=float_to_kcal_mol_angstrom, default=1000 * u.kilocalorie_per_mole/u.angstrom,
               help='If provided, use restraint constant in place of default for positional restraints.')

args = p.parse_args()

receptor_ag = loos.createSystem(str(args.receptor_conf))
ligand_ag = loos.createSystem(str(args.ligand_conf))
with args.ligand_paths.open('rb') as f:
    ligand_paths = pickle.load(f)

complex_simulation, complex_top = get_setup(args.serialized_params, 'complex')
receptor_simulation, receptor_top = get_setup(
    args.serialized_params, 'receptor')
ligand_simulation, ligand_top = get_setup(args.serialized_params, 'ligand')
print('Loaded OpenMM systems. Getting ready to do energy evaluations', flush=True)
traj_z = zip(receptor_traj, ligand_traj)
while next(traj_z, False):
    # always do this receptor first!
    complex_ag = receptor_ag + ligand_ag
    complex_e = get_energy_from_coords(
        complex_simulation, complex_ag, minimize=args.minimize)
    receptor_e = get_energy_from_coords(
        receptor_simulation, receptor_ag, minimize=args.minimize)
    ligand_e = get_energy_from_coords(
        ligand_simulation, ligand_ag, minimize=args.minimize)

    print('complex', complex_e, 'ligand', ligand_e, 'receptor', receptor_e)
    print('Interaction Energy:', complex_e - (receptor_e + ligand_e))

    if args.outconf_prefix:
        save_conf_pdb(receptor_top, receptor_simulation,
                      args.outconf_prefix, '-receptor.pdb')
        save_conf_pdb(ligand_top, ligand_simulation,
                      args.outconf_prefix, '-ligand.pdb')
        save_conf_pdb(complex_top, complex_simulation,
                      args.outconf_prefix, '-complex.pdb')
