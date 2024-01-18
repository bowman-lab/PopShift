from openmm.app import Simulation, PDBFile
from openff.toolkit import Topology
from openmm import unit as u
import openmm as mm
from pathlib import Path
import loos
import pickle
from enspara import ra
import numpy as np
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
    integrator = mm.VerletIntegrator(0.001*mm.unit.picosecond)
    simulation = Simulation(top, system, integrator)
    return simulation, top


def get_setup_restraints(serialize_dir: Path, fn: str, restraint_inds: list[int],
                         restraint_constant=100*u.kilocalorie_per_mole/u.angstrom):
    system_xml_p = (serialize_dir/(fn+'-sys')).with_suffix('.xml')
    system = mm.XmlSerializer.deserialize(system_xml_p.read_text())
    top_json_p = (serialize_dir/(fn+'-top')).with_suffix('.json')
    top = Topology.from_json(top_json_p.read_text()).to_openmm()
    integrator = mm.VerletIntegrator(0.001*mm.unit.picosecond)
    if restraint_inds:
        # non periodic harmonic distance restraint.
        restraint = mm.CustomExternalForce(
            'k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)')
        restraint_forcegroup = system.addForce(restraint)
        restraint.addGlobalParameter('k', restraint_constant)
        restraint.addPerParticleParameter('x0')
        restraint.addPerParticleParameter('y0')
        restraint.addPerParticleParameter('z0')
        particle_term_inds = []
        for atom_ix in restraint_inds:
            # put in placeholder coords since we'll have to overwrite repeadedly later.
            particle_term_inds.append(
                restraint.addParticle(atom_ix, [0.0, 0.0, 0.0]))
    simulation = Simulation(top, system, integrator)
    return simulation, top, restraint, particle_term_inds, restraint_forcegroup


def get_energy_from_coords(simulation: Simulation,
                           traj_ag: loos.AtomicGroup, minimize=True):
    simulation.context.setPositions(traj_ag.getCoords()*u.angstroms)
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(
        u.kilocalories_per_mole) * u.kilocalories_per_mole
    return energy


def get_restrained_energy_from_coords(simulation: Simulation,
                                      traj_ag: loos.AtomicGroup,
                                      restraint_range, restraint_obj,
                                      restraint_group, particle_term_inds,
                                      tolerance=0.001 * u.kilocalories_per_mole):
    positions_angstroms = traj_ag.getCoords()*u.angstroms
    simulation.context.setPositions(positions_angstroms)
    for atom_ix, particle_term_ix in zip(restraint_range, particle_term_inds):
        restraint_obj.setParticleParameters(
            particle_term_ix, atom_ix, positions_angstroms[atom_ix])
        restraint_obj.updateParametersInContext(simulation.context)
    simulation.minimizeEnergy(tolerance=tolerance)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    full_e = state.getPotentialEnergy().value_in_unit(
        u.kilocalories_per_mole) * u.kilocalories_per_mole
    positions = state.getPositions(asNumpy=True)
    restraint_state = simulation.context.getState(
        getEnergy=True, groups=restraint_group)
    restraint_e = restraint_state.getPotentialEnergy().value_in_unit(
        u.kilocalories_per_mole) * u.kilocalories_per_mole
    return full_e - restraint_e, positions


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
p.add_argument('param_dir', type=Path,
               help='Path to directory holding parameterized topology and systems. '
               'Expects dir to contain six files with names: {complex,receptor,ligand}-{sys.xml,top.json}.')
p.add_argument('receptor_conf', type=Path,
               help='Path to coordinates of receptor. Ordering of atoms must match tops and systems.')
p.add_argument('receptor_dir', type=Path,
               help='Path to directory containing the conformations to use as receptor conformations (expected to be PDBs).')
p.add_argument('ligand_conf', type=Path,
               help='Path to coordinates of ligand. Order of atoms must match top and sys files.')
p.add_argument('pose_paths', type=Path,
               help='Pickle file containing coordinates of ligand poses (from extract_scores.py).')
p.add_argument('out_scores', type=Path,
               help="h5 file with enspara RA containing scores in the same order as the docking scores extracted with popshift.")
p.add_argument('--minimize', action=ap.BooleanOptionalAction, default=True,
               help='Minimize before calculating GB energy.')
p.add_argument('--restrain', action=ap.BooleanOptionalAction, default=True,
               help='Restrain receptor heavy atoms during minimization. Ignored if "--no-minimize" is thrown.')
p.add_argument('--outconf-prefix', type=Path, default=None,
               help='If provided, write pdb of each finished structure with this prefix.')
p.add_argument('--restraint-k', type=float_to_kcal_mol_angstrom, default=100 * u.kilocalorie_per_mole/u.angstrom,
               help='If provided, use restraint constant in place of default for positional restraints.')

args = p.parse_args()

receptor_ag = loos.createSystem(str(args.receptor_conf))
ligand_ag = loos.createSystem(str(args.ligand_conf))
with args.pose_paths.open('rb') as f:
    ligand_paths = pickle.load(f)

param_dir = args.param_dir
# prepare to deal with restraints
if args.restrain:
    # generate restraint indices for receptor heavies.
    receptor_heavies = loos.selectAtoms(receptor_ag, '!hydrogen')
    rinds = [at.index() for at in receptor_heavies]
    complex_simulation, complex_top, cplx_res, cplx_part_term_inds, cplx_res_fg = get_setup_restraints(
        param_dir, 'complex', restraint_inds=rinds)
    receptor_simulation, receptor_top, rec_res, rec_part_term_inds, rec_res_fg = get_setup_restraints(
        param_dir, 'receptor', restraint_inds=rinds)
else:
    complex_sim, complex_top = get_setup(param_dir, 'complex')
    receptor_sim, receptor_top = get_setup(param_dir, 'receptor')


ligand_sim, ligand_top = get_setup(param_dir, 'ligand')
print('Loaded OpenMM systems. Getting ready to do energy evaluations', flush=True)
# initialize empty RA with correct shape
lengths = [len(state_ps) for state_ps in ligand_paths]
total = sum(lengths)
score_array = ra.RaggedArray(np.zeros(total), lengths=lengths)
for i, ligand_poses in enumerate(ligand_paths):
    ligand_traj = vtraj_by_filename(ligand_poses, ligand_ag)
    # change the paths to get receptor dir paths, from ligand paths
    receptor_paths = (args.receptor_dir.joinpath(
        *pose[-2:]) for pose in ligand_poses)
    receptor_traj = vtraj_by_filename(receptor_paths, receptor_ag)
    traj_zip = zip(receptor_traj, ligand_traj, receptor_paths)
    # Next will call next on the trajes within the zip object, which will update the atomic group coordinates.
    scores = []
    for _, _, receptor_path in traj_zip:
        # always do this receptor first!
        complex_ag = receptor_ag + ligand_ag
        if args.minimize:
            complex_e, complex_crds = get_restrained_energy_from_coords(
                complex_sim,
                complex_ag,
                rinds,
                cplx_res,
                cplx_res_fg,
                cplx_part_term_inds
            )
            # get just the ligand coordinates out of the complex
            ligand_ag.setCoords(complex_coords[-len(ligand_ag):])
            # use the minimized coords to estimate the ligand alone energy
            ligand_e = get_energy_from_coords(ligand_sim, ligand_ag)
        else:
            complex_e = get_energy_from_coords(complex_sim, complex_ag)
        
        receptor_e = get_energy_from_coords(receptor_sim, receptor_ag)

        interaction_e = complex_e - (receptor_e + ligand_e)
        scores.append(interaction_e)
        print(Path.joinpath(receptor_path.parts[-2:]), 'complex', complex_e, 'ligand', ligand_e,
              'receptor', receptor_e, 'Interaction Energy:', interaction_e)

        if args.outconf_prefix:
            save_conf_pdb(receptor_top, receptor_sim,
                          args.outconf_prefix, '-receptor.pdb')
            save_conf_pdb(ligand_top, ligand_sim,
                          args.outconf_prefix, '-ligand.pdb')
            save_conf_pdb(complex_top, complex_sim,
                          args.outconf_prefix, '-complex.pdb')
    score_array[i] = np.array(scores)

# save the results
ra.save(args.out, score_array)
