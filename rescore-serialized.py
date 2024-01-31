from openmm.app import Simulation, PDBFile
from openmm import Platform
from openff.toolkit import Topology
from openmm import unit as u
import openmm as mm
from pathlib import Path
import loos
from loos import pyloos
import pickle
from enspara import ra
import argparse as ap

import spyrmsd as sp



def map_via_graph(ap: Path, bp: Path):
    a_obmol = sp.io.load(str(ap))
    a_spymol = sp.io.to_molecule(a_obmol, adjacency=True)
    b_obmol = sp.io.load(str(bp))
    b_spymol = sp.io.to_molecule(b_obmol, adjacency=True)
    a_noh = a_spymol.strip()
    b_noh = b_spymol.strip()

    # Convert molecules to graphs
    G1 = a_noh.to_graph()
    G2 = b_noh.to_graph()

    # Get all the possible graph isomorphisms
    isomorphisms = sp.graph.match_graphs(G1, G2)
    return isomorphisms


def float_to_kcal_mol_angstrom(number):
    return float(number) * u.kilocalorie_per_mole/u.angstrom


def get_setup(serialize_dir: Path, fn: str,
              restraint_range=None,
              restraint_constant=1000*u.kilocalorie_per_mole/u.angstrom):
    platform = Platform.getPlatformByName('CPU')
    system_xml_p = (serialize_dir/(fn+'-sys')).with_suffix('.xml')
    system = mm.XmlSerializer.deserialize(system_xml_p.read_text())
    top_json_p = (serialize_dir/(fn+'-top')).with_suffix('.json')
    top = Topology.from_json(top_json_p.read_text()).to_openmm()
    integrator = mm.VerletIntegrator(0.001*mm.unit.picosecond)
    simulation = Simulation(top, system, integrator, platform=platform)
    return simulation, top


# Need different inputs if we are prepping restraints.
def get_setup_restraints(serialize_dir: Path, fn: str, restraint_inds: list[int],
                         restraint_constant=100*u.kilocalorie_per_mole/u.angstrom):
    platform = Platform.getPlatformByName('CPU')
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
    simulation = Simulation(top, system, integrator, platform=platform)
    return simulation, top, restraint, particle_term_inds, restraint_forcegroup


# Don't need to return coords here because they don't change.
def get_energy_from_coords(simulation: Simulation,
                           traj_ag: loos.AtomicGroup):
    simulation.context.setPositions(traj_ag.getCoords()*u.angstroms)
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(
        u.kilocalories_per_mole) * u.kilocalories_per_mole
    return energy


# Need to return energy and coords, because the coords changed.
def get_minimized_energy(simulation: Simulation,
                         traj_ag: loos.AtomicGroup,
                         tolerance=0.001*u.kilocalories_per_mole):
    simulation.context.setPositions(traj_ag.getCoords()*u.angstroms)
    simulation.minimizeEnergy(tolerance=tolerance)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    energy = state.getPotentialEnergy().value_in_unit(
        u.kilocalories_per_mole) * u.kilocalories_per_mole
    positions = state.getPositions(asNumpy=True)
    return energy, positions


# Need different inputs for restraints; also need to return coords.
def get_restrained_energy_from_coords(simulation: Simulation,
                                      traj_ag: loos.AtomicGroup,
                                      restraint_range, restraint_obj,
                                      restraint_group, particle_term_inds,
                                      tolerance=0.001 * u.kilocalorie/(u.mole * u.angstrom)):
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
    postmin_positions = state.getPositions(
        asNumpy=True).value_in_unit(u.angstroms)
    restraint_state = simulation.context.getState(
        getEnergy=True, groups=restraint_group)
    restraint_e = restraint_state.getPotentialEnergy().value_in_unit(
        u.kilocalories_per_mole) * u.kilocalories_per_mole
    return full_e - restraint_e, postmin_positions


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
p.add_argument('receptor_dir', type=Path,
               help='Path to directory containing the conformations to use as receptor conformations (expected to be PDBs).')
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

param_dir = args.param_dir
receptor_ag = loos.createSystem(str(param_dir/'receptor-top.pdb'))
ligand_ag = loos.createSystem(str(param_dir/'ligand-top.pdb'))
with args.pose_paths.open('rb') as f:
    ligand_paths = pickle.load(f)

ligand_paths = [[pose.with_suffix('.pdb') for pose in state_poses]
                for state_poses in ligand_paths]
# Set up simulations, potentially with restraints.
if args.restrain:
    # get restraind indices for ligand heavies.
    ligand_heavies = loos.selectAtoms(ligand_ag, '!hydrogen')
    lig_rest_inds = [at.index() for at in ligand_heavies]

    # generate restraint indices for receptor heavies.
    receptor_heavies = loos.selectAtoms(receptor_ag, '!hydrogen')
    rec_rest_inds = [at.index() for at in receptor_heavies]

    complex_sim, complex_top, cplx_res, cplx_part_term_inds, cplx_res_fg = get_setup_restraints(
        param_dir, 'complex', restraint_inds=rec_rest_inds, restraint_constant=args.restraint_k)
    receptor_sim, receptor_top, rec_res, rec_part_term_inds, rec_res_fg = get_setup_restraints(
        param_dir, 'receptor', restraint_inds=rec_rest_inds, restraint_constant=args.restraint_k)
    # Because ligand bond/angle geometry will be off, the restraint needs to be much lighter to relax pose.
    # ligand_sim, ligand_top, lig_res, lig_part_term_inds, lig_res_fg = get_setup_restraints(
    #     param_dir, 'ligand' , restraint_inds=lig_rest_inds,
    #     restraint_constant=10*u.kilocalorie_per_mole/u.angstrom)
    # ligand_sim, ligand_top = get_setup(param_dir, 'ligand')
else:
    complex_sim, complex_top = get_setup(param_dir, 'complex')
    receptor_sim, receptor_top = get_setup(param_dir, 'receptor')
ligand_sim, ligand_top = get_setup(param_dir, 'ligand')

print('Loaded OpenMM systems. Getting ready to do energy evaluations', flush=True)
# initialize empty RA with correct shape
scores = []
for i, state_pose_ps in enumerate(ligand_paths):
    ligand_traj = vtraj_by_filename(state_pose_ps, ligand_ag)
    # change the paths to get receptor dir paths, from ligand paths
    receptor_paths = list(args.receptor_dir.joinpath(
        *pose_p.parts[-2:]) for pose_p in state_pose_ps)
    receptor_traj = vtraj_by_filename(receptor_paths, receptor_ag)
    traj_zip = zip(receptor_traj, ligand_traj, receptor_paths)
    # Next will call next on the trajes within the zip object, which will update the atomic group coordinates.
    for _, _, receptor_path in traj_zip:
        # always do this receptor first!
        complex_ag = receptor_ag + ligand_ag
        if args.minimize:
            if args.restrain:
                complex_e, complex_crds = get_restrained_energy_from_coords(
                    complex_sim,
                    complex_ag,
                    rec_rest_inds,
                    cplx_res,
                    cplx_res_fg,
                    cplx_part_term_inds
                )
                receptor_e, receptor_crds = get_restrained_energy_from_coords(
                    receptor_sim,
                    receptor_ag,
                    rec_rest_inds,
                    rec_res,
                    rec_res_fg,
                    rec_part_term_inds
                )
                # get just the ligand coordinates out of the complex
                ligand_ag.setCoords(complex_crds[-len(ligand_ag):])

            else:
                complex_e, complex_crds = get_minimized_energy(
                    complex_sim, complex_ag)
                receptor_e, receptor_crds = get_minimized_energy(
                    receptor_sim, receptor_ag)
                # get just the ligand coordinates out of the complex
                ligand_ag.setCoords(complex_crds[-len(ligand_ag):])
                # use the minimized coords to estimate the ligand alone energy
                ligand_e = get_minimized_energy(ligand_sim, ligand_ag)
        else:
            complex_e = get_energy_from_coords(complex_sim, complex_ag)
            receptor_e = get_energy_from_coords(receptor_sim, receptor_ag)
        ligand_e = get_energy_from_coords(ligand_sim, ligand_ag)
        interaction_e = complex_e - (receptor_e + ligand_e)
        # Save and report the scores.
        scores.append(interaction_e.value_in_unit(u.kilocalories_per_mole))
        rec_rel_path = Path().joinpath(*receptor_path.parts[-2:])
        print(rec_rel_path, 'complex', complex_e, 'ligand', ligand_e,
              'receptor', receptor_e, 'Interaction Energy:', interaction_e, flush=True)

        if args.outconf_prefix:
            outdir = args.outconf_prefix/rec_rel_path.with_suffix('')
            outdir.parent.mkdir(parents=True, exist_ok=True)
            save_conf_pdb(receptor_top, receptor_sim, outdir/'receptor.pdb')
            save_conf_pdb(ligand_top, ligand_sim, outdir/'ligand.pdb')
            save_conf_pdb(complex_top, complex_sim, outdir/'complex.pdb')
lengths = [len(state_ps) for state_ps in ligand_paths]
score_array = ra.RaggedArray(scores, lengths=lengths)
# save the results
ra.save(str(args.out_scores), score_array)
