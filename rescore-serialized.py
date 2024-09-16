from openmm.app import Simulation, PDBFile
from openmm import Platform
from openff.toolkit import Topology
from openmm import unit as u
import openmm as mm
import numpy as np
from pathlib import Path
import loos
from loos import pyloos
import pickle
from enspara import ra
import argparse as ap
from rdkit.Chem import AllChem as Chem


# takes a hydrogenated and correct template and a target molecule, returns target's coords
def add_hs_get_coords(template: Chem.Mol, mol: Chem.Mol):
    matched = Chem.AssignBondOrdersFromTemplate(template, mol)
    molh = Chem.AddHs(matched, addCoords=True)
    return molh.GetConformer().GetPositions()


def get_mol_sdf(sdf: Path, remove_hs=False):
    return next(Chem.SDMolSupplier(str(sdf), removeHs=remove_hs))


def get_multiposes_sdf(sdf: Path, remove_hs=False):
    return Chem.SDMolSupplier(str(sdf), removeHs=remove_hs)


def float_to_kcal_mol_angstrom(number):
    return float(number) * u.kilocalorie_per_mole/u.angstrom


# specifically for the ligand, since we'll also need an openforcefield molecule
def get_ligand_setup(serialize_dir: Path, fn: str):
    platform = Platform.getPlatformByName('CPU')
    system_xml_p = (serialize_dir/(fn+'-sys')).with_suffix('.xml')
    system = mm.XmlSerializer.deserialize(system_xml_p.read_text())
    top_json_p = (serialize_dir/(fn+'-top')).with_suffix('.json')
    top = Topology.from_json(top_json_p.read_text())
    ommtop = top.to_openmm()
    # generate an RDKit molecule, for reading conformations
    mol = next(top.molecules).to_rdkit()
    integrator = mm.VerletIntegrator(0.001*mm.unit.picosecond)
    simulation = Simulation(top, system, integrator, platform=platform)
    return simulation, ommtop, mol


# Read serialized files, build an openmm simulation and topology
def get_setup(serialize_dir: Path, fn: str):
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


# expects coordinates to be dimensionless floats that correspond to Angstroms,
# or dimensioned quantities openMM knows how to convert to Angstroms with multiplication.
# Don't need to return coords here because they don't change.
def get_energy_from_coords(simulation: Simulation,
                           coords):
    simulation.context.setPositions(coords)
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(
        u.kilocalories_per_mole) * u.kilocalories_per_mole
    return energy


# expects coordinates to be dimensionless floats that correspond to Angstroms,
# or dimensioned quantities openMM knows how to convert to Angstroms with multiplication.
# Need to return energy and coords, because the coords changed.
def get_minimized_energy(simulation: Simulation,
                         coords,
                         tolerance=0.001*u.kilocalories_per_mole):
    simulation.context.setPositions(coords * u.angstroms)
    simulation.minimizeEnergy(tolerance=tolerance)
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(
        u.kilocalories_per_mole) * u.kilocalories_per_mole
    return energy


# expects coordinates to be dimensionless floats that correspond to Angstroms,
# or dimensioned quantities openMM knows how to convert to Angstroms with multiplication.
# Need to return energy and coords, because the coords changed.
def get_min_energy_coords(simulation: Simulation,
                          coords,
                          tolerance=0.001*u.kilocalories_per_mole):
    simulation.context.setPositions(coords * u.angstroms)
    simulation.minimizeEnergy(tolerance=tolerance)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    post_min_coords = state.getPositions(asNumpy=True)
    energy = state.getPotentialEnergy().value_in_unit(
        u.kilocalories_per_mole) * u.kilocalories_per_mole
    return energy, post_min_coords


# expects coordinates to be dimensionless floats that correspond to Angstroms,
# or dimensioned quantities openMM knows how to convert to Angstroms with multiplication.
# Need different inputs for restraints; also need to return coords.
def get_restrain_min_energy(simulation: Simulation,
                            coords,
                            restraint_range, restraint_obj,
                            restraint_group, particle_term_inds,
                            tolerance=0.001 * u.kilocalorie/(u.mole * u.angstrom)):
    positions_angstroms = coords * u.angstroms
    simulation.context.setPositions(positions_angstroms)
    for atom_ix, particle_term_ix in zip(restraint_range, particle_term_inds):
        restraint_obj.setParticleParameters(
            particle_term_ix, atom_ix, positions_angstroms[atom_ix])
        restraint_obj.updateParametersInContext(simulation.context)
    simulation.minimizeEnergy(tolerance=tolerance)
    state = simulation.context.getState(getEnergy=True)
    full_e = state.getPotentialEnergy().value_in_unit(
        u.kilocalories_per_mole) * u.kilocalories_per_mole
    restraint_state = simulation.context.getState(
        getEnergy=True, groups=restraint_group)
    restraint_e = restraint_state.getPotentialEnergy().value_in_unit(
        u.kilocalories_per_mole) * u.kilocalories_per_mole
    return full_e - restraint_e


# expects coordinates to be dimensionless floats that correspond to Angstroms,
# or dimensioned quantities openMM knows how to convert to Angstroms with multiplication.
# Need different inputs for restraints; also need to return coords.
def get_restrain_min_energy_coords(simulation: Simulation,
                                   coords,
                                   restraint_range, restraint_obj,
                                   restraint_group, particle_term_inds,
                                   tolerance=0.001 * u.kilocalorie/(u.mole * u.angstrom)):
    positions_angstroms = coords * u.angstroms
    simulation.context.setPositions(positions_angstroms)
    for atom_ix, particle_term_ix in zip(restraint_range, particle_term_inds):
        restraint_obj.setParticleParameters(
            particle_term_ix, atom_ix, positions_angstroms[atom_ix])
        restraint_obj.updateParametersInContext(simulation.context)
    simulation.minimizeEnergy(tolerance=tolerance)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    post_min_coords = state.getPositions(asNumpy=True)
    full_e = state.getPotentialEnergy().value_in_unit(
        u.kilocalories_per_mole) * u.kilocalories_per_mole
    restraint_state = simulation.context.getState(
        getEnergy=True, groups=restraint_group)
    restraint_e = restraint_state.getPotentialEnergy().value_in_unit(
        u.kilocalories_per_mole) * u.kilocalories_per_mole
    return full_e - restraint_e, post_min_coords


# Openmm topology and simulation saved to PDB using openmm utilities.
# Coordinates come from current context in simulation.
def save_conf_pdb(omt: mm.app.Topology, simulation: Simulation, outpath: Path):
    # get positions out of postmin state, get convert from nanometers to angstroms.
    positions = simulation.context.getState(getPositions=True).getPositions()
    with outpath.open('w') as f:
        PDBFile.writeFile(omt, positions, f)


def sort_by_filename(path_iterator):
    return sorted(path_iterator, key=lambda x: x.stem)


def vtraj_by_filename(traj_path_iterator, atomic_group):
    return pyloos.VirtualTrajectory(
        *map(lambda fp: pyloos.Trajectory(str(fp), atomic_group),
             sort_by_filename(traj_path_iterator))
    )


"""
The following are functors that use the openMM simulation objects 
this scrip preps to compute and report energies for complex, 
receptor, and ligand components, for calculation of the 
interaction energy.

Because this can either be done with no minimization, minimization,
or minimization with restraints on the receptor, InterEnergy is 
subclassed to give a uniform interface for the other two modes of 
calculation. If other approaches are desired, the hope is that they
could also be implemented by subclassing InterEnergy.
"""


class InterEnergy:
    def __init__(self, ligand_sim, receptor_sim, complex_sim):
        self.ligand_sim = ligand_sim
        self.receptor_sim = receptor_sim
        self.complex_sim = complex_sim

    def __call__(self, frame_coords, pose_coords):
        # always do this receptor first!
        posed_complex_coords = np.concatenate((frame_coords, pose_coords))
        complex_e = get_energy_from_coords(
            self.complex_sim, posed_complex_coords)
        receptor_e = get_energy_from_coords(self.receptor_sim, frame_coords)
        ligand_e = get_energy_from_coords(self.ligand_sim, pose_coords)
        return complex_e, receptor_e, ligand_e


class MinimizedInterEnergy(InterEnergy):
    def __init__(self, ligand_sim, receptor_sim, complex_sim):
        super().__init__(ligand_sim, receptor_sim, complex_sim)

    def __call__(self, frame_coords, pose_coords):
        # always do this receptor first!
        posed_complex_coords = np.concatenate((frame_coords, pose_coords))
        complex_e, min_complex_coords = get_min_energy_coords(
            self.complex_sim, posed_complex_coords)
        receptor_e = get_minimized_energy(self.receptor_sim, frame_coords)
        min_lig_coords = min_complex_coords[-len(pose_coords):]
        ligand_e = get_energy_from_coords(self.ligand_sim, min_lig_coords)
        return complex_e, receptor_e, ligand_e


class RestrainedInterEnergy(InterEnergy):
    def __init__(self, ligand_sim, receptor_sim, complex_sim,
                 rec_rest_inds, rec_res, rec_res_fg, rec_part_term_inds,
                 cplx_res, cplx_res_fg, cplx_part_term_inds):
        super().__init__(ligand_sim, receptor_sim, complex_sim)
        self.rec_rest_inds = rec_rest_inds
        self.rec_res = rec_res
        self.rec_res_fg = rec_res_fg
        self.rec_part_term_inds = rec_part_term_inds
        self.cplx_res = cplx_res
        self.cplx_res_fg = cplx_res_fg
        self.cplx_part_term_inds = cplx_part_term_inds

    def __call__(self, frame_coords, pose_coords):
        # always do this receptor first!
        posed_complex_coords = np.concatenate((frame_coords, pose_coords))
        complex_e, min_complex_coords = get_restrain_min_energy_coords(
            self.complex_sim,
            posed_complex_coords,
            self.rec_rest_inds,
            self.cplx_res,
            self.cplx_res_fg,
            self.cplx_part_term_inds
        )
        receptor_e = get_restrain_min_energy(
            self.receptor_sim,
            frame_coords,
            self.rec_rest_inds,
            self.rec_res,
            self.rec_res_fg,
            self.rec_part_term_inds
        )
        min_lig_coords = min_complex_coords[-len(pose_coords):]
        ligand_e = get_energy_from_coords(self.ligand_sim, min_lig_coords)
        return complex_e, receptor_e, ligand_e


# Using the simulation objects contained in an interaction energy calculator,
# write structures for the three component systems. `pose_index` indicates
# which pose in a multi-pose file is being written, and defaults to zero if
# multi-pose analysis is not requested.
def write_pdbs_from_calculator(calculator, out_prefix, receptor_relative_path,
                               receptor_top, ligand_top, complex_top, pose_index=0):
    outdir = out_prefix / (receptor_relative_path.parent /
                           receptor_relative_path.stem)
    outdir.mkdir(parents=True, exist_ok=True)
    save_conf_pdb(receptor_top, calculator.receptor_sim,
                  outdir / f'receptor-{pose_index:03}.pdb')
    save_conf_pdb(ligand_top, calculator.ligand_sim,
                  outdir / f'ligand-{pose_index:03}.pdb')
    save_conf_pdb(complex_top, calculator.complex_sim,
                  outdir / f'complex-{pose_index:03}.pdb')


p = ap.ArgumentParser(formatter_class=ap.ArgumentDefaultsHelpFormatter)
p.add_argument('param_dir', type=Path,
               help='Path to directory holding parameterized topology and systems. '
               'Expects dir to contain six files with names: {complex,receptor,ligand}-{sys.xml,top.json}.')
p.add_argument('receptor_dir', type=Path,
               help='Path to directory containing the conformations to use as receptor conformations (expected to be PDBs).')
p.add_argument('pose_paths', type=Path,
               help='File containing coordinates of ligand poses (from extract_scores.py). '
               'Assumes file is a pickle unless extension is .txt, in which case it assumes text.')
p.add_argument('out_scores', type=Path,
               help="h5 file with enspara RA containing scores in the same order as the docking scores extracted with popshift.")
p.add_argument('--minimize', action=ap.BooleanOptionalAction, default=True,
               help='Minimize before calculating GB energy.')
p.add_argument('--rel-to', '-d', type=Path, default=None,
               help='Directory to look within for ligand pose paths. If none, then will look at parent of receptor dir.')
p.add_argument('--restrain', action=ap.BooleanOptionalAction, default=True,
               help='Restrain receptor heavy atoms during minimization. Ignored if "--no-minimize" is thrown.')
p.add_argument('--outconf-prefix', type=Path, default=None,
               help='If provided, write pdb of each finished structure with this prefix.')
p.add_argument('--restraint-k', type=float_to_kcal_mol_angstrom, default=100 * u.kilocalorie_per_mole/u.angstrom,
               help='If provided, use restraint constant in place of default for positional restraints.')
p.add_argument('--add-hydrogens', action=ap.BooleanOptionalAction, default=False,
               help='If thrown, assumes poses are stored in SDFs that may be missing all or some hydrogens. '
               'Tries to add hydrogens back to each pose using the RDKit.')
p.add_argument('--multi-pose', action=ap.BooleanOptionalAction, default=False,
               help='If thrown, interpret each pose file as a multi-conformer file, and rescore each pose. '
               'Otherwise, just rescore the first conformation. Only implemented for SDFs at present.')

args = p.parse_args()

param_dir = args.param_dir
receptor_ag = loos.createSystem(str(param_dir/'receptor-top.pdb'))
ligand_ag = loos.createSystem(str(param_dir/'ligand-top.pdb'))
if args.rel_to:
    top_dir = args.rel_to
else:
    top_dir = args.receptor_dir.parent
if args.pose_paths.suffix == '.txt':
    ligand_pl = args.pose_paths.read_text().strip().split()
    ligand_paths = [[top_dir / ligand_path] for ligand_path in ligand_pl]
else:
    with args.pose_paths.open('rb') as f:
        ligand_pll = pickle.load(f)
    ligand_paths = [[top_dir / pose for pose in state_poses]
                    for state_poses in ligand_pll]
# Set up simulations, potentially with restraints.
ligand_sim, ligand_top, ligand_rdkit_mol = get_ligand_setup(
    param_dir, 'ligand')

if args.minimize:
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

        ie_calculator = RestrainedInterEnergy(ligand_sim, receptor_sim, complex_sim,
                                              rec_rest_inds, rec_res, rec_res_fg, rec_part_term_inds,
                                              cplx_res, cplx_res_fg, cplx_part_term_inds)
    else:
        complex_sim, complex_top = get_setup(param_dir, 'complex')
        receptor_sim, receptor_top = get_setup(param_dir, 'receptor')
        ie_calculator = MinimizedInterEnergy(
            ligand_sim, receptor_sim, complex_sim)
else:
    complex_sim, complex_top = get_setup(param_dir, 'complex')
    receptor_sim, receptor_top = get_setup(param_dir, 'receptor')
    ie_calculator = InterEnergy(ligand_sim, receptor_sim, complex_sim)

if ligand_paths[0][0].suffix == '.sdf':
    print('Reading SDFs with the RDKit.')
    do_ligand_updates_ag = False
else:
    print('Reading ligand PDBs with LOOS.')
    do_ligand_updates_ag = True

if args.add_hydrogens:
    if do_ligand_updates_ag:
        raise ValueError('This script can only add Hydrogens on the fly to ligand SDFs,'
                         ' but you asked to add hydrogens and supplied files with '
                         f'the extension: "{args.ligand_paths[0][0].suffix}"')
    if args.multi_pose:
        def get_pose_iter(x): return get_multiposes_sdf(x, remove_hs=True)
    else:
        def get_pose(x): return get_mol_sdf(x, remove_hs=True)
    ligand_rdkit_mol_noh = Chem.rdmolops.RemoveHs(ligand_rdkit_mol)

    def get_pose_coords(pose_mol):
        return add_hs_get_coords(ligand_rdkit_mol_noh, pose_mol)
else:
    if do_ligand_updates_ag:
        if args.multi_pose:
            def get_pose_iter(pose_p): return pyloos.Trajectory(
                str(pose_p), ligand_ag)
        else:
            def get_pose(ligand_traj): return next(ligand_traj)

        def get_pose_coords(ag): return ag.getCoords()
    else:
        if args.multi_pose:
            get_pose_iter = get_multiposes_sdf
        else:
            get_pose = get_mol_sdf

        def get_pose_coords(mol): return mol.GetConformer().GetPositions()


print('Loaded OpenMM systems. Getting ready to do energy evaluations', flush=True)
# initialize empty lists to retain scores, and track lengths.
scores = []
lengths = []
# loop over ligand and receptor poses, get coords, do energy calx, optionally save poses.
for i, state_pose_ps in enumerate(ligand_paths):
    print('Loaded ligand paths for state', i, flush=True)
    lengths.append(len(state_pose_ps))
    # change the paths to get receptor dir paths from ligand paths
    receptor_paths = list(args.receptor_dir.joinpath(
        *pose_p.parts[-2:]).with_suffix('.pdb')
        for pose_p in state_pose_ps)
    receptor_traj = vtraj_by_filename(receptor_paths, receptor_ag)
    print('Loaded receptor paths for state', i, flush=True)
    if do_ligand_updates_ag and not args.multi_pose:
        ligand_traj = vtraj_by_filename(state_pose_ps, ligand_ag)
    traj_zip = zip(receptor_traj, state_pose_ps, receptor_paths)
    # for-loop will call next on the trajes within the zip object,
    # which will update the atomic group coordinates.
    for _, state_pose, receptor_path in traj_zip:
        # If docking failed, there will be an empty file. 
        # Check for this, and proceed by recording zero.
        if state_pose.stat().st_size == 0:
            print(state_pose, 'is empty; recording 0.0.')
            if args.multi_pose:
                scores.append(np.array([0]))
            else:
                scores.append(0)
        # proceed with getting the coordinates and doing the calculation
        else:
            frame_coords = receptor_ag.getCoords()
            if args.multi_pose:
                pose_iter = get_pose_iter(state_pose)
                pose_scores = []
                rec_rel_path = Path().joinpath(*receptor_path.parts[-2:])
                for pose_index, pose in enumerate(pose_iter):
                    pose_coords = get_pose_coords(pose)
                    complex_e, receptor_e, ligand_e = ie_calculator(
                        frame_coords, pose_coords)
                    interaction_e = complex_e - (receptor_e + ligand_e)
                    # Save and report the scores.
                    pose_scores.append(
                        interaction_e.value_in_unit(u.kilocalories_per_mole))
                    print(rec_rel_path, 'pose-index', pose_index, 'complex', complex_e, 'ligand', ligand_e,
                        'receptor', receptor_e, 'Interaction Energy:', interaction_e, flush=True)
                    if args.outconf_prefix:
                        write_pdbs_from_calculator(ie_calculator, args.outconf_prefix, rec_rel_path,
                                                ligand_top, complex_top, pose_index=pose_index)
                scores.append(np.array(pose_scores))
            else:
                pose = get_pose(state_pose)
                pose_coords = get_pose_coords(pose)
                print(pose.GetNumAtoms(), len(ligand_ag),
                    len(pose_coords), len(receptor_ag))
                complex_e, receptor_e, ligand_e = ie_calculator(
                    frame_coords, pose_coords)
                interaction_e = complex_e - (receptor_e + ligand_e)
                scores.append(interaction_e.value_in_unit(u.kilocalories_per_mole))
                rec_rel_path = Path().joinpath(*receptor_path.parts[-2:])
                print(rec_rel_path, 'complex', complex_e, 'ligand', ligand_e,
                    'receptor', receptor_e, 'Interaction Energy:', interaction_e, flush=True)
                if args.outconf_prefix:
                    write_pdbs_from_calculator(ie_calculator, args.outconf_prefix, rec_rel_path, receptor_top,
                                            ligand_top, complex_top)


score_array = ra.RaggedArray(scores, lengths=lengths)
# save the results
ra.save(str(args.out_scores), score_array)
