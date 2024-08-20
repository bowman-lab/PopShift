from openmm.app import Simulation, PDBFile, DCDReporter, StateDataReporter
from openmm import Platform
from openff.toolkit import Topology
from openmm import unit as u
import openmm as mm
from pathlib import Path
import argparse as ap
import numpy as np


def get_ligand_setup(serialize_dir: Path, fn: str, temperature: float, platform_name='CPU',
                     timestep=0.002, friction=2):
    platform = Platform.getPlatformByName(platform_name)
    system_xml_p = (serialize_dir/f'{fn}-hconstrained-sys').with_suffix('.xml')
    system = mm.XmlSerializer.deserialize(system_xml_p.read_text())
    top_json_p = (serialize_dir/(fn+'-top')).with_suffix('.json')
    top = Topology.from_json(top_json_p.read_text()).to_openmm()
    pdb = PDBFile(str(serialize_dir/f'{fn}-top.pdb'))
    integrator = mm.LangevinIntegrator(
        temperature, friction, timestep)
    simulation = Simulation(top, system, integrator, platform=platform)
    simulation.context.setPositions(pdb.positions)
    return simulation, top


def float_to_int(input_str: str):
    return int(float(input_str))

def commas_to_float_tuple(input_str: str):
    return tuple(map(int, input_str.split(',')))


p = ap.ArgumentParser(formatter_class=ap.ArgumentDefaultsHelpFormatter)
p.add_argument('param_dir', type=Path,
               help='Path to directory containing ligand simulation and topology '
               '(as made by parameterize_system.py)')
p.add_argument('length', type=float_to_int,
               help="Number of steps to include in the simulation. "
               "Can be any string recognized by python's 'int' builtin, including strings of the form YeX")
p.add_argument('--outdir', type=Path, default=None,
               help='If provided, write outfiles to this directory. '
               'If not provided, will write to "param_dir".')
p.add_argument('-T', '--temperature', type=float, default=300,
               help='Run at this temperature. Note that the GB parameters are '
               'temperature dependent, so you should use whatever you used '
               'for system parameterization here.')
p.add_argument('--save-rate', type=float_to_int, default=100,
               help='Save a frame and write output every "save-rate" steps.')
p.add_argument('--burn-in', type=float_to_int, default=1e4,
               help='Number of frames to run to burn in the trajectory. '
               'No trajectory is written for this portion of the simulation, '
               'and state data is written to a separate file.')
p.add_argument('--platform-name', type=str, default='CPU',
               help='Platform name for openmm simulation. CPU is recommended because of normal ligand size.')
p.add_argument('--burn-in-timestep', type=float, default=0.001,
               help='burn-in timestep to use, in picoseconds')
p.add_argument('--time-step', '-t', type=float, default=0.002,
               help='Production timestep to use, in picoseconds.')
p.add_argument('--collision-freq', '-g', type=float, default=2.0,
               help='Collision frequency for LangevinIntegrator, in 1/ps.')
p.add_argument('--burn-in-ramp', type=commas_to_float_tuple, default=None,
               help='If provided 3 comma separated values to warm system in stages by. '
               'Add one interval to the final value to finish burn in at that temperature;'
               'The triple will be treaded like the positional arguments to numpy.linspace.')
p.add_argument('--debug-forces', type=Path, default=None,
               help='If a prefix is provided, write post minimized state data to a '
               '.npy of forces and a .pdb showing positions.')
p.add_argument('--ligand-tag', type=str, default='ligand',
               help='Use this tag to find system and topology files within a parameter dir.')

args = p.parse_args()
ligand_tag = args.ligand_tag
sim, top = get_ligand_setup(args.param_dir, ligand_tag, args.temperature, 
                            timestep=args.burn_in_timestep, 
                            friction=args.collision_freq)
print('Minimizing ligand conf read from pdb in ligand parameter directory.', flush=True)
sim.minimizeEnergy()
if args.outdir:
    outdir = args.outdir
else:
    outdir = args.param_dir

outdir.mkdir(exist_ok=True, parents=True)
ligand_out_prefix = outdir/f'{ligand_tag}-sim'
burn_in_fn = str(outdir/f'{ligand_tag}-burnin.out')
traj_fn = str(ligand_out_prefix.with_suffix('.dcd'))
out_fn = str(ligand_out_prefix.with_suffix('.out'))
sim.reporters.append(StateDataReporter(burn_in_fn, args.save_rate, step=True,
                                    potentialEnergy=True, temperature=True,
                                    speed=True, elapsedTime=True,
                                    progress=True, totalSteps=args.burn_in))

if args.debug_forces:
    print('saving debug files with prefix:', args.debug_forces )
    state_postmin = sim.context.getState(getPositions=True, getForces=True)
    forcevals = state_postmin.getForces(asNumpy=True)
    forcemags = np.sqrt((forcevals * forcevals).sum(axis=1))
    np.save(args.debug_forces.with_suffix('.npy'), forcevals)
    mags_filename = args.debug_forces.parent / \
                    (args.debug_forces.stem + '-magnitudes.npy')
    np.save(mags_filename, forcemags)
    with args.debug_forces.with_suffix('.pdb').open('w') as f:
        PDBFile.writeFile(top, state_postmin.getPositions(), file=f)
    
print(f'Running burn-in. Check output in: {burn_in_fn}', flush=True)
if args.burn_in_ramp:
    temp_steps = np.linspace(*args.burn_in_ramp)
    burnin_phase_length = int(args.burn_in / len(temp_steps))
    for temp in temp_steps:
        print('temperature', temp)
        sim.context.getIntegrator().setTemperature(temp)
        sim.step(burnin_phase_length)
else:
    sim.context.setVelocitiesToTemperature(args.temperature)
    sim.step(args.burn_in)
sim_time = args.length * args.time_step * u.picosecond
print(f'Finished burn-in; Starting {sim_time} production run. Saving output and trajectory to:'
      f'\n{out_fn}, {traj_fn}', 
      flush=True)
sim.reporters[-1] = StateDataReporter(out_fn, args.save_rate, step=True,
                                    potentialEnergy=True, temperature=True,
                                    speed=True, elapsedTime=True,
                                    progress=True, totalSteps=args.length)
sim.reporters.append(DCDReporter(traj_fn, args.save_rate))
sim.context.getIntegrator().setStepSize(args.time_step)
sim.context.getIntegrator().setTemperature(args.temperature)
sim.context.setStepCount(0)  # reset steps to zero to elide burn-in from the count.
sim.step(args.length)
print(f'Finished main simulation', flush=True)