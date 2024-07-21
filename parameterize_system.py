import argparse as ap
from openff.toolkit import Molecule, Topology
try:
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
    from openff.nagl_models import get_models_by_type
    nagl_model_path = get_models_by_type("am1bcc")[-1]
except ModuleNotFoundError:
    nagl_model_path = None
    print('Could not find NAGLToolkitWrapper; is openff-nagl installed?',
          'Falling back to single point AM1BCC for charges.')
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openmm.app import ForceField
from openmm import XmlSerializer, CustomExternalForce
from openmm import unit as u
from pathlib import Path
from math import sqrt
from rdkit.Chem import AllChem as Chem


def off_serialize(outdir: Path, name, openff_obj, sdf=False):
    json_str = openff_obj.to_json()
    outp = (outdir/name).with_suffix('.json')
    outp.write_text(json_str)
    openff_obj.to_file(outp.with_suffix('.pdb'))
    if sdf:
        openff_obj.to_file(outp.with_suffix('.sdf'))
    return outp


def omm_serialize(outdir: Path, name, omm_obj):
    xml_str = XmlSerializer.serialize(omm_obj)
    outp = (outdir/name).with_suffix('.xml')
    outp.write_text(xml_str)
    return outp


p = ap.ArgumentParser(formatter_class=ap.ArgumentDefaultsHelpFormatter)
p.add_argument('receptor_pdb', type=Path, 
               help='Receptor PDB file to parameterize.')
p.add_argument('ligand', type=str,
               help='SDF, or Smiles string of ligand to parameterize.')
p.add_argument('out_dir', type=Path,
               help='Name of directory to write parameterized jsons and systems to.')
p.add_argument('--ligand-sdf', '-s', action=ap.BooleanOptionalAction, default=False,
               help='If thrown, interpret ligand as path to SDF file.')
p.add_argument('--ligand-ff', type=str, default='openff_unconstrained-2.2.0.offxml',
               help='Name of force field to use as an argument to SMIRNOFFTemplateGenerator.')
p.add_argument('--receptor-ff', type=str, default='amber/protein.ff14SB.xml',
               help='Name of force field xml to use as argument to openmm.ForceField.')
p.add_argument('--implicit', type=str, default='gbn2.xml', 
               help='Name of implicit solvent model to use as argument to openmm ForceField. '
               'Will have "implicit/" prepended.')
p.add_argument('--restraint-k', type=float, default=None,
               help='Restrain receptor to the starting conformation with the provided force constant;'
                ' for minimized energy evaluations. Assumes k is in  kcal/(mol * Angstrom).')
p.add_argument('--write-sdf', action=ap.BooleanOptionalAction, default=True,
               help='If cancelled, do not write ligand SDF.')
p.add_argument('--temperature', type=float, default=300,
               help='Temperature for GB parameters, in Kelvin.')
p.add_argument('--solvent-dielectric', type=float, default=78.5,
               help='Dielectric to use for GB parameters.')
p.add_argument('--salt-conc', type=float, default=0.150,
               help='Monovalent salt concentration for ionic strength of GB parameters, in Molar.')
p.add_argument('--kappa', type=float, default=None,
               help='Add screening parameters in as kappa directly, '
               'as opposed to using solvent condition inputs to calculate it.')

args = p.parse_args()

# Compute Kappa for implicit solvent Ionic Strength
temperature = args.temperature
solv_dielectric = args.solvent_dielectric
conc = args.salt_conc
if args.kappa:
    kappa = args.kappa
else:
    kappa = 367.434915*sqrt(conc/(solv_dielectric*temperature))

if not args.out_dir.is_dir():
    args.out_dir.mkdir(parents=True)

if args.ligand_sdf:
    suppl = Chem.SDMolSupplier(args.ligand, removeHs=False)
    some_h_mol = next(suppl)  # just grabs the first conf off the supplier.
    Chem.rdmolops.Kekulize(some_h_mol)
    mol = Chem.rdmolops.AddHs(some_h_mol, addCoords=True)
    ligand = Molecule.from_rdkit(mol)
else:
    ligand = Molecule.from_smiles(args.ligand)
ligand.generate_conformers(n_conformers=1)
if nagl_model_path:
    ligand.assign_partial_charges(nagl_model_path, toolkit_registry=NAGLToolkitWrapper())
else:
    print('Generated ligand conformers from SMILES. Prepping charges.', flush=True)
    ligand.assign_partial_charges(partial_charge_method='am1bcc', 
                                use_conformers=ligand.conformers)
print('Finished getting ligand charges.', flush=True)
receptor = Topology.from_pdb(args.receptor_pdb)
print('Got receptor topology.')
receptor_count = len(list(receptor.atoms))
# make topologies, and serialize them
lig_top = ligand.to_topology()
if args.write_sdf:
    # atom order not preserved by Molecule.to_topology().
    top_ligand = Molecule.from_topology(lig_top)
    top_ligand.to_file(str(args.out_dir.with_suffix('.sdf')), 'sdf')
# always do this receptor first!
rl_complex = receptor + lig_top
# rl_complex = lig_top + receptor 
off_serialize(args.out_dir, 'receptor-top', receptor)
off_serialize(args.out_dir, 'ligand-top', lig_top)
off_serialize(args.out_dir, 'complex-top', rl_complex)

# Create the SMIRNOFF template generator with the default installed force field
smirnoff = SMIRNOFFTemplateGenerator(molecules=ligand, forcefield=args.ligand_ff)
# Create an OpenMM ForceField object with AMBER ff14SB and TIP3P with compatible ions
forcefield = ForceField(args.receptor_ff, 'implicit/' + args.implicit)
# Register the SMIRNOFF template generator
forcefield.registerTemplateGenerator(smirnoff.generator)
print('Getting ready to make force fields.')
# make systems from each of the topologies above
receptor_sys = forcefield.createSystem(receptor.to_openmm(), implicitSolventKappa=kappa)
ligand_sys = forcefield.createSystem(lig_top.to_openmm(), implicitSolventKappa=kappa)
rl_complex_ommt = rl_complex.to_openmm()
complex_sys = forcefield.createSystem(rl_complex_ommt, implicitSolventKappa=kappa)
# optionally ad receptor restraints
if args.restraint_k:
    restraint = CustomExternalForce('k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)')
    restraint_ix = complex_sys.addForce(restraint)
    restraint.addGlobalParameter('k', args.restraint_k * u.kilocalories_per_mole / u.angstrom)
    restraint.addPerParticleParameter('x0')
    restraint.addPerParticleParameter('y0')
    restraint.addPerParticleParameter('z0')
    print('Made positional restraints.') 
    # apply positional restraints to all heavy atoms
    for atom in rl_complex_ommt.atoms():
        if atom.element != 'H' and atom.index < receptor_count:
            restraint.addParticle(atom.index, )
# serialize them into the same output directory
omm_serialize(args.out_dir, 'complex-sys', complex_sys)
omm_serialize(args.out_dir, 'receptor-sys', receptor_sys)
omm_serialize(args.out_dir, 'ligand-sys', ligand_sys)



