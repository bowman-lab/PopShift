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


# takes an SDF fn, returns an openFF molecule that has had hydrogens and stereo added
# Adapted from openff.toolkit.utils.rdkit_wrapper._assign_aromaticity_and_stereo_from_3d
def rdkit_sanitize_and_stereo(coords_fp):
    suppl = Chem.SDMolSupplier(str(coords_fp), removeHs=False)
    pose_rdkmol = next(suppl)  # just grabs the first conf off the supplier.
    pose_rdkmol = Chem.rdmolops.AddHs(pose_rdkmol, addCoords=True)
    Chem.SanitizeMol(
        pose_rdkmol,
        Chem.SANITIZE_ALL ^ Chem.SANITIZE_ADJUSTHS
    )
    Chem.AssignStereochemistryFrom3D(pose_rdkmol)
    Chem.rdmolops.Kekulize(pose_rdkmol, clearAromaticFlags=True)
    Chem.SetAromaticity(pose_rdkmol, Chem.AromaticityModel.AROMATICITY_MDL)
    offmol_w_stereo_and_aro = Molecule.from_rdkit(
        pose_rdkmol, allow_undefined_stereo=True, hydrogens_are_explicit=True
    )
    return offmol_w_stereo_and_aro


# ADAPTED FROM OpenFF-Toolkit's Molecule.from_pdb_and_smiles
def gen_offmol_conf_and_smiles(
    file_path: Path,
    smiles: str,
    allow_undefined_stereo: bool = False,
    name: str = "",
    remove_h = False,
    remove_h_smiles=False
):
    """
    Create a Molecule from a pdb file and a SMILES string using RDKit.

    Requires RDKit to be installed.

    The molecule is created and sanitised based on the SMILES string, we then find a mapping
    between this molecule and one from the PDB based only on atomic number and connections.
    The SMILES molecule is then reindexed to match the PDB, the conformer is attached, and the
    molecule returned.

    Note that any stereochemistry in the molecule is set by the SMILES, and not the coordinates
    of the PDB.

    Parameters
    ----------
    file_path
        PDB or SDF file path--should be a conformer with an atom order it is desired to match.
    smiles
        a valid smiles string for the pdb, used for stereochemistry, formal charges, and bond order
    allow_undefined_stereo
        If false, raises an exception if SMILES contains undefined stereochemistry.
    name
        An optional name for the output molecule.
    remove_h
        Remove hydrogens from the conformer.
    remove_h_smiles
        Remove hydrogens from the smiles used to assign bond orders to the conformer.

    Returns
    --------
    molecule
        An OFFMol instance with ordering the same as used in the PDB file.

    Raises
    ------
    InvalidConformerError
    """
    # Make the molecule from smiles
    offmol = Molecule.from_smiles(
        smiles,
        allow_undefined_stereo=allow_undefined_stereo,
    )
    smi_rdkmol = Chem.MolFromSmiles(smiles)
    if remove_h_smiles:
        smi_rdkmol = Chem.rdmolops.RemoveHs(smi_rdkmol)
    if file_path.suffix == '.pdb':
        conf_rdkmol = Chem.MolFromPDBFile(str(file_path), removeHs=remove_h)
    elif file_path.suffix == '.sdf':
        conf_rdkmol = next(Chem.SDMolSupplier(str(file_path), removeHs=remove_h))
    assigned_rdk_conf = Chem.AssignBondOrdersFromTemplate(smi_rdkmol, conf_rdkmol)
    hydro_rdk_conf = Chem.rdmolops.AddHs(assigned_rdk_conf, addCoords=True)
    conf_mol = Molecule.from_rdkit(
        hydro_rdk_conf,
        allow_undefined_stereo=True,
        hydrogens_are_explicit=True
    )
    # check isomorphic and get the mapping if true the mapping will be
    # dict[offmol_index, pdbmol_index] sorted by offmol index
    isomorphic, mapping = Molecule.are_isomorphic(
        offmol,
        conf_mol,
        return_atom_map=True,
        aromatic_matching=False,
        formal_charge_matching=False,
        bond_order_matching=False,
        atom_stereochemistry_matching=False,
        bond_stereochemistry_matching=False,
    )
    # return offmol, pdbmol, pdb_rdkmol, isomorphic, mapping
    if mapping is None:
        from openff.toolkit.topology.molecule import InvalidConformerError

        raise InvalidConformerError(
            "The PDB and SMILES structures do not match.")

    new_mol = offmol.remap(mapping)

    # the pdb conformer is in the correct order so just attach it here
    new_mol._add_conformer(conf_mol.conformers[0])

    # Take residue info from PDB
    for confatom, newatom in zip(conf_mol.atoms, new_mol.atoms):
        newatom.metadata.update(confatom.metadata)
        newatom.name = confatom.name
    new_mol.add_default_hierarchy_schemes()

    if name:
        new_mol.name = name
    else:
        new_mol.name = file_path.stem
    return new_mol


p = ap.ArgumentParser(formatter_class=ap.ArgumentDefaultsHelpFormatter)
p.add_argument('receptor_pdb', type=Path, 
               help='Receptor PDB file to parameterize.')
p.add_argument('ligand', type=str,
               help='SDF, or Smiles string of ligand to parameterize.')
p.add_argument('out_dir', type=Path,
               help='Name of directory to write parameterized jsons and systems to.')
p.add_argument('--ligand-sdf', '-s', action=ap.BooleanOptionalAction, default=False,
               help='If thrown, interpret ligand as path to SDF file.')
p.add_argument('--ligand-from-conf', '-L', type=Path, default=None,
               help='If provided, use the smiles and the path to a pose to produce '
               'parameters with atom order matching poses from docking.')
p.add_argument('--ligand-ff', type=str, default='openff_unconstrained-2.2.0.offxml',
               help='Name of force field to use as an argument to SMIRNOFFTemplateGenerator.')
p.add_argument('--receptor-ff', type=str, default='amber/protein.ff14SB.xml',
               help='Name of force field xml to use as argument to openmm.ForceField.')
p.add_argument('--implicit', type=str, default='implicit/gbn2.xml', 
               help='Name of implicit solvent model to use as argument to openmm ForceField. '
               'If a path is provided, will read that file instead.')
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
p.add_argument('--name', type=str, default=None,
               help='If provided, will use to name the parameterized ligand' 
               'as part of openff-toolkit.Molecule metadata.')

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
    ligand = rdkit_sanitize_and_stereo(args.ligand_sdf)
elif args.ligand_from_conf:
    ligand = gen_offmol_conf_and_smiles(args.ligand_from_conf, args.ligand, name=args.name)
else:
    ligand = Molecule.from_smiles(args.ligand)
if nagl_model_path:
    ligand.assign_partial_charges(nagl_model_path, toolkit_registry=NAGLToolkitWrapper())
    print('Assigned Charges using NAGL.', flush=True)
else:
    ligand.generate_conformers(n_conformers=1)
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
forcefield = ForceField(args.receptor_ff, args.implicit)
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



