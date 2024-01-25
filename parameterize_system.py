import argparse as ap
from openff.toolkit import Molecule, Topology
from openff.toolkit import ForceField as off_ForceField
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openmm.app import ForceField
from openmm import XmlSerializer, CustomExternalForce
from openmm import unit as u
from pathlib import Path
from math import sqrt


def off_serialize(outdir: Path, name, openff_obj):
    json_str = openff_obj.to_json()
    outp = (outdir/name).with_suffix('.json')
    outp.write_text(json_str)
    openff_obj.to_file(outp.with_suffix('.pdb'))
    return outp


def omm_serialize(outdir: Path, name, omm_obj):
    xml_str = XmlSerializer.serialize(omm_obj)
    outp = (outdir/name).with_suffix('.xml')
    outp.write_text(xml_str)
    return outp


# excerpted from 
# https://github.com/jeff231li/HG-GBSA-Optimization/blob/main/paper_simulations/01-optimization/result/optimize/openff-2.0.0-GBSA_OBC2-tagged.offxml
gb_snippet_setiadi = \
'''<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <Author>The Open Force Field Initiative AND J. Chodera, (MSKCC); J. Wagner (Open Force Field)</Author>
    <Date>2021-08-16 AND 2019-08-08</Date>
    <GBSA version="0.3" gb_model="OBC2" solvent_dielectric="78.5" solute_dielectric="1" sa_model="ACE" surface_area_penalty="5.4 * angstrom**-2 * mole**-1 * calorie" solvent_radius="1.4 * angstrom">
        <Atom smirks="[*:1]" radius="0.15 * nanometer" scale="0.8"/>
        <Atom smirks="[#1:1]" radius="7.068228303295e-02 * nanometer" scale="0.85" parameterize="radius"/>
        <Atom smirks="[#1:1]~[#7]" radius="1.218256113149e-01 * nanometer" scale="0.85" parameterize="radius"/>
        <Atom smirks="[#6:1]" radius="1.662784203155e-01 * nanometer" scale="0.72" parameterize="radius"/>
        <Atom smirks="[#7:1]" radius="5.327800106182e-02 * nanometer" scale="0.79" parameterize="radius"/>
        <Atom smirks="[#8:1]" radius="1.632214397109e-01 * nanometer" scale="0.85" parameterize="radius"/>
        <Atom smirks="[#9:1]" radius="0.15 * nanometer" scale="0.88"/>
        <Atom smirks="[#14:1]" radius="0.21 * nanometer" scale="0.8"/>
        <Atom smirks="[#15:1]" radius="0.185 * nanometer" scale="0.86"/>
        <Atom smirks="[#16:1]" radius="0.18 * nanometer" scale="0.96"/>
        <Atom smirks="[#17:1]" radius="0.17 * nanometer" scale="0.8"/>
    </GBSA>
</SMIRNOFF>'''


p = ap.ArgumentParser(formatter_class=ap.ArgumentDefaultsHelpFormatter)
p.add_argument('receptor_pdb', type=Path, 
               help='Receptor PDB file to parameterize.')
p.add_argument('ligand_smiles', type=str,
               help='Smiles string of ligand to parameterize.')
p.add_argument('out_dir', type=Path,
               help='Name of directory to write parameterized jsons and systems to.')
p.add_argument('--ligand-ff', type=str, default='openff_unconstrained-2.1.0.offxml',
               help='Name of force field to use as an argument to SMIRNOFFTemplateGenerator.')
p.add_argument('--ligand-implicit', type=str, default=gb_snippet_setiadi,
               help='GB parameters to use for atoms parameterized by openff.')
p.add_argument('--receptor-ff', type=str, default='amber/protein.ff14SB.xml',
               help='Name of force field xml to use as argument to openmm.ForceField.')
p.add_argument('--implicit', type=str, default='implicit/obc2.xml', 
               help='Name of implicit solvent model to use as argument to openmm ForceField.')
p.add_argument('--restraint-k', type=float, default=None,
               help='Restrain receptor to the starting conformation with the provided force constant;'
                ' for minimized energy evaluations. Assumes k is in  kcal/(mol * Angstrom).')

args = p.parse_args()

# Compute Kappa for  implicit solvent Ionic Strength
temperature = 300
solv_dielectric = 78.5
conc = 0.150
kappa = 367.434915*sqrt(conc/(solv_dielectric*temperature))

if not args.out_dir.is_dir():
    args.out_dir.mkdir(parents=True)

ligand = Molecule.from_smiles(args.ligand_smiles)
ligand.generate_conformers(n_conformers=10)
ligand.assign_partial_charges(partial_charge_method='am1bcc', 
                              use_conformers=ligand.conformers)
receptor = Topology.from_pdb(args.receptor_pdb)
receptor_count = len(list(receptor.atoms))
# make topologies, and serialize them
lig_top = ligand.to_topology()
# always do this receptor first!
rl_complex = receptor + lig_top
# rl_complex = lig_top + receptor 
off_serialize(args.out_dir, 'receptor-top', receptor)
off_serialize(args.out_dir, 'ligand-top', lig_top)
off_serialize(args.out_dir, 'complex-top', rl_complex)

ligand_off = off_ForceField(args.ligand_ff, args.ligand_implicit)
# Create the SMIRNOFF template generator with the default installed force field
smirnoff = SMIRNOFFTemplateGenerator(molecules=ligand, forcefield=ligand_off)
# Create an OpenMM ForceField object with AMBER ff14SB and TIP3P with compatible ions
forcefield = ForceField(args.receptor_ff, args.implicit)
# Register the SMIRNOFF template generator
forcefield.registerTemplateGenerator(smirnoff.generator)

# make systems from each of the topologies above
receptor_sys = forcefield.createSystem(receptor.to_openmm(), implicitSolventKappa=kappa)
ligand_sys = forcefield.createSystem(lig_top.to_openmm(), implicitSolventKappa=kappa)
rl_complex_ommt = rl_complex.to_openmm()
complex_sys = forcefield.createSystem(rl_complex_ommt, implicitSolventKappa=kappa)
# optionally add receptor restraints
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



