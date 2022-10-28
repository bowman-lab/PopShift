import loos
import argparse
from loos import pyloos
from pathlib import Path
from loos.pyloos import options

p = options.LoosOptions("Read a ligand and a receptor pdb (as from docking); "
                        "Align and compute RMSD to reference.")
p.modelSelectionOptions()
p.add_argument('-H', action=argparse.BooleanOptionalAction,
               help='If thrown, use hydrogens; otherwise add " && ! hydrogen" '
                    'to selections.')
p.add_argument('--translate-selection', type=str, default=None,
               help='If provided, use this selection to subset which atoms are '
                    'translated by alignment. Default uses primary selection.,')
p.add_argument('--rmsd-selection', type=str, default=None,
               help='If provided, use this selection to subset which atoms are '
                    'used for the RMSD calculation. Must be a viable selection '
                    'on ref and on sample pdbs.')
p.add_argument('--receptors', '-R', type=Path, nargs='+', required=True,
               help='Receptors in order they should be matched to ligands.')
p.add_argument('--ligands', '-L', type=Path, nargs='+', required=True,
                help='Ligands in order they should be matched to receptors.')
p.add_argument('--matched-names', '-m', action=argparse.BooleanOptionalAction,
               help='If thrown, assume arguments to --ligands and --receptors '
                    'are directories containing PDBs with matching base-names.')

args = p.parse_args()
