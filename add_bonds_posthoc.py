from pick_align_frames import floatpair, add_bonds_two_cuts
import argparse
import loos
import multiprocessing
from pathlib import Path

p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter())


p.add_argument('model', nargs='+', type=Path,
                help='Model to add bonds to using cutoffs.')
p.add_argument('--cutoffs', '-c', default=(1.95, 1.35), type=floatpair,
                help='Provide a cutoff for heavy atom bonds, and another for bonds to hydrogen.')
pa.add_argument('--clear-bonds', action=argparse.BooleanOptionalAction,
                    help='If thrown, remove connectivity information from model. If incorrect bonds present, but bonds'
                         ' are needed, use this in conjunction with "--find-bonds" to first clear old bonds then assign'
                         ' new ones.')