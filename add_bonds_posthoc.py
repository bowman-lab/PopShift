from pick_align_frames import floatpair, add_bonds_two_cuts
import argparse
import loos
import multiprocessing
from pathlib import Path
from functools import partial

def read_find_write_fused(heavy_cut: float, hydro_cut: float, structure: Path):
    s = loos.createSystem(str(structure))
    add_bonds_two_cuts(s, heavy_cut, hydro_cut)
    with structure.open('w') as f:
        f.write(str(loos.PDB.fromAtomicGroup(s)))

def read_clear_find_write_fused(heavy_cut: float, hydro_cut: float, structure: Path):
    s = loos.createSystem(str(structure))
    s.clearBonds()
    add_bonds_two_cuts(s, heavy_cut, hydro_cut)
    with structure.open('w') as f:
        f.write(str(loos.PDB.fromAtomicGroup(s)))


p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

p.add_argument('models', nargs='+', type=Path,
               help='Models to add bonds to using cutoffs.')
p.add_argument('--cutoffs', '-c', default=(1.95, 1.35), type=floatpair,
               help='Provide a cutoff for heavy atom bonds, and another for bonds to hydrogen.')
p.add_argument('--clear-bonds', action=argparse.BooleanOptionalAction,
               help='If thrown, remove connectivity information from model. If incorrect bonds present, but bonds'
                    ' are needed, use this in conjunction with "--find-bonds" to first clear old bonds then assign'
                    ' new ones.')
p.add_argument('--nprocs', '-n', type=int, default=1,
               help='Use this to control the number of multiprocessing threads deployed. Parallel across structures.')

if __name__ == '__main__':
    args = p.parse_args()
    print(args)
    if args.clear_bonds:
        structure_operator = partial(read_clear_find_write_fused, *args.cutoffs)
    else:
        structure_operator = partial(read_find_write_fused, *args.cutoffs)

    pool = multiprocessing.Pool(args.nprocs)
    pool.map(structure_operator, args.models)
    pool.close()
