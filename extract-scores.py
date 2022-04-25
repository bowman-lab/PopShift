from enspara import ra
import multiprocessing as mp
import argparse
import re
from pathlib import Path

def extract_score_from_vina_pdbqt(pdbqt: Path, search_re=re.compile(r'^REMARK VINA RESULT')):
    with pdbqt.open() as f:
        for line in f:
            if search_re.match(line):
               return float(line.split(maxsplit=4)[3])


extract_methods = {
    'vina': extract_score_from_vina_pdbqt,
    'smina': None
}

parser = argparse.ArgumentParser()
parser.add_argument('docking-runs', nargs='+',
                    help='Path(s) to docked ligand directories, containing PDBQTs.')
parser.add_argument('--nprocs','-n', type=int, default=1,
                    help='number of processors to use for multiprocessing module.')
parser.add_argument('--result-type', '-t', choices=extract_methods.keys(), default='vina',
                    help='Which type of output file you are trying to extract a score from.')

args = parser.parse_args()
pool = mp.Pool(args.nprocs)
em = extract_methods[args.result_type]
for dock_run in args.docking_runs:
    dock_path = Path(dock_run)
    for ligand_path in dock_path.iterdir():
        state_paths = sorted(state_path for state_path in ligand_path.iterdir())
        r = ra.RaggedArray(pool.starmap(em, state_paths))
        r.save(

pool.close()
