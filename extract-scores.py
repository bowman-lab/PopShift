from enspara import ra
import numpy as np
import multiprocessing as mp
import argparse
import re
from pathlib import Path


def extract_score_from_vina_pdbqt(pdbqts: list, search_re=re.compile(r'^REMARK VINA RESULT')):
    outlist = []
    for pdbqt in pdbqts:
        with pdbqt.open() as f:
            for line in f:
                if search_re.match(line):
                    outlist.append(float(line.split(maxsplit=4)[3]))
                    break
    return np.array(outlist, dtype='f4')


extract_methods = {
    'vina': extract_score_from_vina_pdbqt,
    'smina': None  # add this later. should be easy.
}

parser = argparse.ArgumentParser()
parser.add_argument('docking_runs', nargs='+',
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
    out_path = dock_path.parent/'extracted_scores'/dock_path.name
    out_path.mkdir(exist_ok=True, parents=True)
    for ligand_path in dock_path.iterdir():
        state_paths = sorted(sorted(sample_path for sample_path in state_path.iterdir())
                             for state_path in ligand_path.iterdir())
        r = ra.RaggedArray(pool.map(em, state_paths))
        ra.save(str(out_path/ligand_path.stem)+'.h5', r)

pool.close()
