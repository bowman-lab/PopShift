from enspara import ra
import numpy as np
import multiprocessing as mp
import argparse
import re
from pathlib import Path
import pickle

def extract_score_from_vina_pdbqt(pdbqts: list, search_re=re.compile(r'^REMARK VINA RESULT')):
    outlist = []
    for pdbqt in pdbqts:
        with pdbqt.open() as f:
            for line in f:
                if search_re.match(line):
                    outlist.append(float(line.split(maxsplit=4)[3]))
                    break
    return np.array(outlist, dtype='f4')

def extract_score_from_smina_pdbqt(pdbqts: list, search_re=re.compile(r'^REMARK minimizedAffinity')):
    outlist = []
    for pdbqt in pdbqts:
        with pdbqt.open() as f:
            for line in f:
                if search_re.match(line):
                    outlist.append(float(line.split(maxsplit=4)[2]))
                    break
    return np.array(outlist, dtype='f4')



number_pattern = re.compile('\d+')


def rip_last_number(p: Path, numpattern=number_pattern):
    return int(numpattern.findall(p.stem)[-1])


extract_methods = {
    'vina': extract_score_from_vina_pdbqt,
    'smina': extract_score_from_smina_pdbqt
}

parser = argparse.ArgumentParser()
parser.add_argument('docking_runs', nargs='+',
                    help='Path(s) to docked ligand directories, containing PDBQTs.')
parser.add_argument('--centers', action=argparse.BooleanOptionalAction, default=False,
                    help='If thrown, expect that all pdbqts for a given ligand will be '
                         'together in one large directory, with state numbers in file '
                         'names (4x mydockrun/myligand/state000123.pdbqt).')
parser.add_argument('--replicas', action=argparse.BooleanOptionalAction, default=False,
                    help='If thrown, expect that there will be an additional level of depth '
                         'of the directory structure named something like "replica0".')
parser.add_argument('--nprocs', '-n', type=int, default=1,
                    help='number of processors to use for multiprocessing module.')
parser.add_argument('--result-type', '-t', choices=extract_methods.keys(), default='vina',
                    help='Which type of output file you are trying to extract a score from.')
parser.add_argument('--extract-to', type=Path, default=None,
                    help='If provided, write extracted scores to this directory instead of customary one.')

args = parser.parse_args()
pool = mp.Pool(args.nprocs)
em = extract_methods[args.result_type]
for dock_run in args.docking_runs:
    dock_path = Path(dock_run)
    if args.extract_to:
        out_path = args.extract_to
    else:
        out_path = dock_path.parent / 'extracted_scores' / dock_path.name
    out_path.mkdir(exist_ok=True, parents=True)
    for ligand_path in dock_path.iterdir():
        if args.replicas:
            rep_paths = sorted((rep_path for rep_path in ligand_path.iterdir()),
                               key=rip_last_number)
            for rep_ix, rep_path in enumerate(rep_paths):
                if args.centers:
                    result_paths = [sorted((sample_path for sample_path in rep_path.rglob('*.pdbqt')),
                                           key=rip_last_number) for rep_path in rep_paths]

                else:
                    state_paths = sorted((state_path for state_path in rep_path.iterdir()),
                                         key=rip_last_number)
                    result_paths = [sorted((result_path for result_path in state_path.iterdir()),
                                           key=lambda x: tuple(map(int, number_pattern.findall(str(x))[-2:])))
                                    for state_path in state_paths]
                extracted_results = pool.map(em, result_paths)
                for i, result in enumerate(extracted_results):
                    if len(result) != 1:
                        print(result_paths[i], result)
                outpre = out_path / (ligand_path.stem + f'-{rep_ix}')
                r = ra.RaggedArray(extracted_results)
                ra.save(str(outpre.with_suffix('.h5')), r)
                with outpre.with_suffix('.pickle').open('wb') as f:
                    pickle.dump(result_paths, f)
        else:
            if args.centers:
                result_paths = sorted((state_path for state_path in ligand_path.rglob('*.pdbqt')),
                                      key=rip_last_number)
            else:
                state_paths = sorted((state_path for state_path in ligand_path.iterdir()),
                                     key=lambda x: int(x.stem))
                result_paths = [sorted((sample_path for sample_path in state_path.rglob('*.pdbqt')),
                                       key=lambda x: tuple(map(int, number_pattern.findall(str(x))[-2:])))
                                for state_path in state_path]
            extracted_results = pool.map(em, result_paths)
            for i, result in enumerate(extracted_results):
                if len(result) != 1:
                    print(result_paths[i], result)
            outpre = out_path/ligand_path.stem
            r = ra.RaggedArray(extracted_results)
            ra.save(str(outpre.with_suffix('.h5')), r)
            with outpre.with_suffix('.pickle').open('wb') as f:
                pickle.dump(result_paths, f)

pool.close()
