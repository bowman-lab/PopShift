import json
import multiprocessing as mp
import argparse
import re
from pathlib import Path

def extract_score_from_vina_pdbqt(pdbqt: Path, search_re=re.compile(r'^REMARK VINA RESULT')):
    


parser = argparse.ArgumentParser()
parser.add_argument('ligand_paths', nargs='+',
                    help='Path(s) to docked ligand directories, containing PDBQTs.')

args = parser.parse_args()
for ligand_path in args.ligand_paths:

