from jug import TaskGenerator, Task
from vina import Vina
import argparse
import subprocess as sp
from pathlib import Path


def coordreader(s, delim=','):
    s = s.replace('\\', '')
    try:
        x, y, z = map(float, s.split(delim))
        return x, y, z
    except:
       raise argparse.ArgumentTypeError("Center-coordinates and dimensions must be specified as 'x,y,z'")


def intrange(s, delim=','):
    try:
        start_ind, end_ind = map(int, s.split(delim))
        return start_ind, end_ind
    except:
        raise argparse.ArgumentTypeError('Index ranges must be specified as "start_ind,end_ind"')


@TaskGenerator
def dock_vina(box_center, box_size, exhaustiveness, receptor_path, ligand_path, output_path):
    v = Vina(sf_name='vina', cpu=exhaustiveness)
    v.set_receptor(str(receptor_path))
    v.set_ligand_from_file(str(ligand_path))
    v.compute_vina_maps(center=box_center, box_size=box_size)
    v.dock(exhaustiveness=exhaustiveness)
    v.write_poses(str(output_path), n_poses=1, overwrite=True)
    return True


@TaskGenerator
def dock_smina(box_center, box_size, exhaustiveness, receptor_path, ligand_path, output_path):
    return sp.run(['smina', '--receptor', str(receptor_path), '--ligand', str(ligand_path), \
        '--center_x', f'{box_center[0]}', \
        '--center_y', f'{box_center[1]}', \
        '--center_z', f'{box_center[2]}', \
        '--size_x', f'{box_size[0]}', \
        '--size_y', f'{box_size[1]}', \
        '--size_z', f'{box_size[2]}', \
        '--exhaustiveness', f'{exhaustiveness}', \
        '--cpu', '1', \
        '--num_modes','1', \
        '--out', str(output_path)])


docking_methods = {
    'vina': dock_vina,
    'smina': dock_smina,
}


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('receptor_dir',
                    help='Path to protein directory')
parser.add_argument('out_dir',
                    help='Path to the output. By convention, the name of the docking run including info like box size '
                         'if multiple are being tested.')
parser.add_argument('box_center', type=coordreader,
                    help='Comma delimited string listing x,y,z of box center. If x is negative start with backslash ("\-5,5,3")')
parser.add_argument('box_size', type=coordreader,
                    help='Comma delimited string listing lx,ly,lz as the lengths of the x, y and z box-sides.')
parser.add_argument('ligand_list', nargs="+",
                    help='Path(s) to ligand pdbqts. Alternatively a txt file with a ligand path on each line.')
parser.add_argument('-r', '--replicas', type=intrange, default=None,
                    help='Number of replica docking runs to perform.')
parser.add_argument('-e', '--exhaustiveness', type=int, default=32,
                    help='AutoDock-Vina exhaustiveness parameter. Threads used proportional to this value.')
parser.add_argument('--protein-prefix', type=str, default='frame00',
                    help='String to prefix output pdbqts with.')
parser.add_argument('-d', '--docking_algorithm', default='vina',
                    choices=docking_methods.keys(),
                    help='Pick which docking algorithm to use.')
parser.add_argument('-s', '--symlink-receptors', action=argparse.BooleanOptionalAction,
                    help='Create relative symlinks for receptor PDBs into each docked ligand dir.')
parser.add_argument('-t', '--top-dir', type=Path, default=Path.cwd(),
                    help='Set a top directory for relative symlinks.')
parser.add_argument('--dry-run', action=argparse.BooleanOptionalAction,
                    help="If thrown, don't actually run docking; just create directories and (optionally) symlinks.")


args = parser.parse_args()
if len(args.ligand_list) == 1:
    ligand_list_path = Path(args.ligand_list[0])
    if ligand_list_path.suffix != '.pdbqt':
        ligand_paths = ligand_list_path.read_text().split()
    else:
        ligand_paths = args.ligand_list
else:
    ligand_paths = args.ligand_list
path_receptor = Path(args.receptor_dir)

# Figure out whether we need to index them by replica count.
if args.replicas:
    start_ind, end_ind = args.replicas
    output_paths = [Path(args.out_dir + '-{}'.format(r)) for r in range(start_ind, end_ind)]
else:
    output_paths = [Path(args.out_dir)]

# Make output dirs.
for p in output_paths:
    p.mkdir(exist_ok=True, parents=True)

# uses recursive glob. Must be sorted to get same order across runs.
frame_paths = sorted(path_receptor.rglob('*.pdbqt'))


for run_path in output_paths:
    for ligand in ligand_paths:
        lig_path = Path(ligand)
        ligand_name = lig_path.stem
        lig_output_path = run_path / ligand_name
        for frame_path in frame_paths:
            docked_lig_path = lig_output_path.joinpath(*frame_path.parts[-2:])
            docked_dir_path = docked_lig_path.parent
            if not docked_dir_path.is_dir():
                docked_dir_path.mkdir(exist_ok=True, parents=True)
            if not args.dry_run:
                docking_methods[args.docking_algorithm](
                    args.box_center,
                    args.box_size,
                    args.exhaustiveness,
                    frame_path,
                    lig_path,
                    docked_lig_path
                )

