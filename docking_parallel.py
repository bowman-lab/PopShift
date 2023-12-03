"""
Using Jug, run docking calculations across a collection of picked and prepped frames.
    Copyright (C) 2023 Louis G. Smith and Borna Novak

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""
import jug
from vina import Vina
import argparse
import subprocess as sp
from pathlib import Path
from functools import partial


def coordreader(s, delim=','):
    s = s.replace('\\', '')
    try:
        x, y, z = map(float, s.split(delim))
        return x, y, z
    except:
        raise argparse.ArgumentTypeError(
            "Center-coordinates and dimensions must be specified as 'x,y,z'")


def intrange(s, delim=','):
    try:
        start_ind, end_ind = map(int, s.split(delim))
        return start_ind, end_ind
    except:
        raise argparse.ArgumentTypeError(
            'Index ranges must be specified as "start_ind,end_ind"')


plantsconfig_template = """
# scoring function and search settings
scoring_function chemplp
search_speed speed1


# input
protein_file {receptor}
ligand_file {ligand}
aco_ants 40

# output
output_dir {output_dir}

# write single mol2 files (e.g. for RMSD calculation)
write_multi_mol2 0

# binding site definition
bindingsite_center {binding_center}
bindingsite_radius {binding_radius}


# cluster algorithm
cluster_structures 10
cluster_rmsd 2.0
"""


@jug.TaskGenerator
def dock_vina(box_center, box_size, exhaustiveness, receptor_path, ligand_path, output_path):
    v = Vina(sf_name='vina', cpu=exhaustiveness)
    v.set_receptor(str(receptor_path))
    v.set_ligand_from_file(str(ligand_path))
    v.compute_vina_maps(center=box_center, box_size=box_size)
    v.dock(exhaustiveness=exhaustiveness)
    v.write_poses(str(output_path), n_poses=1, overwrite=True)
    return True


@jug.TaskGenerator
def dock_smina(box_center, box_size, exhaustiveness, receptor_path, ligand_path, output_path):
    return sp.run(['smina', '--receptor', str(receptor_path), '--ligand', str(ligand_path),
                   '--center_x', f'{box_center[0]}',
                   '--center_y', f'{box_center[1]}',
                   '--center_z', f'{box_center[2]}',
                   '--size_x', f'{box_size[0]}',
                   '--size_y', f'{box_size[1]}',
                   '--size_z', f'{box_size[2]}',
                   '--exhaustiveness', f'{exhaustiveness}',
                   '--cpu', '1',
                   '--num_modes', '1',
                   '--out', str(output_path)])


# make a functor to hold plants exe and plants template.
class plants_docker:
    def __init__(self, plants_exe_path, plants_template, smina=False):
        self.plants_exe_path = plants_exe_path
        self.plants_template = plants_template
        self.smina = smina

    def __call__(self, binding_center, box_size, receptor_path, ligand_path, output_path, **kwargs):
        plants_dir = output_path.parent/'plants'
        plants_dir.mkdir(parents=True, exist_ok=True)
        plants_conf_p = plants_dir/'plantsconfig'
        plants_out = plants_dir/'results'
        binding_center_str = ' '.join(map(str, binding_center))
        plants_conf = self.plants_template.format(binding_center=binding_center_str, binding_radius=box_size[0]/2,
                            receptor=receptor_path, ligand=ligand_path, output_dir=plants_out, **kwargs)
        plants_conf_p.write_text(plants_conf)
        plants_exit = sp.run(f"{self.plants_exe_path} --mode screen {plants_conf_p}".split())
        # This gets the highest ranking pose, of 10
        plants_pose = next(plants_out.glob('*entry_*_conf_01.mol2'))
        if self.smina:
            if plants_exit.returncode != 0:
                print('Plants failed:', plants_dir)
                return plants_exit
            smina_out = output_path.parent/'smina-min.out'
            smina_exit = sp.run(f'smina -r {receptor_path} -l {plants_pose} -o {output_path} --minimize --accurate_line | tee {smina_out}'.split())
            return smina_exit
        return plants_exit



docking_methods = {
    'vina': dock_vina,
    'smina': dock_smina,
    'plants': None
}

if __name__ == '__main__' or jug.is_jug_running():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

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
    parser.add_argument('-d', '--docking-algorithm', default='vina',
                        choices=docking_methods.keys(),
                        help='Pick which docking algorithm to use.')
    parser.add_argument('--plants-path', default=None, type=Path,
                        help='If docking with plants, provide plants path.')
    parser.add_argument('--rescore-smina', default=True, action=argparse.BooleanOptionalAction,
                        help='If docking with PLANTS, rescore the best pose with SMINA.')
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
        output_paths = [Path(args.out_dir + '-{}'.format(r))
                        for r in range(start_ind, end_ind)]
    else:
        output_paths = [Path(args.out_dir)]

    # Make output dirs.
    for p in output_paths:
        p.mkdir(exist_ok=True, parents=True)

    # if docking using SMINA, get pdbqts; otherwise, use mol2s, and prep docking function
    if args.docking_algorithm == 'plants':
        if args.plants_path.is_file():
            plants_functor = plants_docker(plants_exe_path=args.plants_path, plants_template=plantsconfig_template, smina=args.rescore_smina)
            docking_methods['plants'] = plants_functor
        else:
            print(f'ERROR: Plants path is {args.plants_path}, which is not a file.')
            exit(1)
    # uses recursive glob. Must be sorted to get same order across runs.
    frame_paths = sorted(path_receptor.rglob('*.pdbqt'))

    for run_path in output_paths:
        for ligand in ligand_paths:
            lig_path = Path(ligand)
            ligand_name = lig_path.stem
            lig_output_path = run_path / ligand_name
            for frame_path in frame_paths:
                docked_lig_path = lig_output_path.joinpath(
                    *frame_path.parts[-2:])
                docked_dir_path = docked_lig_path.parent
                if not docked_dir_path.is_dir():
                    docked_dir_path.mkdir(exist_ok=True, parents=True)
                if not args.dry_run:
                    jug.Task(docking_methods[args.docking_algorithm](
                        args.box_center,
                        args.box_size,
                        args.exhaustiveness,
                        frame_path,
                        lig_path,
                        docked_lig_path
                    ))
