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
import argparse as ap
import subprocess as sp
from pathlib import Path
from functools import partial


def coordreader(s, delim=','):
    s = s.replace('\\', '')
    try:
        x, y, z = map(float, s.split(delim))
        return x, y, z
    except:
        raise ap.ArgumentTypeError(
            "Center-coordinates and dimensions must be specified as 'x,y,z'")


def intrange(s, delim=','):
    try:
        start_ind, end_ind = map(int, s.split(delim))
        return start_ind, end_ind
    except:
        raise ap.ArgumentTypeError(
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
def dock_vina(box_center, box_size, receptor_path, ligand_path, output_path, exhaustiveness=32):
    from vina import Vina
    v = Vina(sf_name='vina', cpu=exhaustiveness)
    v.set_receptor(str(receptor_path))
    v.set_ligand_from_file(str(ligand_path))
    v.compute_vina_maps(center=box_center, box_size=box_size)
    v.dock(exhaustiveness=exhaustiveness)
    v.write_poses(str(output_path), n_poses=1, overwrite=True)
    return True


@jug.TaskGenerator
def dock_smina(box_center, box_size,  receptor_path, ligand_path, output_path, exhaustiveness=32,cpu=1):
    return sp.run(['smina', '--receptor', str(receptor_path), '--ligand', str(ligand_path),
                   '--center_x', f'{box_center[0]}',
                   '--center_y', f'{box_center[1]}',
                   '--center_z', f'{box_center[2]}',
                   '--size_x', f'{box_size[0]}',
                   '--size_y', f'{box_size[1]}',
                   '--size_z', f'{box_size[2]}',
                   '--exhaustiveness', f'{exhaustiveness}',
                   '--cpu', f'{cpu}',
                   '--num_modes', '1',
                   '--out', str(output_path)])

@jug.TaskGenerator
def dock_gnina(box_center, box_size, receptor_path, ligand_path, output_path, 
               num_modes=1, cnn_scoring='rescore', addH=0, exhaustiveness=32,
               cnn_freeze_receptor='--cnn_freeze_receptor', 
               pose_sort_order='CNNaffinity', cnn=None, cpu=1):
    if cnn:
        return sp.run(['gnina', '--receptor', str(receptor_path), '--ligand', str(ligand_path),
                        '--center_x', f'{box_center[0]}',
                        '--center_y', f'{box_center[1]}',
                        '--center_z', f'{box_center[2]}',
                        '--size_x', f'{box_size[0]}',
                        '--size_y', f'{box_size[1]}',
                        '--size_z', f'{box_size[2]}',
                        '--cnn_scoring', f'{cnn_scoring}',
                        '--exhaustiveness', f'{exhaustiveness}',
                        '--num_modes', f'{num_modes}',
                        '--pose_sort_order', pose_sort_order,
                        '--cnn', cnn,
                        '--addH', f'{addH}',
                        '--cpu', f'{cpu}',
                        f'{cnn_freeze_receptor}',
                        '--out', str(output_path)])
    else:
        return sp.run(['gnina', '--receptor', str(receptor_path), '--ligand', str(ligand_path),
                   '--center_x', f'{box_center[0]}',
                   '--center_y', f'{box_center[1]}',
                   '--center_z', f'{box_center[2]}',
                   '--size_x', f'{box_size[0]}',
                   '--size_y', f'{box_size[1]}',
                   '--size_z', f'{box_size[2]}',
                   '--cnn_scoring', f'{cnn_scoring}',
                   '--exhaustiveness', f'{exhaustiveness}',
                   '--num_modes', f'{num_modes}',
                   '--pose_sort_order', pose_sort_order,
                   '--addH', f'{addH}',
                   '--cpu', f'{cpu}',
                   f'{cnn_freeze_receptor}',
                   '--out', str(output_path)])


# make a functor to hold plants exe and plants template.
class plants_docker:
    def __init__(self, plants_exe_path, plants_template, smina=False, overwrite=True):
        self.plants_exe_path = plants_exe_path
        self.plants_template = plants_template
        self.smina = smina
        self.overwrite = overwrite

    def __call__(self, binding_center, box_size, receptor_path, ligand_path, output_path, **kwargs):
        plants_dir = output_path.parent/'plants'
        plants_dir.mkdir(parents=True, exist_ok=True)
        plants_conf_p = plants_dir/'plantsconfig'
        plants_out = plants_dir/'results'
        if plants_out.is_dir():
            if self.overwrite:
                rm_exit = sp.run(f'rm -r {plants_out}'.split())
                print(rm_exit, flush=True)
        binding_center_str = ' '.join(map(str, binding_center))
        # write conf such that plants can be run from plants_out as cwd.
        plants_conf = self.plants_template.format(binding_center=binding_center_str, 
                                                  binding_radius=box_size[0]/2, 
                                                  receptor=receptor_path.resolve(), 
                                                  ligand=ligand_path.resolve(), 
                                                  output_dir=plants_out.name, **kwargs)
        print(plants_conf, flush=True)
        plants_conf_p.write_text(plants_conf)
        plants_exit = sp.run(f"{self.plants_exe_path} --mode screen {plants_conf_p.name}".split(), cwd=plants_conf_p.parent)
        # This gets the highest ranking pose, of 10
        plants_pose = next(plants_out.glob('*entry_*_conf_01.mol2'))
        if self.smina:
            print('re-minimizing with smina:', plants_pose)
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
    'gnina': dock_gnina,
    'plants': None
}

if __name__ == '__main__' or jug.is_jug_running():
    parser = ap.ArgumentParser(
        formatter_class=ap.ArgumentDefaultsHelpFormatter)

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
                        help='Path(s) to ligand pdbqts, pdbs, or sdfs. Alternatively a txt file with a ligand path on each line.')
    parser.add_argument('--ligand-name-from-dir', action=ap.BooleanOptionalAction, default=False,
                        help='If thrown, assume that the ligand name for outfiles should come from the directory' 
                        ' above the ligand specifying file--for GB rescoring.')
    parser.add_argument('-r', '--replicas', type=intrange, default=None,
                        help='Number of replica docking runs to perform.')
    parser.add_argument('-c', '--centers', action=ap.BooleanOptionalAction, default=False,
                        help='If thrown, assume receptor structures are all in one large directory, rather '
                        'than subdirs corresponding to state identity. Write docked poses in the same fashion.')
    parser.add_argument('-e', '--exhaustiveness', type=int, default=32,
                        help='AutoDock-Vina exhaustiveness parameter. Threads used proportional to this value.')
    parser.add_argument('--protein-prefix', type=str, default='frame00',
                        help='String to prefix output pdbqts with.')
    parser.add_argument('-d', '--docking-algorithm', default='vina',
                        choices=docking_methods.keys(),
                        help='Pick which docking algorithm to use.')
    parser.add_argument('--cpu', type=int, default=1,
                        help='Number of CPUs _per docking job_. Leave at 1 unless using a method bottlenecked '
                        'by the availability of an accelerator like a GPU.')
    parser.add_argument('--num_modes', type=int, default=1,
                        help='Number of binding modes to ask for.')
    parser.add_argument('--cnn-scoring', type=str, default='rescore',
                        help='If using GNINA for docking, pass this argument through to call to gnina.')
    parser.add_argument('--cnn-freeze-receptor', default=True, action=ap.BooleanOptionalAction,
                        help='If using GNINA, pass this onto the call.')
    parser.add_argument('--plants-path', default=None, type=Path,
                        help='If docking with plants, provide plants path.')
    parser.add_argument('--rescore-smina', default=True, action=ap.BooleanOptionalAction,
                        help='If docking with PLANTS, rescore the best pose with SMINA.')
    parser.add_argument('-s', '--symlink-receptors', action=ap.BooleanOptionalAction,
                        help='Create relative symlinks for receptor PDBs into each docked ligand dir.')
    parser.add_argument('-t', '--top-dir', type=Path, default=Path.cwd(),
                        help='Set a top directory for relative symlinks.')
    parser.add_argument('--dry-run', action=ap.BooleanOptionalAction,
                        help="If thrown, don't actually run docking; just create directories and (optionally) symlinks.")

    args = parser.parse_args()
    if len(args.ligand_list) == 1:
        ligand_list_path = Path(args.ligand_list[0])
        if ligand_list_path.suffix == '.pdbqt' or ligand_list_path.suffix == '.sdf' or ligand_list_path.suffix == '.pdb':
            ligand_paths = args.ligand_list
        else:
            ligand_paths = ligand_list_path.read_text().split()
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
    # get dock task generator set up.
    dock_algo_name = args.docking_algorithm
    dock_algo = docking_methods[args.docking_algorithm]

    # if docking using SMINA, get pdbqts; otherwise, use mol2s, and prep docking function
    if dock_algo_name == 'plants':
        if args.plants_path.is_file():
            plants_functor = plants_docker(plants_exe_path=args.plants_path, plants_template=plantsconfig_template, smina=args.rescore_smina)
            docking_methods['plants'] = plants_functor
        else:
            print(f'ERROR: Plants path is {args.plants_path}, which is not a file.')
            exit(1)


    # uses recursive glob. Must be sorted to get same order across runs.
    if dock_algo_name == 'gnina':
        frame_paths = sorted(map(lambda x: x.with_suffix('.pdb'), path_receptor.rglob('*.pdb')))
        if args.cnn_freeze_receptor:
            dock_algo = partial(dock_algo, num_modes=args.num_modes, cnn_scoring=args.cnn_scoring, 
                                cnn_freeze_receptor='--cnn_freeze_receptor', cpu=args.cpu, exhaustiveness=args.exhaustiveness)
        else:
            dock_algo = partial(dock_algo, num_modes=args.num_modes, cnn_scoring=args.cnn_scoring, 
                                cpu=args.cpu, exhaustiveness=args.exhaustiveness)
    else: 
        frame_paths = sorted(path_receptor.rglob('*.pdbqt'))
    # if docking to centers, or files without directory structure, 
    # only use file name for pose writing
    if args.centers:  
        child_include_lvl = -1 
    else: # include the directory above the frame file name, as this will record state id.
        child_include_lvl = -2
    for run_path in output_paths:
        for ligand in ligand_paths:
            lig_path = Path(ligand)
            if args.ligand_name_from_dir:
                ligand_name = lig_path.parent.stem
            else:
                ligand_name = lig_path.stem
            lig_output_path = run_path / ligand_name
            for frame_path in frame_paths:
                if dock_algo_name == 'gnina': 
                    docked_lig_path = lig_output_path.joinpath(*frame_path.parts[child_include_lvl:]).with_suffix('.sdf')
                    docked_dir_path = docked_lig_path.parent
                elif dock_algo_name == 'smina' or dock_algo_name == 'vina': 
                    docked_lig_path = lig_output_path.joinpath(*frame_path.parts[child_include_lvl:]).with_suffix('.pdbqt')
                    docked_dir_path = docked_lig_path.parent
                if not docked_dir_path.is_dir():
                    docked_dir_path.mkdir(exist_ok=True, parents=True)
                if not args.dry_run:
                    dock_algo(
                        args.box_center,
                        args.box_size,
                        frame_path,
                        lig_path,
                        docked_lig_path
                    )
