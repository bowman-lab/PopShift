from jug import TaskGenerator
from vina import Vina
from functools import partial
import os
import argparse
import glob
import subprocess as sp
import multiprocessing as mp
from pathlib import Path

@TaskGenerator
def coordreader(s, delim=','):
    try:
        x, y, z = map(float, s.split(delim))
        return x, y, z
    except:
       raise argparse.ArgumentTypeError("Center-coordinates and dimensions must be specified as 'x,y,z'")


@TaskGenerator
def dock_vina(box_center, box_size, exhaustiveness, receptor, ligand, ligand_name, output):
    receptor_name = Path(receptor).stem
    v = Vina(sf_name='vina',cpu=exhaustiveness)
    v.set_receptor(receptor)
    v.set_ligand_from_file(ligand)
    v.compute_vina_maps(center=box_center, box_size=box_size)  # Don't forget to set this
    v.dock(exhaustiveness=exhaustiveness)
    fname = '{}/{}-{}-{:.2f}_{:.2f}_{:.2f}_{}_{}_{}_vina.pdbqt'.format(
        output,
        receptor_name,
        ligand_name,
        box_center[0],
        box_center[1],
        box_center[2],
        box_size[0],
        box_size[1],
        box_size[2]
    )
    v.write_poses(fname, n_poses=1,overwrite=True)

@TaskGenerator
def dock_smina(box_center, box_size, exhaustiveness, receptor, ligand, ligand_name, output):
    receptor_name = Path(receptor).stem
    fname = '{}/{}-{}-{:.2f}_{:.2f}_{:.2f}_{}_{}_{}_smina.pdbqt'.format(
        output,
        receptor_name,
        ligand_name,
        box_center[0],
        box_center[1],
        box_center[2],
        box_size[0],
        box_size[1],
        box_size[2]
    )
    sp.run(['smina', '--receptor', receptor, '--ligand', ligand, \
        '--center_x', '%s' % box_center[0], \
        '--center_y', '%s' % box_center[1], \
        '--center_z', '%s' % box_center[2], \
        '--size_x', '%s' % box_size[0], \
        '--size_y', '%s' % box_size[1], \
        '--size_z', '%s' % box_size[2], \
        '--exhaustiveness', '%s' % exhaustiveness, \
        '--cpu', '%s' % exhaustiveness, \
        '--num_modes','1', \
        '--out', fname])

docking_methods = {
    'vina': dock_vina,
    'smina': dock_smina,
}

parser = argparse.ArgumentParser()

parser.add_argument('protein_dir',
                    help='Path to protein directory')
parser.add_argument('out_dir',
                    help='Path to the output. By convention, the name of the docking run including info like box size '
                         'if multiple are being tested.')
parser.add_argument('box_center', type=coordreader,
                    help='Comma delimited string listing x,y,z of box center')
parser.add_argument('box_size', type=coordreader,
                    help='Comma delimited string listing lx,ly,lz as the lengths of the x, y and z box-sides.')
parser.add_argument('ligand_list', nargs="+",
                    help='Path(s) to ligand pdbqts. Alternatively a txt file with a ligand path on each line.')
parser.add_argument('-r', '--replicas', type=int, default=1,
                    help='Number of replica docking runs to perform. Default: 1')
parser.add_argument('-e', '--exhaustiveness', type=int, default=32,
                    help='AutoDock-Vina exhaustiveness parameter. Threads used proportional to this value. Default: 32')
parser.add_argument('--protein-prefix', type=str, default='frame00',
                    help='String to prefix output pdbqts with. Default: frame00')
parser.add_argument('-d','--docking_algorithm', default='vina',
                    choices=docking_methods.keys(),
                    help='Pick which docking algorithm to use. Default: vina')


args = parser.parse_args()
if len(args.ligand_list) == 1:
    ligand_list_path = Path(args.ligand_list[0])
    if ligand_list_path.suffix != '.pdbqt':
        ligand_paths = ligand_list_path.read_txt().split()
    else:
        ligand_paths = args.ligand_list
else:
    ligand_paths = args.ligand_list
path_prot = args.protein_dir
path_output = args.output
n_procs = mp.cpu_count() #TODO Need to check whether this takes cpu count of the node, or the number of cpus requested
pool = mp.Pool(processes=int(n_procs/args.exhaustiveness)) #if the above function gets the # of CPUs requested, this should work

try:
    os.makedirs('%s' % path_output)
except FileExistsError:
    pass

for ligand in ligand_paths:
    ligand_name = Path(ligand).stem
    for replica in range(args.replicas):
        try:
            os.makedirs('%s/%s/replica%s' % (path_output,ligand_name,replica))
        except FileExistsError:
            pass
        frames = sorted(glob.glob('%s/*pdbqt' % path_prot))
        ligand_files = [i for i in [ligand] for l in range(len(frames))]
        ligand_names = [i for i in [ligand_name] for l in range(len(frames))]
        output_path = '%s/%s/replica%s' % (path_output,ligand_name,replica)
        output_paths = [i for i in [output_path] for l in range(len(frames))]
        loaded_dock = partial(docking_methods[args.docking_algorithm], args.box_center, args.box_size, args.exhaustiveness)
        arguments = zip(frames, ligand_files, ligand_names, output_paths)
        list(pool.starmap(loaded_dock, arguments))

pool.close()

