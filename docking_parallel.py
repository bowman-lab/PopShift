
from vina import Vina
from functools import partial
import os
import argparse
import glob
import subprocess as sp
import multiprocessing as mp


def coordreader(s, delim=','):
    try:
        x, y, z = map(float, s.split(delim))
        return x, y, z
    except:
       raise argparse.ArgumentTypeError("Center-coordinates and dimensions must be specified as 'x,y,z'")

def dock_vina(box_center, box_size, exhaustiveness, receptor, ligand, ligand_name, output):
    receptor_name = receptor.split('/')[-1]
    receptor_name = receptor_name.split('.')[0]
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

def dock_smina(box_center, box_size, exhaustiveness, receptor, ligand, ligand_name, output):
    receptor_name = receptor.split('/')[-1]
    receptor_name = receptor_name.split('.')[0]
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
parser.add_argument('ligand_dir',
                    help='Path to ligand directory')
parser.add_argument('protein_dir',
                    help='Path to protein directory')
parser.add_argument('output',
                    help='Path to the output')
parser.add_argument('box_center', type=coordreader,
                    help='Comma delimited string listing x,y,z of box center')
parser.add_argument('box_size', type=coordreader,
                    help='Comma delimited string listing lx,ly,lz as the lengths of the x, y and z box-sides.')
parser.add_argument('-r', '--replicas', type=int, default=1,
                    help='Number of replica docking runs to perform')
parser.add_argument('-e', '--exhaustiveness', type=int, default=32,
                    help='AutoDock-Vina exhaustiveness parameter. Threads used proportional to this value.')
parser.add_argument('--protein-prefix', type=str, default='frame00',
                    help='String to prefix output pdbqts with.')
parser.add_argument('-d','--docking_algorithm', default='vina',
                    choices=docking_methods.keys(),
                    help='Pick which docking algorithm to use')


args = parser.parse_args()
path_lig = args.ligand_dir
path_prot = args.protein_dir
path_output = args.output
protein_name = path_output.split('/')[-1]
docking_alg = args.docking_algorithm
n_procs = mp.cpu_count()
pool = mp.Pool(processes=4)

try:
    os.makedirs('%s' % path_output)
except FileExistsError:
    pass

for ligand in sorted(glob.glob('%s/*pdbqt' % path_lig)):
    ligand_name = (ligand.split('/')[-1]).split('.')[0]
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
        if docking_alg == 'vina':
            loaded_dock = partial(dock_vina, args.box_center, args.box_size, args.exhaustiveness)
        else:
            loaded_dock = partial(dock_smina, args.box_center, args.box_size, args.exhaustiveness)
        #loaded_dock = partial(args.docking_algorithm, args.box_center, args.box_size, args.exhaustiveness)
        #print(output_paths)
        arguments = zip(frames, ligand_files, ligand_names, output_paths)
        list(pool.starmap(loaded_dock, arguments))

pool.close()

