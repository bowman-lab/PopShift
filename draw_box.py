import loos
import argparse
from sys import stdout

p = argparse.ArgumentParser(argument_default=True)
p.add_argument('center_x', type=float, help='Center x coordinate')
p.add_argument('center_y', type=float, help='Center y coordinate')
p.add_argument('center_z', type=float, help='Center z coordinate')
p.add_argument('edge_length_x', type=int, help='X edge length')
p.add_argument('edge_length_y', type=int, help='Y edge length')
p.add_argument('edge_length_z', type=int, help='Z edge length')
p.add_argument('-g', '--grid-spacing', type=float, default=1.0,
               help='Grid spacing for autodock-style tools; will be multiplied '
                    'through each edge length.')
args = p.parse_args('0 0 0 8 8 8'.split())


box = loos.AtomicGroup()
center_crd = loos.GCoord(args.center_x, args.center_y, args.center_z)
center = loos.Atom(1, 'center', center_crd)
gs = args.grid_spacing
half_x = (args.edge_length_x / 2) * gs
half_y = (args.edge_length_y / 2) * gs
half_z = (args.edge_length_z / 2) * gs
longest = max(half_x, half_y, half_z) * 2
updown = dict(u=1, d=-1)
count = 1
for updown_x in updown:
    for updown_y in updown:
        for updown_z in updown:
            count += 1
            atom_name_fstring = f'{updown_x}{updown_y}{updown_z}'
            xf, yf, zf = updown[updown_x], updown[updown_y], updown[updown_z]
            corner_crd = loos.GCoord(xf * half_x, yf * half_y, zf * half_z) \
                         + center_crd
            corner = loos.Atom(count, atom_name_fstring, corner_crd)

            box.append(corner)
box.findBonds(longest + gs/2)
box.append(center)
pdb = loos.PDB.fromAtomicGroup(box)
stdout.write(str(pdb))
