#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import numpy as np
import os
import sys
import argparse
import argcomplete
from tqdm import tqdm

import MDAnalysis as mda
import myMDAnalysis.ContactMapAnalysis.ContactMapAnalysisAPI as cma
import scatter.scatterMDA as scatter

description="""Component Intermediate Incoherent structure factors.
Ouputs HDF5 files:
- sf.h5  Incoherent structure factor (ISF)
- sfYY.h5 ISF for bound solvent to solute.
- sfNN.h5 ISF for unbound of bulk solvent.
- sfNY.h5 ISF for binding processes.
- sfYN.h5 ISF for unbinnding processes.
"""
#
#  Required arguments
parser = argparse.ArgumentParser(description=description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('topology_file', type=str,
                    help='topology file. PDB format is OK')
parser.add_argument('trajectory_file', type=str,
                    help='Trajectory file, wrapped coordinates')
#
#  Optional arguments
parser.add_argument('--solvent_sel', default='resname XXXX and name H*',
                    help='selection for solvent atoms')
#  Default solute is polymer plus bound Na atoms
parser.add_argument('--solute_sel',
                    default='(not resname XXXX) or (resname XXXX and name Na*)',
                    help='selection for solute atoms')
help = """maximum distance from solvent atom to solute to consider solute
atom as bound"""
parser.add_argument('--cut_off', default=2.5, help=help)
help = """File to save (if files does not exist) or load (if file exists)
the trajectory flagging for each solvent atom whether it's bound to any
solute atoms"""
parser.add_argument('--save_flags_file', default='bound_solute_flags.npy',
                    help=help)
help = """Triad defining the momentum transfer values
min_val, max_val, step. For instance, 0.1, 2.1, 0.2 will yield Q values
from 0.1 up to and included 1.9. Last value 2.1 is excluded."""
parser.add_argument('--qvalues', type=str, default='"0.1, 2.1, 0.2"',
                    help=help)
help="""Elemental elapsed time, t of I(Q,t), in units of number of
trajectory frames"""
parser.add_argument('--dt', type=int, default=10,
                    help='Elemental elapsed time')
parser.add_argument('--nt', type=int, default=1000,
                    help='Number of elapsed times')
parser.add_argument('--prefix', default='',
                    help='prefix each output ISF with this input sting')
argcomplete.autocomplete(parser)
args = parser.parse_args()

a_universe = mda.Universe(args.topology_file, args.trajectory_file)

# Find when are solvent (water Hydrogen) atoms bound to solute (polymer)
print('Finding solvent atoms bound to solute')
if not os.path.isfile(args.save_flags_file):
    flags = cma.solvent_bound_flag(a_universe, args.solvent_sel,
                                   args.solute_sel, args.cut_off)
    flags = np.array(flags).transpose()  # shape = (nframes, natoms)
    np.save(args.save_flags_file, flags)
else:
    flags = np.load(args.save_flags_file)

# Extract coordinates of solvent Hydrogen atoms
print('Extract coordinates of solvent atoms')
solvent_group = a_universe.select_atoms(args.solvent_sel)
positions = list()
for _ in tqdm(a_universe.trajectory, total=len(a_universe.trajectory)):
    positions.append(solvent_group.positions)

# Find scattering
b = np.ones(len(solvent_group))  # scattering lengths set to unity
q_values = np.arange(*[float(q) for q in args.qvalues.strip('"').split(',')])
# frames separated every 0.1ps
# nt = 1000  number of time points
# dt = 10 separation between time consecutive time points (1ps)
print('Calculate intermediate incoherent structure factors')
sqe = scatter.II(positions, b, q_values, dt=args.dt, nt=args.nt, c1=flags)
# sqe = scatter.II(a_universe.trajectory, b, QQ, dt=10, nt=1000, c1=flags)

for scatter_function in 'sf sfNN sfNY sfYN sfYY'.split():
    file_name = '{}{}.h5'.format(args.prefix, scatter_function)
    scatter.saveIISassenaFormat(sqe[scatter_function], QQ, file_name)

