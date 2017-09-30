#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import numpy as np
import os
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
fmt = argparse.ArgumentDefaultsHelpFormatter
parser = argparse.ArgumentParser(description=description, formatter_class=fmt)
parser.add_argument('topology_file', type=str,
                    help='topology file. PDB format is OK')
parser.add_argument('trajectory_file', type=str,
                    help='Trajectory file, unwrapped coordinates')
parser.add_argument('save_flags_file', default='bound_solute_flags.npy',
                    help='File containing the flags trajectory for each'\
                         'solvent atom determining whether it is bound '\
                         'to any solute atoms.')
# Optional arguments
parser.add_argument('--solvent_sel', default='resname XXXX and name H*',
                    help='selection for solvent atoms. Must be the same '\
                         'than the selection used to create the flags '\
                         'trajectory file.')
parser.add_argument('--qvalues', type=str,
                    default='"0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9"',
                    help='Triad defining the momentum transfer values '\
                         'min_val, max_val, step. For instance, '\
                         '0.1, 2.1, 0.2 will yield Q values from 0.1 up '\
                         'to and included 1.9. Last value 2.1 is excluded.')
parser.add_argument('--dt', type=int,
                    default=10,
                    help='Elemental elapsed time, t of I(Q,t), in units '\
                         'of number of trajectory frames')
parser.add_argument('--nt', type=int, default=1000,
                    help='Number of elapsed times')
parser.add_argument('--prefix', default='',
                    help='prefix each output ISF with this input sting')
argcomplete.autocomplete(parser)
args = parser.parse_args()

a_universe = mda.Universe(args.topology_file, args.trajectory_file)

# Find when are solvent (water Hydrogen) atoms bound to solute (polymer)
print('Finding solvent atoms bound to solute')
flags = np.load(args.save_flags_file)

# Extract coordinates of solvent Hydrogen atoms
print('Extract coordinates of solvent atoms')
solvent_group = a_universe.select_atoms(args.solvent_sel)
positions = list()
for _ in tqdm(a_universe.trajectory, total=len(a_universe.trajectory)):
    positions.append(solvent_group.positions)

# Find scattering
b = np.ones(len(solvent_group))  # scattering lengths set to unity
q_values = np.array([float(q) for q in args.qvalues.strip('"').split(',')])
# frames separated every 0.1ps
# nt = 1000  number of time points
# dt = 10 separation between time consecutive time points (1ps)
print('Calculate intermediate incoherent structure factors')
sqe = scatter.II(positions, b, q_values, dt=args.dt, nt=args.nt, c1=flags)

for scatter_function in 'sf sfNN sfNY sfYN sfYY'.split():
    file_name = '{}{}.h5'.format(args.prefix, scatter_function)
    scatter.saveIISassenaFormat(sqe[scatter_function], q_values, file_name)

