#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

from __future__ import print_function

import numpy as np
import os
import sys
import argparse
import argcomplete
#import concurrent.futures as futures

import MDAnalysis as mda
import myMDAnalysis.ContactMapAnalysis.ContactMapAnalysisAPI as cma

#  Required arguments
description="""Calculate the flags trajectory determining
for each solvent atom whether it is bound to any solute atom"""
fmt = argparse.ArgumentDefaultsHelpFormatter
parser = argparse.ArgumentParser(description=description,
                                 formatter_class=fmt)
parser.add_argument('topology_file', type=str,
                    help='topology file. PDB format is OK')
parser.add_argument('trajectory_file', type=str,
                    help='Trajectory file, wrapped coordinates')
parser.add_argument('save_flags_file', default='bound_solute_flags.npy',
                    help='File name to save the flags trajectory.')
#  Optional arguments
parser.add_argument('--solvent_sel', default='resname XXXX and name H*',
                    help='selection for solvent atoms')
parser.add_argument('--solute_sel',
                    # Default solute is polymer plus bound Na atoms
                    default='(not resname XXXX) or (resname XXXX and name Na*)',
                    help='selection for solute atoms')
parser.add_argument('--cut_off', default=2.5,
                    help='maximum distance from solvent atom to solute '\
                         'to consider solute atom as bound')
argcomplete.autocomplete(parser)
args = parser.parse_args()

print('Loading trajectory into MDAnalysis universe...')
a_universe = mda.Universe(args.topology_file, args.trajectory_file)
print('Calculate flags trajectory')
flags = cma.solvent_bound_flag(a_universe, args.solvent_sel,
                               args.solute_sel, args.cut_off)
flags = np.array(flags).transpose()  # shape = (nframes, natoms)
print('Save to file')
np.save(args.save_flags_file, flags)
