#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import numpy as np
import argparse
import argcomplete
import MDAnalysis as mda
from MDAnalysis.analysis.distances import dist as mda_dist
from tqdm import tqdm
import concurrent.futures as futures

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Flag trajectory for surface water hydrogen atoms')
    parser.add_argument('topology', type=str, help='Topology file (PDB is OK)')
    parser.add_argument('trajectory', type=str, help='Atomic trajectory')
    parser.add_argument('solute', type=str,
                        help='selection for solute atoms, in MDAnalysis syntax')
    parser.add_argument('hydrogens', type=str,
                        help='selection for hydrogen atoms, in MDAnalysis syntax')
    parser.add_argument('flag_traj', type=str,
                        help='name of the output flag trajectory file')
    parser.add_argument('cutoff', type=float, default=4.0,
                        help='maximum distance for bound water hydrogens to any solute atom')
    argcomplete.autocomplete(parser)
    args = parser.parse_args()
    u = mda.Universe(args.topology, args.trajectory)

