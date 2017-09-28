#!/bin/bash

datadir="/data/nfs/camm/users/ggz"
# Short runs, wrapped, 2NS, 10000 frames or 0.2ps/frame
for dcd in PAA.T270K.2NS.dcd  PAA.T280K.2NS.dcd  PAA.T290K.2NS.dcd  PAA.T300K.2NS.dcd  PAA.T310K.2NS.dcd;do
    base="${dcd%.*}"  # remove extension from name
    python flag_atoms.py acrylate.pdb $datadir/$dcd ${base}.npy &
done
