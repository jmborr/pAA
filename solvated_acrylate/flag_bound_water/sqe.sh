#!/bin/bash

# IRIS values for momentum transfer
qvalues="0.483619,0.607871,0.729147,0.847027,0.960857,1.07024,1.17452,1.2732,1.36594,1.45216,1.53158,1.60372,1.66829,1.72494,1.77349,1.8136,1.8451"

# Short runs, unwrapped, 2NS, 10000 frames or 0.2ps/frame
for name in PAA.T270K.2NS  PAA.T280K.2NS  PAA.T290K.2NS  PAA.T300K.2NS  PAA.T310K.2NS;do
    # Output 2000 time points separated by 5frame * 0.2ps/frame = 1ps
    python sqe.py acrylate.pdb ${name}_unwrapped.dcd  ${name}.npy --qvalues=$qvalues --prefix="${name}_" --dt=5 --nt=2000 &
    sleep 60s
done
