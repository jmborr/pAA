#!/bin/bash

#!/bin/bash -l
#SBATCH --repo=m1503
#SBATCH --partition=regular
#SBATCH --nodes=8
#SBATCH --time=04:00:00
#SBATCH --job-name=r1
#SBATCH --output r1.o%j
#SBATCH --license=SCRATCH

prefix='relax'
currindex=1

# Load machine-dependent lammps
source $MODULESHOME/init/bash
module load lammps/
lmp=lmp_edison
echo "LAMMPS executable is $lmp"

cd $PBS_O_WORKDIR
srun -np 196  $lmp -in $prefix.$currindex.in # Do not use "<" in place of "-in"
#wait  # only needed when submitting several executables

###############################
# Some useful tidbits:
#  Qeue usage info at https://www.nersc.gov/users/queues/queue-wait-times/
#     qsub -W depend=afterok:123451 2.pbs  chain submissions
#     showbf   shows unallocated Nodes

# Setting number of cores in each cluster, running parallel jobs
#        Edison                Hopper                Carver
#  PBS -l mppwidth=24     PBS -l mppwidth=24    PBS -l nodes=3:ppn=8
#     aprun -np 24          aprun -np 24           mpirun -np 24
