#!/bin/bash

#!/bin/bash -l
#SBATCH --account=m1503
#SBATCH --partition=debug
#SBATCH --nodes=32
#SBATCH --time=00:30:00
#SBATCH --job-name=benchmark
#SBATCH --output benchmark.o%j
#SBATCH --license=SCRATCH

prefix='relax'
currindex=1

# Load machine-dependent lammps
source $MODULESHOME/init/bash
module load lammps/20161117
lmp=lmp_edison
echo "LAMMPS executable is $lmp"

cd $SLURM_SUBMIT_DIR
srun -n 768  $lmp -in benchmark_NVT.in # Do not use "<" in place of "-in"
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
