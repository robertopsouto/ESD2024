#!/bin/bash

#SBATCH --nodes=1                       # here the number of nodes
#SBATCH --ntasks=1                      # here total number of mpi tasks
#SBATCH --cpus-per-task=1               # number of cores per node
#SBATCH -p cpu_dev                      # target partition
#SBATCH --threads-per-core=1
#SBATCH -J NPB_BT-MZ                    # job name
#SBATCH --time=00:05:00                 # time limit
#SBATCH --exclusive                     # to have exclusive use of your nodes

echo "Cluster configuration:"
echo "==="
echo "Partition: " $SLURM_JOB_PARTITION
echo "Number of nodes: " $SLURM_NNODES
echo "Number of MPI processes: " $SLURM_NTASKS " (" $SLURM_NNODES " nodes)"
echo "Number of MPI processes per node: " $SLURM_NTASKS_PER_NODE
echo "Number of threads per MPI process: " $SLURM_CPUS_PER_TASK
echo "NPB Benchmark: " $1
echo "Bechmark class problem: " $2

###################################
#           COMPILER              #
###################################
module load openmpi/gnu/4.1.2+cuda-11.2

DIR=$PWD

bench=${1}
class=${2}
execfile="${bench}.${class}.x"
BIN=$DIR/${execfile}

export OMPI_MCA_opal_warn_on_missing_libcuda=0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd $DIR

srun -n $SLURM_NTASKS $BIN

dirdest="${bench}_${class}_MPI-${SLURM_NTASKS}_OMP-${SLURM_CPUS_PER_TASK}_JOBID-${SLURM_JOBID}"
mkdir $dirdest
cp slurm-${SLURM_JOBID}.out $dirdest/
