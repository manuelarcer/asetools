#!/bin/bash

#PBS -N Battery-NEB
#PBS -l select=1:ncpus=128:mem=480G
#PBS -l walltime=2:00:00
#PBS -j oe
#PBS -o alpha.out
#PBS -e alpha.err
#PBS -P 13003830
#PBS -q normal


cd $PBS_O_WORKDIR
echo Working directory is $PBS_O_WORKDIR

## Especific to ASPIRE2A
module swap PrgEnv-cray PrgEnv-intel/8.3.3
module list

# Calculate the number of processors allocated to this run.
export NPROCS=`wc -l < $PBS_NODEFILE`

# Calculate the number of nodes allocated.
export NNODES=`uniq $PBS_NODEFILE | wc -l`

### Display the job context
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Using ${NPROCS} processors across ${NNODES} nodes

######################################
## From Mike's submission script #####
export OMP_NUM_THREADS=1
myppn=128
mpitasks=$(($NNODES * $myppn / $OMP_NUM_THREADS))   ## Mike's
#mpitasks=$(($NPROCS / $OMP_NUM_THREADS))             ## Mine
echo "mpitasks is $mpitasks"
export MKL_DEBUG_CPU_TYPE=5
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=off
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_STACKSIZE=2G
ulimit -s unlimited
export OMP_WAIT_POLICY=PASSIVE
######################################

## Setting ASE env variables
export VASPDIR='/home/project/13001717/software/vasp.6.3.2-intel-libxc6/bin'
export VASP_EXEC='vasp_std'
export VASP_COMMAND='aprun -n '$mpitasks' -N '$myppn' -d '$OMP_NUM_THREADS' '$VASPDIR'/'$VASP_EXEC' | tee vasp.out >/dev/null'
echo $VASP_COMMAND

echo " "
python calcVASPcalc_battery.py