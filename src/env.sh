#!/bin/bash
# # Source the system default environment
# source /opt/cray/pe/cpe/24.07/restore_lmod_system_defaults.sh

## Load required modules in the right order
module restore
module load PrgEnv-nvidia
module load cudatoolkit
module load cray-mpich/8.1.25
module load cray-fftw
#
# # Set PETSC path if needed
export PETSC_PATH=/global/cfs/cdirs/mp118/software/petsc/install
export MPIPATH=$CRAY_MPICH_DIR
export LD_LIBRARY_PATH=/global/cfs/cdirs/mp118/software/petsc/install/lib:$LD_LIBRARY_PATH
