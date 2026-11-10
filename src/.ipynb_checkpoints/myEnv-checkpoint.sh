#!/bin/bash

# Reset to NERSC defaults
module reset

# Remove MATLAB to prevent its libstdc++ from contaminating GEMX
module unload matlab 2>/dev/null || true

# GPU programming environment
module load PrgEnv-nvidia
module load cudatoolkit/13.2
module load cray-mpich
module load cray-fftw

# PETSc 3.25.5 + CUDA 13
export PETSC_PATH=/global/cfs/cdirs/mp118/software/petsc/install_3.25.5_cuda13
export PETSC_DIR=$PETSC_PATH
export MPIPATH=$CRAY_MPICH_DIR

export LD_LIBRARY_PATH=$PETSC_PATH/lib:${LD_LIBRARY_PATH:-}