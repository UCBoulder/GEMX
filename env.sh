#!/bin/bash

module reset

# Pre-maintenance Cray Programming Environment
module load cpe/24.07

# IMPORTANT: explicitly use the PrgEnv version belonging to CPE 24.07
module load PrgEnv-nvidia/8.5.0

# Match the toolchain used to build PETSc 3.22.3
module load nvidia/24.5
module load cudatoolkit/12.4
module load cray-mpich/8.1.30

# Application dependency
module load cray-fftw

export PETSC_PATH=/global/cfs/cdirs/mp118/software/petsc/install_02242025
export PETSC_DIR=$PETSC_PATH

export MPIPATH=$CRAY_MPICH_DIR

# export LD_LIBRARY_PATH=$PETSC_PATH/lib:${LD_LIBRARY_PATH:-}
export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$PETSC_PATH/lib:${LD_LIBRARY_PATH:-}