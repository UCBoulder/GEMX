#Shell script to load perlmutter environment.
module restore
module load nvhpc
export PETSC_PATH=/global/cfs/cdirs/mp118/software/petsc_3.19.3/install
export LD_LIBRARY_PATH=$PETSC_PATH/lib:$LD_LIBRARY_PATH
