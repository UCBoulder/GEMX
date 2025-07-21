#Shell script to load perlmutter environment.
# module restore
# ##module load nvhpc
# module load cray-fftw
# module load PrgEnv-nvidia
# module load cuda
# export PETSC_PATH=/global/cfs/cdirs/mp118/software/petsc/install
# export LD_LIBRARY_PATH=$PETSC_PATH/lib:$LD_LIBRARY_PATH


module restore
module load nvidia
module load cray-fftw
# export LD_LIBRARY_PATH=/global/cfs/cdirs/mp118/software/petsc/install/lib:$LD_LIBRARY_PATH

export PETSC_PATH=/global/cfs/cdirs/mp118/software/petsc/install
# export FFT_PATH=/global/homes/u/u10198/installed/dfftpack_cray/libdfftpack.a #Calder Edit
export LD_LIBRARY_PATH=$PETSC_PATH/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/global/cfs/cdirs/mp118/software/petsc/install_02242025/lib:$LD_LIBRARY_PATH
~                                              