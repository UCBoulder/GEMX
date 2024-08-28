# GEMX

Gyrokinetic simulation of the tokamak edge including the scrape-off layer

Please add any useful notes to this file. Place it anywhere and we can clean up later.

GEMX.pdf is our equations or user manual.  Consider placing any notes on implementations in this overleaf document: https://www.overleaf.com/project/64456a5749a720661d4d1f0a

Versions:

Sept 13, 2023 
Have GK Poisson working and particle tracers in equilibrium B working.  Zhichen is working getting Alfven waves.  Field solve needs to be parallelized and optimized. Qiheng and Scott are cleaning up wrappers, adding diagnostic routines and .ipynb diagnostic scripts.

## Setup and Running

On perlmutter:

```bash
git clone git@github.com:UCBoulder/GEMX.git
cd src
source env.sh
make
cd ../bin
sbatch <job_script_name>
```

All files necessary to run are copied to the ```bin``` folder. However this is done everytime ```make``` is called so changes to the input file ```gemx.in``` or NERSC job scripts may get overwritten
if compiling code updates.

A temporary fix is to make a new run directory and copy files there from ```bin```. Then update the job script to call ```../bin/gemx``` rather than ```./gemx```. Then you can have specific run files which call the latest version of the code.
${\color{red}\textbf{\textrm{(This should be automated with a script.)}}}$

To clean up run output, use the reset script in the run directory: ```./reset.sh```

Different job scripts will be moved depending on make flags. They currently have different names for clarity.

The code can be run using cpu or gpu. As well as in debug mode as described below.

### Running on CPU

To run serially on cpu use ```make CPU=1```. It's possible to run ```OpenACC``` parallel on CPU with the ```-acc=multicore``` flag rather than ```-acc=host``` but this might not work depending on the CPUs and probably using ```OpenMP``` would be more prudent on CPU.

### Debugging

Create a debug build by calling ```make DEBUG=1```. This will lower the optimization and add check flags.

Also the debug job script will dump core files as well if needed on perlmutter. The core dump files can be read with ```gdb gemx <core_file>```.

Note, bounds checks are disabled by OpenACC. A cpu run would be required.

### Running Locally

The makefile is setup to run locally and on perlmutter without changing any settings. The only difference is the ```env.sh``` doesn't need to be called locally. Instead one should update their local environment in ```~/.bashrc``` manually.

To run locally, the Nvidia HPC compilers and PETSc need to be installed.

At the time of writing this, perlmutter uses version 23.9 of the Nvidia HPC compilers available [here](https://developer.nvidia.com/nvidia-hpc-sdk-releases).

Update ```~/.bashrc``` with the following, per the [installation guide](https://docs.nvidia.com/hpc-sdk/archive/23.9/hpc-sdk-install-guide/index.html):
```bash
NVARCH=`uname -s`_`uname -m`; export NVARCH
NVCOMPILERS=/opt/nvidia/hpc_sdk; export NVCOMPILERS
MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/23.9/compilers/man; export MANPATH
PATH=$NVCOMPILERS/$NVARCH/23.9/compilers/bin:$PATH; export PATH
export PATH=$NVCOMPILERS/$NVARCH/23.9/comm_libs/mpi/bin:$PATH
export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/23.9/comm_libs/mpi/man
```

Then download PETSc to your preferred location (currently v3.19.3 is installed manually on NERSC for mp118) per the [installation instructions](https://petsc.org/release/install/):
```bash
git clone -b release https://gitlab.com/petsc/petsc.git petsc
cd petsc
git checkout v3.19.3
```

Update ```~/.bashrc``` with the following for PETSc as well:
```bash
export PETSC_PATH=<cloned-petsc-location>/install
export LD_LIBRARY_PATH=$PETSC_PATH/lib:$LD_LIBRARY_PATH
```

Reload the terminal environment so PETSc and GEMX can be compiled: ```source ~/.bashrc```.

Then configure PETSc to install to the chosen path:
```bash
./configure --prefix=$PETSC_PATH/install --with-debugging=0 --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 COPTFLAGS='-O3' CXXOPTFLAGS='-O3' FOPTFLAGS='-O3'
```
This should set PETSc to be compiled with the Nvidia compilers given the Nvidia environment paths above.

Finally, follow the instructions output by PETSc to further make and install the library.

Then one can make and run GEMX using the instructions from before, but without needing ```source env.sh``` or ```sbatch```:
```bash
git clone git@github.com:UCBoulder/GEMX.git
cd src
make
cd ../bin
./job_local
```

## Analysing Runs

Scripts to read run output are also copied to the ```bin``` folder. Some are MATLAB scripts some are Jupyter notebook files. 

MATLAB scripts ```readden.m```, ```readphi.m``` and ```phi_gif.m``` can be opened by MATLAB. ```readden.m``` will plot the ion's density and n_i*v_parallel  .   ```readphi.m``` will plot the perturbed fields, perturbed electron density and parallel current. ```phi_gif.m``` will plot a gif figure for the perturbed electrostatic field. Open the file, change the MATLAB working path to your running folder, run the script, select ```add to path```.

Python script ```readphi.ipynb``` should be placed in the running folder to be run. It will plot the traced particles' orbits.

## Cheat Sheet

tar cvf * - When taring up GEMX, do not include top directory, just the files and any subdirectories.

tar xvf file.tar

NERSC user names:
cqh2021 - Qiheng
u10198 - Yang
nonohing - Zhichen
jycheng - Junyi
stirkas - Stefan

give -u <username> file.tar
take -u <username> -a

source env.sh

sbatch submit.cmd
sacct
sqs
scancel $jobid
scancel -u $USER
