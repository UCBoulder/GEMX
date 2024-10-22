# GEMX

Gyrokinetic simulation of the tokamak edge including the scrape-off layer

Please add any useful notes to this file. Place it anywhere and we can clean up later.

GEMX.pdf is our equations or user manual.  Consider placing any notes on implementations in this overleaf document: <https://www.overleaf.com/project/64456a5749a720661d4d1f0a>

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
cd ..
./newRun.py
cd runs/run001
sbatch job_gpu
```

Mostly files related to the build are copied to the ```bin``` folder, including ```gemx.in``` and the ```job``` associated with the make configuration. This folder is meant to just store some of these
build related files, but not meant to be run out of because recompiling will overwrite changes to the job/input files for instance. To get around this issue the ```newRun.py``` script is provided to make run directories which are stored in ```runs``` and not source controlled. Run directories should contain their own versions of job/input files associated with the run, so that originals remain
unchanged.

The ```newRun.py``` script generates new run folders in runs: ```run001```, ```run002```, and so on. You can also pass a run folder name with ```-n``` flag using the run script. Everything from ```bin``` is copied to the run directory except ```gemx``` and ```codeChanges.txt```. The executable ```gemx``` is remade as a symbolic link to the original in ```bin```. That way when recompiling, all run directories will point to the latest version of the code.

The file ```codeChanges.txt``` will be copied over by the job scripts at run time. This file is generated at compile time and stores the most recent info about the code (local changes to source code, git commit, branch, and compile info via the job script name in ```bin``` (gpu/cpu/debug/etc.)) That way old runs can be recreated fairly easily.

All job scripts for the necessary system are copied over to the run directory so changes to them are stored in the run folder. They currently have different names for clarity.

The shell scripts are symbolic links to the ```tools``` folder similar to ```gemx```, since when changes should be made to these they likely should be source controlled and updated for every run.

To clean up run output, use the reset script in the run directory: ```./reset.sh```. ```env.sh``` is for loading the perlmutter environment, and is also provided in the ```src``` directory.

The code will run with GPU acceleration by default, but can be run cpu only as well as mentioned below. It can be run with debug capablities as well as described below.

Again, run-specific files like ```gemx.in``` or job scripts are meant to be changed in the run directory to not affect the source-controlled versions. If changes should be permanent, remember to change the original files as well when pushing. Then probably a new run folder should be made as the old ones are more or less obsolete, but the updated run-specific files could be copied manually to old run directories if preferred.

### Running on CPU

To run serially on cpu use ```make CPU=1```. It's possible to run ```OpenACC``` parallel on CPU with the ```-acc=multicore``` flag rather than ```-acc=host``` but this might not work depending on the CPUs and probably using ```OpenMP``` would be more prudent on CPU.

### Debugging

Create a debug build by calling ```make DEBUG=1```. This will lower the optimization and add check flags.

Also the debug job script will dump core files as well if needed on perlmutter. The core dump files can be read with ```gdb gemx <core_file>```. The stack trace can be looked at using bt or backtrace; the stack trace doesn't seem to be output in run.err even with the current compiler settings.

Note, bounds checks are disabled by OpenACC. A CPU run would be required, or a tool like valgrind could be used to bound check CPU code when OpenACC is using GPU.

There is a ```dbg``` flag in ```gemx.in``` which can be set to ```1``` to enter an infinite while loop. A file ```launch.json``` is included for VSCode users to attach to the process and debug at real time. Set ```dbg=0``` in the debug console once attached to continue. This is meant for local runs, not perlmutter, though it can really work anywhere as long as it is run on only 1 MPI process, otherwise all others will be stuck in an infinite while loop and MPI calls won't work when they need to sync.

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

We also need to set an MPI path in ```~/.bashrc``` to include the header files when compiling (this path could be defined earlier and used above as well):

```bash
export MPIPATH=$NVCOMPILERS/$NVARCH/23.9/comm_libs/mpi
```

Then download PETSc to your preferred location (currently v3.19.3 is installed manually on NERSC for mp118) per the [installation instructions](https://petsc.org/release/install/):

```bash
git clone -b release https://gitlab.com/petsc/petsc.git petsc
cd petsc
git checkout v3.19.3
```

Update ```~/.bashrc``` with the following (and set your cloned path) for PETSc as well:

```bash
export PETSC_PATH=<cloned-petsc-location>/installation
export LD_LIBRARY_PATH=$PETSC_PATH/lib:$LD_LIBRARY_PATH
```

Reload the terminal environment so PETSc and GEMX can be compiled: ```source ~/.bashrc```.

Then configure PETSc to install to the chosen path and for an optimized build:

```bash
./configure --prefix=$PETSC_PATH --with-debugging=0 --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 COPTFLAGS='-O3' CXXOPTFLAGS='-O3' FOPTFLAGS='-O3'
```

This should set PETSc to be compiled with the Nvidia compilers given the Nvidia environment paths above. Then just follow the command line instructions output by PETSc to further make and install the library.

Then one can make and run GEMX using the instructions from before, but without needing ```source env.sh``` or ```sbatch```:

```bash
git clone git@github.com:UCBoulder/GEMX.git
cd src
make
cd ..
./newRun.py
cd runs/run001
./job_local
```

## Analysing Runs

Scripts to read run output are also copied to the ```run``` folders. Some are MATLAB scripts some are Jupyter notebook files at the moment. They are stored in the ```analysis``` folder.

MATLAB scripts ```readden.m```, ```readphi.m``` and ```phi_gif.m``` can be opened by MATLAB. ```readden.m``` will plot the ion's density and n_i*v_parallel.   ```readphi.m``` will plot the perturbed fields, perturbed electron density and parallel current. ```phi_gif.m``` will plot a gif figure for the perturbed electrostatic field. Open the file, change the MATLAB working path to your running folder, run the script, select ```add to path```.

Python script ```readphi.ipynb``` should be placed in the running folder to be run. It will plot the traced particles' orbits.

These could also be made symbolic links to the originals, so they can be updated for all runs, but then personal changes need to be worked around somehow (unique data paths, data plotted in jupyter notebooks) or those specific changes should be removed before updating source controlled files. This is fine with python/matlab though when using a GUI/IDE to set relative run directories for instance
and output plots in said GUI.

## Cheat Sheet

tar cvf * - When taring up GEMX, do not include top directory, just the files and any subdirectories.

tar xvf file.tar

NERSC user names:
cqh2021 - Qiheng
u10198 - Yang
nonohing - Zhichen
jycheng - Junyi
stirkas - Stefan

give -u \<username\> file.tar
take -u \<username\> -a

source env.sh

sbatch submit.cmd
sacct
sqs
scancel $jobid
scancel -u $USER
