#include "gemx_com_c.h"
#include "equil_c.h"
#include "mpi.h"

void ppinit_c(int& idproc, int& nproc, int &ntube,int &imx, int& i3D, MPI_Comm &com1,MPI_Comm &com2, MPI_Comm &com_petsc,int &petsc_color,int &petsc_rank);