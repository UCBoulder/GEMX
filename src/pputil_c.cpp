#include <iostream>
#include <iomanip>

#include "pputil_c.hpp"

using namespace std;

static int me, nvp, npp, GCLR, TCLR, p_color, p_rank;
// static MPI_Comm GRID_COMM, TUBE_COMM, PETSC_COMM;

void ppinit_c(int& idproc, int& nproc, int &ntube,int &imx, int& i3D, MPI_Comm &com1,MPI_Comm &com2, MPI_Comm &com_petsc,int &petsc_color,int &petsc_rank){
    int ierr;
    int n_tor;
    int n_tor_planes;
    
    ierr = MPI_Init(nullptr, nullptr); //Null since Init in c++ can interpret command line input, there are none here
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &npp);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &me);
    nproc = npp;
    idproc = me;

    GCLR = static_cast<int>(me/ntube);
    TCLR = me%ntube;

    ierr = MPI_Comm_split(MPI_COMM_WORLD, GCLR, TCLR, &GRID_COMM);
    ierr = MPI_Comm_split(MPI_COMM_WORLD, TCLR, GCLR, &TUBE_COMM);

    if(i3D != 0){
        n_tor_planes = kmx+1;
    } else {
        n_tor_planes = 1;
    }
    
    p_color = static_cast<int>(me*(n_tor_planes)/npp);
    p_rank = me%(npp/(kmx+1));
    ierr = MPI_Comm_split(MPI_COMM_WORLD, p_color, p_rank, &PETSC_COMM);
    petsc_color = p_color;
    petsc_rank = p_rank;
//         else
//            CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,0,nproc,PETSC_COMM,ierr)
// 
    com1 = TUBE_COMM;
        com2 = GRID_COMM;
        com_petsc = PETSC_COMM;
    nvp = npp/ntube;
}