#include "gemx_com_c.h"

//Define global pointers 1D

//Define global pointers 2D

//Define global pointers 3D

Array2D<double> b0_c;
Array2D<double> t0i_c; 
Array2D<double> xn0i_c; 

void new_gemx_com_c_(){
    //call to functions that will calculate location of data in fortran arrays
    Allocate2dPointerArrays_gemx_com();
}

//Allocation of 2D and 3D arrays below
void Allocate2dPointerArrays_gemx_com(){
    b0_c.CreateArray2D(b0_ptr,nx,nz);
    t0i_c.CreateArray2D(t0i_ptr,nx,nz);
    xn0i_c.CreateArray2D(xn0i_ptr, nx, nz);
}