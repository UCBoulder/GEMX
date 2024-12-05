#include "gemx_com_c.h"

//Define global pointers 1D

//Define global pointers 2D

//Define global pointers 3D

//other pointers
Array2D<double> b0_c;
Array2D<double> t0i_c; 
Array2D<double> xn0i_c; 
Array2D<double> ileft_c;
Array2D<double> xbackw_c;
Array2D<double> zbackw_c;
Array2D<double> jleft_c;
Array2D<double> iright_c;
Array2D<double> xforw_c;
Array2D<double> jright_c;
Array2D<double> zforw_c;
Array2D<double> b0zeta_c;

void new_gemx_com_c_(){
    //call to functions that will calculate location of data in fortran arrays
    Allocate2dPointerArrays_gemx_com();
    Allocate3dPointerArrays_gemx_com();
}

//Allocation of 2D and 3D arrays below
void Allocate2dPointerArrays_gemx_com(){
    b0_c.CreateArray2D(b0_ptr,nx,nz);
    t0i_c.CreateArray2D(t0i_ptr,nx,nz);
    xn0i_c.CreateArray2D(xn0i_ptr, nx, nz);
    ileft_c.CreateArray2D(ileft_ptr,imx,jmx);
    xbackw_c.CreateArray2D(xbackw_ptr,imx,jmx);
    zbackw_c.CreateArray2D(zbackw_ptr,imx,jmx);
    jleft_c.CreateArray2D(jleft_ptr,imx,jmx);
    iright_c.CreateArray2D(iright_ptr,imx,jmx);
    xforw_c.CreateArray2D(xforw_ptr,imx,jmx);
    jright_c.CreateArray2D(jright_ptr,imx,jmx);
    zforw_c.CreateArray2D(zforw_ptr,imx,jmx);
    b0zeta_c.CreateArray2D(b0zeta_ptr, nx, nz);
}

void Allocate3dPointerArrays_gemx_com(){
    
}