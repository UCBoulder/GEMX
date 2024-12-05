#pragma once

#include <complex>
#include <cmath>

#include "MultiArrays.hpp"

extern int numprocs;
extern int cnt;
extern int mmx;
extern int iseed;
extern int MyId;
extern int imx;
extern int jmx;
extern int kmx;
extern int idg;
extern int ierr;

extern int nx;
extern int nz;

extern double dxeq;
extern double xdim;
extern double xctr;

extern double dzeq;
extern double zdim;

extern double dx;
extern double dz;
extern double dzeta;
extern double dt;
extern double totvol;
extern double n0;
extern double tcurr;

extern double bu,tu,nu,xu,frequ,vu,eru;

//1D arrays
extern int *tmm_ptr; //c_loc will grab address in fortran (of given index) and will set a pointer to point at it.
extern int *mm_ptr;

extern double *zeta2_ptr;
extern double *x2_ptr;
extern double *z2_ptr;
extern double *mims_ptr;
extern double *u2_ptr;
extern double *mu_ptr;
extern double *w2_ptr;

extern double *x3_ptr;
extern double *z3_ptr;
extern double *u3_ptr;
extern double *w3_ptr;
extern double *zeta3_ptr;

//extern double *Rgrid_ptr;

//2D Arrays
extern double *b0_ptr;
extern Array2D<double> b0_c; //= Array2D<double>(nullptr, 0,0);

extern double *t0i_ptr;
extern Array2D<double> t0i_c;//= Array2D<double>(nullptr, 0,0);

extern double *xn0i_ptr;
extern Array2D<double> xn0i_c;//= Array2D<double>(nullptr, 0,0);

extern double *ileft_ptr;
extern Array2D<double> ileft_c;

extern double *xbackw_ptr;
extern Array2D<double> xbackw_c;

extern double *zbackw_ptr;
extern Array2D<double> zbackw_c;

extern double *jleft_ptr;
extern Array2D<double> jleft_c;

extern double *iright_ptr;
extern Array2D<double> iright_c;

extern double *xforw_ptr;
extern Array2D<double> xforw_c;

extern double *jright_ptr;
extern Array2D<double> jright_c;

extern double *zforw_ptr;
extern Array2D<double> zforw_c;

extern double *b0zeta_ptr;
extern Array2D<double> b0zeta_c;

//Array c++ translate
extern "C"
{
    void new_gemx_com_c_();
}

void Allocate2dPointerArrays_gemx_com();
void Allocate3dPointerArrays_gemx_com();
