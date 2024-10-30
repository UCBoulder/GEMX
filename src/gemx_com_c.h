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

//1D arrays
extern int *tmm_ptr;
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

//2D Arrays
extern double *b0_ptr;
extern Array2D<double> b0_c; //= Array2D<double>(nullptr, 0,0);

extern double *t0i_ptr;
extern Array2D<double> t0i_c;//= Array2D<double>(nullptr, 0,0);

extern double *xn0i_ptr;
extern Array2D<double> xn0i_c;//= Array2D<double>(nullptr, 0,0);


//Array c++ translate
extern "C"
{
    void new_gemx_com_c_();
}

void Allocate2dPointerArrays_gemx_com();
