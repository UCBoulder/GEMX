#ifndef EQUIL_EXTERNS
#define EQUIL_EXTERNS

#include <complex>
#include <cmath>

#include "MultiArrays.hpp"

extern int nx;
extern int nz;
extern int nr, nr2, ntheta, isgnf, isgnq, isupae, tor_n;
extern int num_lines;
extern int lines;

extern double dxeq;
extern double xdim;
extern double xctr;
extern double dzeq;
extern double zdim;
extern double bu;
extern double tu;
extern double nu;
extern double xu;
extern double frequ;
extern double vu;
extern double eru;
extern double dR,dth,mu0,e,proton;


//1D arrays
extern double *psitab_ptr;
extern double *weight00_ptr;
extern double *weight10_ptr;
extern double *weight01_ptr;
extern double *weight11_ptr;
extern double *jacobian_ptr;
extern double *deno_ptr;

extern int *gindex_ptr;
extern int *iarray_ptr;
extern int *jarray_ptr;
extern int *priv_ptr;

//extern double *Rgrid_ptr;
//extern double *Zgrid_ptr;

//2D arrays
extern double *b0_ptr;
extern Array2D<double> b0_c; 

extern double *t0i_ptr;
extern Array2D<double> t0i_c;

extern double *xn0i_ptr;
extern Array2D<double> xn0i_c;

extern double *b0zeta_ptr;
extern Array2D<double> b0zeta_c;

extern double *phiavg_ptr;
extern Array2D<double> phiavg_c;

extern double *t0s_ptr;
extern Array2D<double> t0s_c;

extern double *xn0s_ptr;
extern Array2D<double> xn0s_c;

extern double *capts_ptr;
extern Array2D<double> capts_c;

extern double *capns_ptr;
extern Array2D<double> capns_c;

extern double *vpars_ptr;
extern Array2D<double> vpars_c;

extern double *vparsp_ptr;
extern Array2D<double> vparsp_c;

extern double *psi_p_ptr;
extern Array2D<double> psi_p_c;

extern double *mask_ptr;
extern Array2D<double> mask_c;

extern double *mask2_ptr;
extern Array2D<double> mask2_c;

extern double *mask3_ptr;
extern Array2D<double> mask3_c;

extern double *mask4_ptr;
extern Array2D<double> mask4_c;

// extern double *c2_over_vA2_ptr; //Currently Broken, not sure why
// extern Array2D<double> c2_over_vA2_c;

extern double *xn0e_ptr;
extern Array2D<double> xn0e_c;

extern double *t0e_ptr;
extern Array2D<double> t0e_c;

extern "C"
{
    void new_gemx_com_c_();
}

#endif //EQUIL_EXTERNS