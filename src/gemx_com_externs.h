#pragma once

#include <complex>
#include <cmath>

#include "MultiArrays.hpp"


//Array c++ translate
extern "C"
{
    void new_gemx_com_c_();
}

extern int numprocs;
extern int Last, myid, cnt , ierr;
extern int mmx;
extern int iseed;
extern const int imx;
extern const int jmx;
extern const int kmx;
extern int idg;
extern int nm,nsm,ncurr,iflr,ifield_solver,ntracer,i3d,icollision;
extern int timestep, iez;
extern int nonlin,nonline,iflut,ifluid,ipara;
extern int iput,iget,ision,isham,peritr,iadi;
extern int nmx,nsmx,nsubd,ntube,petsc_color,petsc_rank,iBoltzmann,globle_integer,eBoltzmann,eAdiabatic,iterations,dbg;
extern int lx,lz;
extern int nplot,xnplt;
//extern MPI_Comm GRID_COMM,TUBE_COMM, PETSC_COMM;

extern double dx;
extern double dz;
extern double dzeta;
extern double dt;
extern double totvol;
extern double n0;
extern double tcurr;
extern double bu,tu,nu,xu,frequ,vu,eru;
extern double vcut;
extern double starttm,lasttm,tottm;
extern double start_total_tm, end_total_tm, start_integ_tm, end_integ_tm, start_ppush_tm, end_ppush_tm, start_cpush_tm, end_cpush_tm;
extern double total_tm, integ_tm, ppush_tm, cpush_tm;
extern double cut,amp,tor,amie,emass,qel,rneu;


//1D arrays
extern int *tmm_ptr;
extern int *mm_ptr;
extern int *rand_table_ptr;
extern double *zeta2_ptr;

//extern double *x2_ptr;
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
extern double *q_ptr;
extern int *lr_ptr;
extern double *jac_ptr;
extern double *xg_ptr;

//2D Arrays
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

extern double *den2d2_ptr;
extern Array2D<double> den2d2_c;

extern double *dden2d_ptr;
extern Array2D<double> dden2d_c;

extern double *den2d1_ptr;
extern Array2D<double> den2d1_c;

extern double *gn0e_ptr;
extern Array2D<double> gn0e_c;

extern double *gbtor_ptr;
extern Array2D<double> gbtor_c;

extern double *bmag_ptr;
extern Array2D<double> bmag_c;

extern double *gcpnex_ptr;
extern Array2D<double> gcpnex_c;

extern double *gcpnez_ptr;
extern Array2D<double> gcpnez_c;

//3D arrays
extern double *phi_ptr;
extern Array3D<double> phi_c;

extern double *ex_ptr;
extern Array3D<double> ex_c;

extern double *ez_ptr;
extern Array3D<double> ez_c;

extern double *ezeta_ptr;
extern Array3D<double> ezeta_c;

extern double *phi_k_ptr;
extern Array3D<double> phi_k_c;

extern double *dphidr_ptr;
extern Array3D<double> dphidr_c;

extern double *dphi_kdr_ptr;
extern Array3D<double> dphi_kdr_c;

extern double *dphidz_ptr;
extern Array3D<double> dphidz_c;

extern double *dphi_kdz_ptr;
extern Array3D<double> dphi_kdz_c;

extern double *d2phidr2_ptr;
extern Array3D<double> d2phidr2_c;

extern double *d2phi_kdr2_ptr;
extern Array3D<double> d2phi_kdr2_c;

extern double *d2phidz2_ptr;
extern Array3D<double> d2phidz2_c;

extern double *d2phi_kdz2_ptr;
extern Array3D<double> d2phi_kdz2_c;

extern double *oppphi_ptr; 
extern Array3D<double> OPPphi_c;

extern double *oppphik_ptr;
extern Array3D<double> OPPphik_c;

extern double *l_hand_ptr;
extern Array3D<double> l_hand_c;

extern double *r_hand_ptr;
extern Array3D<double> r_hand_c;

extern double *upar_ptr;
extern Array3D<double> upar_c;

extern double *apars_ptr;
extern Array3D<double> apars_c;

extern double *apar_ptr;
extern Array3D<double> apar_c;

extern double *jpar_ptr;
extern Array3D<double> jpar_c;

extern double *rho_ptr;
extern Array3D<double> rho_c;

extern double *dene_ptr;
extern Array3D<double> dene_c;

extern double *delbx_ptr;
extern Array3D<double> delbx_c;

extern double *delby_ptr;
extern Array3D<double> delby_c;

extern double *delbz_ptr;
extern Array3D<double> delbz_c;

extern double *phis_ptr;
extern Array3D<double> phis_c;

extern double *denes_ptr;
extern Array3D<double> denes_c;

extern double *upars_ptr;
extern Array3D<double> upars_c;

//4D arrays
extern double *den_ptr;
extern Array4D<double> den_c;

void Allocate2dPointerArrays_gemx_com();
void Allocate3dPointerArrays_gemx_com();
void Allocate4dPointerArrays_gemx_com();