//actual gemx_com_c file to convert
//include this header to access all variables and arrays allocated in gemx_com_c
#include "MultiArraysC.hpp"
#include "mpi.h"
#include <complex>
#pragma once

void cleanupCom();
void new_gemx_com();
//Variables and Constants
extern int imx, jmx, kmx, mmx;
extern int nmx,nsmx,nsubd,ntube,petsc_color,petsc_rank,iBoltzmann,globle_integer=0,eBoltzmann,eAdiabatic,iterations,dbg;
extern int rand_table[]; //10007
extern int timestep,iez;
extern int iseed;
extern int nm,nsm,ncurr,iflr,ifield_solver,ntracer,i3D,icollision;
extern int iput,iget,ision,isham,peritr,iadi;
extern int idg;
extern int nonlin,nonline,iflut,ifluid,ipara;
extern int nplot, xnplt;
extern int numprocs;
extern int last,myid, cnt , ierr;
extern char outname[]; //71
extern double endtm,begtm,pstm;
extern double starttm, lasttm, tottm;
extern double start_total_tm, end_total_tm, start_integ_tm, end_integ_tm, start_ppush_tm, end_ppush_tm, start_cpush_tm, end_cpush_tm;
extern double total_tm, integ_tm, ppush_tm, cpush_tm;
extern double dx,dz,dzeta,pi,pi2,dt,totvol,n0,tcurr;
extern double etaohm;
extern double lx,lz;
extern double cut,amp,tor,amie,emass,qel,rneu;
extern double vcut;
extern const int master;
extern MPI_Comm TUBE_COMM, GRID_COMM, PETSC_COMM;
//1D arrays
extern int *mm;
extern int *tmm;
extern int *lr;

extern double *mims;
extern double *q;
//extern double *time;
extern double *xg;
extern double *zg;
extern double *jac;
extern double *fe;
extern double *te;
extern double *rmsphi;
extern double *rmsez;
extern double *rmsapa;
extern double *avewi;
extern double *vol;
extern double *mu;
extern double *x2;
extern double *zeta2;
extern double *z2;
extern double *u2;
extern double *x3;
extern double *zeta3;
extern double *z3;
extern double *u3;
extern double *w2;
extern double *w3;

//2D Arrays
extern CArray2D<double> den2d1;
extern CArray2D<double> den2d2;
extern CArray2D<double> dden2d;

extern CArray2D<double> bmag;
extern CArray2D<double> gbtor;
extern CArray2D<double> gbx;
extern CArray2D<double> gbz;

extern CArray2D<double> gn0i;
extern CArray2D<double> gn0e;
extern CArray2D<double> gt0i;
extern CArray2D<double> gt0e;
extern CArray2D<double> xforw;
extern CArray2D<double> zforw;
extern CArray2D<double> zbackw;
extern CArray2D<double> xbackw;

extern CArray2D<double> gcpnex;
extern CArray2D<double> gcpnez;
extern CArray2D<double> gcptex;
extern CArray2D<double> gcptez;

extern CArray2D<double> gnuobx;
extern CArray2D<double> gnuoby;
extern CArray2D<double> gupae0;

extern CArray2D<int> ileft;
extern CArray2D<int> jleft;
extern CArray2D<int> iright;
extern CArray2D<int> jright;

extern CArray2D<double> ke; 
extern CArray2D<double> nos; 

extern CArray2D<double> efle; 
extern CArray2D<double> pfle; 
extern CArray2D<double> pfl; 
extern CArray2D<double> efl;

//3D Arrays
extern CArray3D<double> rho;
extern CArray3D<double> phi;

extern CArray3D<double> ex;
extern CArray3D<double> ez; 
extern CArray3D<double> ezeta;

extern CArray3D<double> delbx; 
extern CArray3D<double> delbz; 
extern CArray3D<double> delby;

extern CArray3D<double> apar;
extern CArray3D<double> dene;

extern CArray3D<double> upar;
extern CArray3D<double> jpar;
extern CArray3D<double> upars;
extern CArray3D<double> phis;
extern CArray3D<double> denes;
extern CArray3D<double> apars;

extern CArray3D<double> dnedx;
extern CArray3D<double> dnedy;
extern CArray3D<double> dupadx;
extern CArray3D<double> dupady;
extern CArray3D<double> phi_k;
extern CArray3D<double> dphidr;
extern CArray3D<double> dphi_kdr;
extern CArray3D<double> d2phidr2;
extern CArray3D<double> d2phi_kdr2;
extern CArray3D<double> d2phidz2;
extern CArray3D<double> d2phi_kdz2;
extern CArray3D<double> dphidz;
extern CArray3D<double> dphi_kdz;
extern CArray3D<double> OPPphi;
extern CArray3D<double> OPPphik;
extern CArray3D<double> l_hand;
extern CArray3D<double> r_hand;

//4D arrays
extern CArray4D<double> den;