//Actual equil.h file
#include "MultiArraysC.hpp"
#include "gemx_com_c.hpp" 
#include "readDatFiles.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#pragma once
void new_equil_c();
void cleanUpEquil();

//variables extern to create one copy used everywhere

extern double mimp, chgi;
extern double betaVal,rmaj0,a,q0,r0,q0p,q0abs,shat0;
extern double R,dth,mu0,e,proton;
extern int nr,nr2,ntheta,isgnf,isgnq,isupae0,tor_n;
extern double psi_max, psi_min ,R_min, Z_min, Z_internal, psi_div,psi_a;
extern double phi_diag, phi_diag_freq;

extern int nzeta;
extern int nx, nz;
extern double zctr;
extern double  dxeq, xdim, xctr, zdim, dzeq;
extern double pi,pi2;

extern CArray2D<double> b0,b0x,b0z,b0zeta,dbdx,dbdz,c2_over_vA2,q_grid;
extern CArray2D<double> t0i,t0e,xn0i,xn0e,captix,captex,capnix,capnex,captiz,captez,capniz,capnez;

extern double *psi;
extern double *psip;
extern double *sf;
extern double *vpari;
extern double *vparip;
extern double *zeff;
extern double *nue0;
extern double *phinc;
extern double *phincp;
extern double *er;
extern double *upari;
extern double *Rgrid;
extern double *Zgrid;

extern CArray2D<double> t0s;
extern CArray2D<double> xn0s;
extern CArray2D<double> capts;
extern CArray2D<double> capns;
extern CArray2D<double> vpars ;
extern CArray2D<double> vparsp;
extern CArray2D<double> psi_p;
extern CArray2D<double> mask;
extern CArray2D<double> mask2;
extern CArray2D<double> mask3;
extern CArray2D<double> mask4;
extern double bu,tu,nu,xu,frequ,vu,eru;

extern CArray2D<double> bdcrvb;
extern CArray2D<double> upae0,nuob,dnuobdr,dnuobdt;
extern CArray3D<double> curlb;
extern CArray2D<double> rho_i;
extern CArray2D<double> dpsi_dr;
extern CArray2D<double> dpsi_dz;

extern int num_lines, line;
extern CArray2D<double> phiavg;
extern double *psitab;
extern double *weight00;
extern double *weight01;
extern double *weight10;
extern double *weight11;
extern double *jacobian;
extern double *deno;
extern int *gindex;
extern int *iarray;
extern int *jarray;
extern int *priv;