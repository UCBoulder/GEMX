#include "MultiArraysC.hpp"
#include <petscsys.h>
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscerror.h>
#pragma once

// extern "C"
// {
void parperp_c_(double& vpar, double& vperp2, const int& m, const int& cnt);
void loadi_c_();
double ran2_c_(int& idum);
void gradu_c_(CArray3D<double> u, CArray3D<double> ux, CArray3D<double> uz);
void fluxavg_c_(CArray3D<double> phi, CArray2D<double> phiavg_in);
void efieldcalc_c_(CArray3D<double> phi_input);
void growthdiag_c_(CArray3D<double> &input_phi); 
void BoltzSolve_c_(CArray3D<double> &input_phi);
void smooth_c_(CArray3D<double> matrix, int &mk);
void get_jpar_(CArray3D<double> matrix);
void get_apar_(const int &flagnumber);
void integ_c_(int iflag);
void grid1_c_(int &ip, int &n, int &MyId);
void pintef_c_();
void gradparz_c_(double *matrix);
void get_ne_c_(int flagnumber);
void initialize_c_();
// }
void gradz_c_(CArray3D<double> &u, CArray3D<double> &uz);
void gradpar_c_(CArray3D<double> &matrix, CArray3D<double> &gradPar);
void init();
PetscErrorCode ComputeInitialGuess(KSP ksp, Vec init_guess, void* ctx_void);
PetscErrorCode ComputeMatrix(KSP ksp, Mat AA, Mat BB, void* dummy);
PetscErrorCode ComputeRHS(KSP ksp, Vec bbb, void* k);