#include "MultiArrays.hpp"
#pragma once

extern "C"
{
   void parperp_c_(double& vpar, double& vperp2, const int& m, const int& cnt);
   void loadi_c_();
   double ran2_c_(int& idum);
   void gradu_c_(double* u, double* ux, double* uz);
   void fluxavg_c_(double *phi_in, double *phiavg_in);
   void efieldcalc_c_(double *phi_input);
   void growthdiag_c_(double *input_phi); 
   void BoltzSolve_c_(double *input_phi);
   void smooth_c_(double* matrix, int &mk);
   void get_jpar_c_(double* matrix);
   void get_apar_c_(const int &flagnumber);
   void integ_c_(int &iflag);
   void grid1_c_(int &ip, int &n, int &MyId);
   void pintef_c_();
   void gradparz_c_(double *matrix);
   void get_ne_c_(int &flagnumber);
}
void gradz_c_(Array3D<double> &u, Array3D<double> &uz);
void gradpar_c_(Array3D<double> &matrix, Array3D<double> &gradPar);