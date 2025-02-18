#pragma once

extern "C"
{
   void parperp_c_(double& vpar, double& vperp2, const int& m, const int& cnt);
   void loadi_c_();
   double ran2_c_(int& idum);
   void gradu_c_(double* u, double* ux, double* uz);
   void gradz_c_(double* u, double* uz);
   void fluxavg_c_(double *phi_in, double *phiavg_in);
   void efieldcalc_c_(double *phi_input);
   void growthdiag_c_(double *input_phi); 
   void BoltzSolve_c_(double *input_phi);
   void smooth_c_(double* matrix, int &mk);
   void get_jpar_c_(double* matrix);
   void integ_c_(int &iflag);
   void grid1_c_(int &ip, int &n, int &MyId);
}