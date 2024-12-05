#pragma once

extern "C"
{
   void parperp_c_(double* vpar,double* vperp2,int& m,int& cnt,int& MyId);
   void loadi_c_();
   double ran2_c_(int& idum);
   void gradu_c_(double* u, double* ux, double* uz);
   void gradz_c_(double* u, double* uz, double* Rgrid); 
}