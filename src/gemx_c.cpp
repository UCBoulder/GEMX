#include "gemx.hpp"
#include "gemx_com_c.h"
#include "fcnt.hpp"
#include "mpi.h"

#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>


using namespace std;

void parperp_c_(double* vpar,double* vperp2, const int& m, const int& cnt, const int& MyId){
   double r1 = 0.0;
   double r2 = 0.0;
   double t = 0.0;

   const double c0 = 2.515517;
   const double c1 = 0.802853;
   const double c2 = 0.010328;

   const double d1 = 1.432788;
   const double d2 = 0.189269;
   const double d3 = 0.001308;
   double temp = 0.0;
   int iflag = 1;

   //const auto start = std::chrono::high_resolution_clock::now(); //timing stuff
   r1 = revers_c_(m+MyId*cnt, 7); //sets values for r1 and r2 (random numbers)
   r2 = revers_c_(m+MyId*cnt, 11); 

   //.....quiet start---see denavit pf '71(?) & abramowitz hand book
   //.....fibonacci start---see denavit comm. pla. phy. & con. fus. '81
   // warning: we have g1=1 in the x-direction. This surpresses all odd
   //          modes in the x-direction!!!

   if (r1 <= 0.5) goto jump;

   r1 = 1.0 - r1;
   iflag = -1;

   jump:
   if(r1 >= 1.0e-6){
      t = sqrt(log(1.0/(r1*r1)));
   }
   else{
      t = 5.0;
      cout << "parperp2 warning m= " << m << endl;
   }

   temp = t-(c0+c1*t+c2*(t*t))/(1.+d1*t+d2*(t*t)+d3*(t*t*t));
   *vpar = temp*iflag;

   *vperp2 = -2.0*log(r2); 

   //timing stuff
   // const auto end = std::chrono::high_resolution_clock::now();
   // const std::chrono::duration<double> diff = end - start;
   // printf("%.10f\n",diff);

   return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double ran2_c_(int& idum){ 
   const int IM1 = 2147483563;
   const int IM2 = 2147483399;
   const int IMM1 = IM1-1;
   const int IA1 = 40014;
   const int IA2 = 40692;
   const int IQ1 = 53668;
   const int IQ2 = 52774;
   const int IR1 = 12211;
   const int IR2 = 3791;
   const int NTAB = 32; 
   
   const double AM = 1.0/IM1; //4.6566131e-10 
   const double NDIV = 1+IMM1/NTAB;
   const double EPS = 1.2e-7;
   const double RNMX = 1.0-EPS;

   int j = 0;
   int k = 0;

   static long idum2 = 123456789;
   static long iy = 0;
   static long iv[NTAB];

   double temp;
   double retVal;

   if(idum <= 0){
      if((-1 * idum) < 1)
      { 
         idum = 1;
      }
      else
      { 
         idum = (-1 * idum);
      }
      idum2 = idum;
      for(j = NTAB+7; j >= 0; j--){
         k = idum/IQ1;
         idum = IA1*(idum-k*IQ1)-k*IR1;

         if(idum < 0)
         {
            idum = idum+IM1;
         } 
         if(j < NTAB)
         {
            iv[j] = idum;
         } 
      }
      iy = iv[0];
   }

   k = idum/IQ1;
   idum = IA1*(idum-k*IQ1)-k*IR1;
   if(idum < 0)
   {
      idum = idum+IM1;
   } 
   k = idum2/IQ2;
   idum2 = IA2*(idum2-k*IQ2)-k*IR2;
   if(idum2 < 0)
   {
     idum2 = idum2+IM2; 
   } 
   j = iy/NDIV;
   iy = iv[j]-idum2; 
   iv[j] = idum;
   if(iy<1)
   {
      iy = iy+IMM1;
   } 
   temp = AM*iy;

   if(temp>RNMX){
      retVal = RNMX;
   }else{
      retVal = temp;
   }
   return retVal;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void loadi_c_(){ 
   int MyId = 0; //TODO - Ok for now - when running on multiple MPI processes Mpi must be global
   int i = 0; 
   int j = 0;
   int k = 0;
   int m = 1;
   //double debug1 = 0;
   //double debug2 = 0;

   double vpar = 0.0;
   double vperp2 = 0.0;
   double r = 0;
   double x = 0;
   double z = 0;
   // double b = 0;
   double ter = 0;
   double bfldp = 0;

   double avgv = 0;
   double myavgv = 0;
   // double avgw = 0;
   double myavgw = 0;
   
   double dumx, dumy, dumz, jacp;
   double wx0, wx1, wz0, wz1;

   const long double pi2 = M_PI*2;

   cnt = static_cast<int>(tmm_ptr[0]/numprocs);
   cnt = mmx; 

 //needed I think in order to properly index through the arrays. Most 1D arrays in fortran start at 1 instead of 0, offsetting everything by 1
   while(m <= mm_ptr[0]){
      dumx=2*dxeq+(xdim-4*dxeq)*ran2_c_(iseed); 
      dumy=2*dzeq+(zdim-4*dzeq)*ran2_c_(iseed);
      dumz=pi2*ran2_c_(iseed);   

      r = xctr-xdim/2+dumx;   
      jacp = r/(xctr+xdim/2);

      zeta2_ptr[m-1] = dumz;
      x2_ptr[m-1] = dumx;
      z2_ptr[m-1] = dumy;

      parperp_c_(&vpar, &vperp2, m, cnt, MyId);

      x = x2_ptr[m-1];
      i = static_cast<int>(x/dxeq);
      wx0 = ((i+1)*dxeq-x)/dxeq;
      wx1 = 1.0 - wx0;

      z = z2_ptr[m-1];
      k = static_cast<int>(z/dzeq);
      wz0 = ((k+1)*dzeq-z)/dzeq;
      wz1 = 1-wz0;

      bfldp = wx0*wz0*b0_c(i,k)+wx0*wz1*b0_c(i,k+1)+wx1*wz0*b0_c(i+1,k)+wx1*wz1*b0_c(i+1,k+1);
      //debug1 = b0_c(i,k);
      ter = wx0*wz0*t0i_c(i,k)+wx0*wz1*t0i_c(i,k+1)+wx1*wz0*t0i_c(i+1,k)+wx1*wz1*t0i_c(i+1,k+1);
      //debug2 = t0i_c(i,k);
      u2_ptr[m-1] = vpar/sqrt(mims_ptr[0]/ter);
      mu_ptr[m-1] = 0.5*vperp2/bfldp*ter;

      myavgv = myavgv+u2_ptr[m-1];
//    LINEAR: perturb w(m) to get linear growth...
//       w2(m)=2.*amp*ran2(iseed)
      w2_ptr[m-1] = (wx0*wz0*xn0i_c(i,k)+wx0*wz1*xn0i_c(i,k+1) 
               +wx1*wz0*xn0i_c(i+1,k)+wx1*wz1*xn0i_c(i+1,k+1))*r/xctr*((imx-3)*(jmx-3)*(kmx+1))/(numprocs*mmx);
//    w2(m) = r/xctr*((imx-1)*(jmx-1)*(kmx+1))/(numprocs*mmx)

      myavgw += w2_ptr[m-1];
      m++;
   }

   if(MyId==0){ //opens file and writes position to file
      ofstream myFile;
      string fileName = "testdepo_posi"; //not sure why I needed to do this - .open('blah.txt') was reading an integer for 'blah.txt'. 
      myFile.open(fileName, ios::app);

      j = mmx-10001;
      while(j <= mmx){
         myFile << x2_ptr[j] << "      ";
         myFile << z2_ptr[j] << "      ";
         myFile << zeta2_ptr[j] << "      " << endl;
         j++;
      }
      myFile.close();
   }

   myavgw = myavgw/mm_ptr[0];

   ierr = MPI_Allreduce(&myavgv, &avgv, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD);

   if(idg == 1) cout << "all reduce" << endl;
   avgv = avgv/static_cast<float>(tmm_ptr[0]);

   m = 0;
   do{
      u2_ptr[m] = u2_ptr[m]-avgv;
      x3_ptr[m] = x2_ptr[m];
      z3_ptr[m] = z2_ptr[m];
      zeta3_ptr[m] = zeta2_ptr[m];
      u3_ptr[m] = u2_ptr[m];
//    w2(m) = w2(m)-myavgw
      w3_ptr[m] = w2_ptr[m];
      m++;
   }while(m <= mm_ptr[0]);

   return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void gradu_c_(double* u, double* ux, double* uz){
   Array3D<double> u_;
   Array3D<double> ux_;
   Array3D<double> uz_;

   u_.CreateArray3D(u, imx,jmx,1);
   ux_.CreateArray3D(ux, imx,jmx,kmx);
   uz_.CreateArray3D(uz, imx, jmx, kmx);

   int ju = 0;
   int jl = 0;
   double ul = 0;

   for(int j = 0; j < imx; j++){
      ju = j+1;
      jl = j-1;
      if(j == 0) jl = jmx-1;
      for (int i = 0; i < imx; i++){
         for (int k = 0; k <= kmx; i++){
            uz_(i,j,k) = (u_(i,ju,k)-u_(i,jl,k))/(2.*dz);
         }
      }  
   }

   for(int i = 1; i < imx; i++){
      for(int j = 0; j < jmx; j++){
         for(int k = 0; k <= kmx; k++){
            ux_(i,j,k) = (u_(i+1,ju,k)-u_(i-1,j,k))/(2.*dx);
         }
      }
   }

   for(int j = 0; j < jmx; j++){
      for(int k = 0; k <= kmx; k++){
         ul = u_(imx-1, j, k);
         ux_(0,j,k) = (u_(1,j,k)-ul)/(2.*dx);
      }
   }
   return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void gradz_c_(double* u, double* uz, double* Rgrid_c){
   int kleft, kright;
   double wx0, wx1, wz0, wz1, uleft, uright;

   Array3D<double> u_;
   Array3D<double> uz_;

   u_.CreateArray3D(u, imx, jmx, kmx);
   uz_.CreateArray3D(u, imx, jmx, kmx);
   //Rgrid is a 1D pointer declared in equil. I tried externing it to C but it wasn't working. For now, I'm passing a pointer to the first element via arguments

   *uz = 0.;
   for(int k = 0; k <= kmx; k++)
   {
      kleft = k-1;
      if(k==0) kleft = kmx;
      kright = k + 1;
      if(k==kmx) kright = 0;
      for(int i = 1; i < imx; i++){
         for(int j = 1; j < jmx; j++)
         {
            wx0 = ((ileft_c(i,j)+1)*dx-xbackw_c(i,j))/dx;
               wx1 = 1.0-wx0;
               wz0 = ((jleft_c(i,j)+1)*dz-zbackw_c(i,j))/dz;
               wz1 = 1.0-wz0;
               uleft = wx0*wz0*u_(ileft_c(i,j),jleft_c(i,j),kleft) 
                      +wx1*wz0*u_(ileft_c(i,j)+1,jleft_c(i,j),kleft) 
                      +wx0*wz1*u_(ileft_c(i,j),jleft_c(i,j)+1,kleft) 
                      +wx1*wz1*u_(ileft_c(i,j)+1,jleft_c(i,j)+1,kleft);
               wx0 = ((iright_c(i,j)+1)*dx-xforw_c(i,j))/dx;
               wx1 = 1.0-wx0;
               wz0 = ((jright_c(i,j)+1)*dz-zforw_c(i,j))/dz;
               wz1 = 1.0-wz0;
               uright = wx0*wz0*u_(iright_c(i,j),jright_c(i,j),kright) 
                      +wx1*wz0*u_(iright_c(i,j)+1,jright_c(i,j),kright) 
                      +wx0*wz1*u_(iright_c(i,j),jright_c(i,j)+1,kright) 
                      +wx1*wz1*u_(iright_c(i,j)+1,jright_c(i,j)+1,kright);

               uz_(i,j,k)=(uright-uleft)/(2.*b0_c(i,j)/b0zeta_c(i,j)*dzeta*(Rgrid_c[i])/xu);
         }
      }
   }
   return;
}