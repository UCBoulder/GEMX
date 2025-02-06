#include "gemx_c.hpp"
#include "gemx_com_externs.h"
#include "equil_externs.h"
//#include "pputil_c.hpp" //WIP
#include "fcnt.hpp"
#include "mpi.h"

#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <iomanip>


using namespace std;

//const auto start = std::chrono::high_resolution_clock::now(); //timing stuff
//timing stuff
   // const auto end = std::chrono::high_resolution_clock::now();
   // const std::chrono::duration<double> diff = end - start;
   // printf("%.10f\n",diff);


void parperp_c_(double& vpar,double& vperp2, const int& m, const int& cnt, const int& MyId){ //Working
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
   vpar = temp*iflag;

   vperp2 = -2.0*log(r2); 
   return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double ran2_c_(int& idum){ //Working
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
   
   const double AM = 1.0/IM1; 
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
      if((-1 * idum) < 1){ 
         idum = 1;
      }
      else{ 
         idum = (-1 * idum);
      }
      idum2 = idum;
      for(j = NTAB+7; j >= 0; j--){
         k = idum/IQ1;
         idum = IA1*(idum-k*IQ1)-k*IR1;

         if(idum < 0) idum = idum+IM1;
         if(j < NTAB) iv[j] = idum;
      }
      iy = iv[0];
   }

   k = idum/IQ1;
   idum = IA1*(idum-k*IQ1)-k*IR1;
   if(idum < 0) idum = idum+IM1;
   k = idum2/IQ2;
   idum2 = IA2*(idum2-k*IQ2)-k*IR2;
   if(idum2 < 0) idum2 = idum2+IM2;  
   j = iy/NDIV;
   iy = iv[j]-idum2; 
   iv[j] = idum;
   if(iy<1) iy = iy+IMM1;
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
void loadi_c_(){ //Working
   int MyId = 0; //TODO - Ok for now - when running on multiple MPI processes Mpi must be global
   int i = 0; 
   int j = 0;
   int k = 0;
   int m = 0;

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
   while(m < mm_ptr[0]){
   //load a slab of ions...

      //dumx=xdim*(ran2(iseed)+0.01)*0.9
      //dumy=zdim*(ran2(iseed)+0.01)*0.9
      //revers(MyId*cnt+j,2) !ran2(iseed)
      dumx=2*dxeq+(xdim-4*dxeq)*ran2_c_(iseed); 
      dumy=2*dzeq+(zdim-4*dzeq)*ran2_c_(iseed);
      dumz=pi2*ran2_c_(iseed);   

      //dumx=dxeq+(xdim-2*dxeq)*m/((mm(1)))

      r = xctr-xdim/2+dumx;   
      jacp = r/(xctr+xdim/2);
//    if(ran2(iseed)<jacp){
//       x2(m)=min(dumx,xdim-dxeq)
//       z2(m)=min(dumy,zdim-dzeq)
//       x2(m)=max(dumx,dxeq)
//       z2(m)=max(dumz,dzeq)
//       }
      zeta2_ptr[m] = dumz;
      x2_ptr[m] = dumx;
      z2_ptr[m] = dumy;

      parperp_c_(vpar, vperp2, m+1, cnt, MyId);

      x = x2_ptr[m];
      i = static_cast<int>(x/dxeq);
      wx0 = ((i+1)*dxeq-x)/dxeq;
      wx1 = 1.0 - wx0;

      z = z2_ptr[m];
      k = static_cast<int>(z/dzeq);
      wz0 = ((k+1)*dzeq-z)/dzeq;
      wz1 = 1-wz0;

      bfldp = wx0*wz0*b0_c(i,k)+wx0*wz1*b0_c(i,k+1)+wx1*wz0*b0_c(i+1,k)+wx1*wz1*b0_c(i+1,k+1);
      ter = wx0*wz0*t0i_c(i,k)+wx0*wz1*t0i_c(i,k+1)+wx1*wz0*t0i_c(i+1,k)+wx1*wz1*t0i_c(i+1,k+1);
      u2_ptr[m] = vpar/sqrt(mims_ptr[0]/ter);
      mu_ptr[m] = 0.5*vperp2/bfldp*ter;

      myavgv = myavgv+u2_ptr[m];
//    LINEAR: perturb w(m) to get linear growth...
//       w2(m)=2.*amp*ran2(iseed)
      w2_ptr[m] = (wx0*wz0*xn0i_c(i,k)+wx0*wz1*xn0i_c(i,k+1) 
               +wx1*wz0*xn0i_c(i+1,k)+wx1*wz1*xn0i_c(i+1,k+1))*r/xctr*((imx-3)*(jmx-3)*(kmx+1))/(numprocs*mmx);
//    w2(m) = r/xctr*((imx-1)*(jmx-1)*(kmx+1))/(numprocs*mmx)

      myavgw += w2_ptr[m];
      m++;
   }
   //do i=1,mmx
   //avex=avex+x2(i)
   //end do
   //write(*,*)avex/mmx
   if(MyId==0){ 
      ofstream myFile;
      string fileName = "testdepo_posi"; 
      myFile.open(fileName, ios::app);

      j = mmx-10001;
      while(j < mmx){
         myFile << setprecision(16) << x2_ptr[j] << "      ";
         myFile << setprecision(16) << z2_ptr[j] << "      ";
         myFile << setprecision(16) << zeta2_ptr[j] << "      " << endl;
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

   u_.CreateArray3D(u, imx, jmx, kmx); 
   ux_.CreateArray3D(ux, imx, jmx, kmx);
   uz_.CreateArray3D(uz, imx, jmx, kmx);

   int ju = 0;
   int jl = 0;
   double ul = 0;
   //const auto start = std::chrono::high_resolution_clock::now(); //timing stuff
   for(int j = 0; j < jmx; ++j){
      ju = j+1;
      jl = j-1;
      if(j == 0) jl = jmx-1;
      for (int i = 0; i < imx; ++i){
         for (int k = 0; k <= kmx; ++k){
            uz_(i,j,k) = (u_(i,ju,k)-u_(i,jl,k))/(2.*dz);
         }
      }  
   }

   for(int i = 1; i < imx; ++i){
      for(int j = 0; j < jmx; ++j){
         for(int k = 0; k <= kmx; ++k){
            ux_(i,j,k) = (u_(i+1,ju,k)-u_(i-1,j,k))/(2.*dx);
         }
      }
   }

   for(int j = 0; j < jmx; ++j){
      for(int k = 0; k <= kmx; ++k){
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
   uz_.CreateArray3D(uz, imx, jmx, kmx);
   //Rgrid is a 1D pointer declared in equil. I tried externing it to C but it wasn't working. For now, I'm passing a pointer to the first element via arguments

   for(int k = 0; k <= kmx; ++k)
   {
      kleft = k-1;
      if(k==0) kleft = kmx;
      kright = k+1;
      if(k==kmx) kright = 0;
      for(int i = 1; i < imx; ++i){
         for(int j = 1; j < jmx; ++j)
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

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALDER Flux Average SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void fluxavg_c_(int &i3D, double *phi_in, double *phiavg_in){ //Working

    //TODO READ FROM GEMXIN FILE - i3D SET THERE

   //Input is 3D array to be flux averaged 
   //input is phi in all cases, passed as externed 3D array as of now
   //Output is 2D interpolated array
   //Output is phiavg in all cases, passed as extrned 2D array as of now

   //Local Variables
   double phiavg1d[102];//index 0 - 101
   double psi1d[102];
   double phiavg1d_private[102];
   int gi = 0, xix, yjy, miw, psi_zero, store, k; 
   double weightinput,weightinput3D, phiavggi, psival, wmx0, wmx1;

   Array3D<double> phi_c;
   phi_c.CreateArray3D(phi_in, imx, jmx, kmx);
   Array2D<double> phiavg_c;
   phiavg_c.CreateArray2D(phiavg_ptr, nx, nz);

   //set all arrays to zero here.
   std::fill(std::begin(phiavg1d), std::end(phiavg1d), 0);
   std::fill(std::begin(psi1d), std::end(psi1d), 0);
   std::fill(std::begin(phiavg1d_private), std::end(phiavg1d_private), 0);

   psi_zero = 1;

   //Computation
   for(int line = 1; line <= 101; ++line){
      phiavggi = 0;
      store = 0;
      for(gi = 0; gi < num_lines; ++gi){
         if(gindex_ptr[gi] == line-1){
            weightinput = (weight00_ptr[gi]*phi_c(iarray_ptr[gi], jarray_ptr[gi], 0) + 
                           weight10_ptr[gi]*phi_c(iarray_ptr[gi]+1, jarray_ptr[gi], 0) + 
                           weight01_ptr[gi]*phi_c(iarray_ptr[gi], jarray_ptr[gi]+1, 0) + 
                           weight11_ptr[gi]*phi_c(iarray_ptr[gi]+1, jarray_ptr[gi]+1, 0));

            if(i3D == 0){
               phiavggi = phiavggi + (weightinput*jacobian_ptr[gi])/deno_ptr[gi]; 
            }else{
               for(k = 1; k <= kmx; ++k){
                  weightinput3D = (weight00_ptr[gi]*phi_c(iarray_ptr[gi], jarray_ptr[gi],k) + 
                           weight10_ptr[gi]*phi_c(iarray_ptr[gi]+1, jarray_ptr[gi],k) + 
                           weight01_ptr[gi]*phi_c(iarray_ptr[gi], jarray_ptr[gi]+1,k) + 
                           weight11_ptr[gi]*phi_c(iarray_ptr[gi]+1, jarray_ptr[gi]+1,k));
               }
               weightinput = weightinput3D + weightinput;
               phiavggi = phiavggi +(weightinput*jacobian_ptr[gi])/(deno_ptr[gi]*(kmx+1)); 
            }

            if(priv_ptr[gi] == 0){
               store = gi;
            }
         }
         //Remove redundancy from closed loop integration process
         if(phiavggi != 0){
            if(gindex_ptr[gi] != line-1){
               if(i3D == 0){
                  phiavggi = phiavggi - (weightinput*jacobian_ptr[gi-1])/deno_ptr[gi-1]; 
               }else{
                  phiavggi = phiavggi - (weightinput*jacobian_ptr[gi-1])/(deno_ptr[gi-1]*(kmx+1));               
               }
               break;
            }
         }
      }
      if(store != 0){
            phiavg1d[line] = phiavggi;
            psi1d[psi_zero] = psitab_ptr[store];
            psi_zero += 1;
         }else{
            phiavg1d_private[line] = phiavggi; 
         }
   }
   phiavg1d[0] = phiavg1d[1];
   //  phiavg1d(0) = input(268,254,0) //phiavg1d[0] = phi_c(268,254,0);
   
   // Save psi1d and timesteps of phiavg1d to understand convergence
   //  if (timestep == 10) then
   //     open(unit=11, file = 'psi1d',status='unknown',action='write')
   //                 write(11,*) psi1d(:)
   //              close(11)
   //  endif

   //  open(unit=11, file = 'phiavg1d',status='unknown',position='append')                
   //  write(11,*) phiavg1d(:)
   //  close(11)

   // Initialize output to zero

   //INTERPOLATION
   for(xix = 0; xix <= nx; ++xix){
      for(yjy = 0; yjy <= nz; ++yjy){
         psival = psi_p_c(xix, yjy);
         if(mask_c(xix,yjy) < 0.99){ 
            phiavg_c(xix,yjy) = 0;
         }else{
            miw  = int(psival/(psi1d[2]-psi1d[1]));
            wmx0 = ((miw+1)*(psi1d[2]-psi1d[1])-psival)/(psi1d[2]-psi1d[1]);
            wmx1 = 1-wmx0;
            if (yjy < 75 && xix < 150 && psival > 0.29 && psival < 0.31){ //Private region under X-point
               phiavg_c(xix,yjy) = wmx0*phiavg1d_private[miw] + wmx1*phiavg1d_private[miw+1];
            }else{
               phiavg_c(xix,yjy) = wmx0*phiavg1d[miw] + wmx1*phiavg1d[miw+1]; 
            }
         }
      }
   }
}