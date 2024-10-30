#include "gemx.hpp"
#include "gemx_com_c.h"
#include "fcnt.hpp"
#include "mpi.h"

#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>


using namespace std;

void parperp_c_(double* vpar,double* vperp2,int& m,int& cnt,int& MyId){
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
   r1 = revers_c(m+MyId*cnt, 7); //sets values for r1 and r2 (random numbers)
   r2 = revers_c(m+MyId*cnt, 11); 

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

void ran2_c_(int& idum, double *r2Val){ //double *r2Val
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

   static int idum2 = 123456789;
   static int iy = 0;
   static int iv[NTAB];

   // static int idum2 = 300128977;
   // static int iy = 935447419; //(Stick with me here) Produces results that agree with fortran, but there may be an un-intended bug here: iy is set to equal 0 in fortran, but it's not actually 0
   // static int iv[NTAB] = {937875999, 1663474829, 1028772421, 1588498714, 376286980, 584250848, 865121955, 739957527, 455331591, 1717783064, 1358293836, 514315536, 1721350347, 1601205264, 1593124450, 1938785586, 707365054, 499070822, 849022444, 969044322, 1410893119, 387733782, 758231585, 1818821612, 987322158, 110032498, 2001988383, 1936233189, 900462814, 2131290336, 187234747, 754377697};
   
   double temp;
   double retVal = 0.0;

   if(idum <= 0){
      if(-idum < 1) idum = 1;
      else idum = -idum;

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

   *r2Val = retVal;
   //return retVal;

}

// void loadi_c_(){ 
//    int MyId = 0; //TODO - Ok for now - when running on multiple MPI processes Mpi must be global
//    int i = 0; 
//    int j = 0;
//    int k = 0;
//    int m = 1;
//    int idum = 0;
//    int ns  = 0;
//    int m1 = 0;

//    double vpar = 0.0;
//    double vperp2 = 0.0;
//    double r = 0;
//    double x = 0;
//    double z = 0;
//    double b = 0;
//    double ter = 0;
//    double bfldp = 0;

//    double avgv = 0;
//    double myavgv = 0;
//    double avgwl = 0;
//    double myavgw = 0;

//    double dumx, dumy, dumz, jacp, rand[4];
//    double wx0, wx1, wz0, wz1, avex = 0;
//    double r2Val1, r2Val2, r2Val3;

//    double pi2 = M_PI*2;

//    cnt = int(tmm_ptr[1]/numprocs);
//    cnt = mmx; 

//    //double ran2Val = ran2_c_(iseed);

//    while(m<=mm_ptr[1]){
//       ran2_c_(iseed, &r2Val1);
//       ran2_c_(iseed, &r2Val2); //need to do this in order to have ran2_c called in fortran and c++ properly for now
//       ran2_c_(iseed, &r2Val3);

//       dumx = 2*dxeq+(xdim-4*dxeq)*r2Val1; 
//       dumy=2*dzeq+(zdim-4*dzeq)*r2Val2;
//       dumz=pi2*r2Val3;   

//       r = xctr-xdim/2+dumx;   
//       jacp = r/(xctr+xdim/2);

//       zeta2_ptr[m] = dumz;
//       x2_ptr[m] = dumx;
//       z2_ptr[m] = dumy;

//       parperp_c_(&vpar, &vperp2, m, cnt, MyId);

//       x = x2_ptr[m];
//       i = int(x/dxeq);
//       wx0 = ((i+1)*dxeq-x)/dxeq;
//       wx1 = 1.0 - wx0;

//       z = z2_ptr[m];
//       k = int(z/dzeq);
//       wz0 = ((k+1)*dzeq-z)/dzeq;
//       wz1 = 1-wz0;

//       bfldp = wx0*wz0*b0_c(i,k)+wx0*wz1*b0_c(i,k+1) 
//          +wx1*wz0*b0_c(i+1,k)+wx1*wz1*b0_c(i+1,k+1);
//       ter = wx0*wz0*t0i_c(i,k)+wx0*wz1*t0i_c(i,k+1) 
//          +wx1*wz0*t0i_c(i+1,k)+wx1*wz1*t0i_c(i+1,k+1);

//       u2_ptr[m] = vpar/sqrt(mims_ptr[1]/ter);
//       mu_ptr[m] = 0.5*vperp2/bfldp*ter;

//       myavgv = myavgv+u2_ptr[m];
// //    LINEAR: perturb w(m) to get linear growth...
// //       w2(m)=2.*amp*ran2(iseed)
//       w2_ptr[m] = (wx0*wz0*xn0i_c(i,k)+wx0*wz1*xn0i_c(i,k+1) 
//                +wx1*wz0*xn0i_c(i+1,k)+wx1*wz1*xn0i_c(i+1,k+1))*r/xctr*((imx-3)*(jmx-3)*(kmx+1))/(numprocs*mmx);
// //    w2(m) = r/xctr*((imx-1)*(jmx-1)*(kmx+1))/(numprocs*mmx)

//       myavgw = myavgw+w2_ptr[m];
//       m = m+1;
//    }
   
      

//    if(MyId==0){ //opens file and writes position to file
//       ofstream myFile;
//       string fileName = "testdepo_posi"; //not sure why I needed to do this - .open('blah.txt') was reading an integer for 'blah.txt'. 
//       myFile.open(fileName, ios::app);

//       int j = mmx-100000;
//       do{
//          myFile << x2_ptr[j];
//          myFile << z2_ptr[j];
//          myFile << zeta2_ptr[j] << endl;
//          j++;
//       }while(j <= mmx);
//       myFile.close();
//    }

//    myavgw = myavgw/mm_ptr[1];

//    ierr = MPI_Allreduce(&myavgv, &avgv, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD);

//    if(idg == 1) cout << "all reduce" << endl;
//    avgv = avgv/float(tmm_ptr[1]);

//    m = 1;
//    do{
//       u2_ptr[m] = u2_ptr[m]-avgv;
//       x3_ptr[m] = x2_ptr[m];
//       z3_ptr[m] = z2_ptr[m];
//       zeta3_ptr[m] = zeta2_ptr[m];
//       u3_ptr[m] = u2_ptr[m];
// //    w2(m) = w2(m)-myavgw
//       w3_ptr[m] = w2_ptr[m];
//       m++;
//    }while(m <= mm_ptr[1]);

//    return;
// }