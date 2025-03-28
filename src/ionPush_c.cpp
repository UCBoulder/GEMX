#include "ionPush_c.h" 
#include "gemx_com_externs.h"
#include "equil_externs.h"
#include "mpi.h"

#include <cmath>
#include <iostream>

using namespace std;

//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//       Ion pre-push
//
void ppush_c_(int &n){
    double exp1,ezp,ezetap,delbxp,delbzp,energy, energy0,nudi0,nudi,T_center,ni_temp;
    double wx0,wx1,wy0,wy1,wz0,wz1,dum1;
    int m,i,j,k,l,k_plus_1;
    double rhog,vfac,kapxp,kapzp,vpar,kaptxp,kapnxp,kaptzp,kapnzp,xnp;
    double b,enerb,ter,z,zeta,bstar;
    double x;
    double xt,zt,xdot,zdot,zetadot,pzdot,edot;
    double dbdxp,dbdzp,bfldp,bfldxp,bfldzp,bfldzetap,dbdzetap=0;
    double rhox[4], rhoy[4], curlbp[3], Bstar3[3]; //Realize these are 0 indexed: 0, 1, 2, 3 not 1, 2, 3, 4
    //real(8),dimension(3)::curlbp,Bstar3
    start_ppush_tm = MPI_Wtime();

    nudi0 = 1/sqrt(2)*18.4*pow(e,1.5)*(4.7140*pow(10,-8))*(1*pow(10,-6));
    //write(*,*)t0i(200,201);
    T_center = t0i_c(imx/2, jmx/2);

    #pragma acc parallel loop gang vector private(rhoy,bstar3,rhox) copy(rand_table_ptr)
    for(m = 0; m < mm_ptr[0]; ++m){
        x = x2_ptr[m];
        i = static_cast<int>(x/dxeq);
        i = min(i,nx-1);
        wx0 = (i+1)-x/dxeq;
        wx1 = 1-wx0;

        z = z2_ptr[m];
        k = static_cast<int>(z/dzeq);
        k = min(k,nz-1);
        wz0 = (k+1)-z/dzeq;
        wz1 = 1-wz0;




        //bdcurlbp =wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
        //                 +wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1)
        for(j = 0; j < 3; ++j){
            curlbp[j]= wx0*wz0*curlb_c(i,k,j)+wx0*wz1*curlb_c(i,k+1,j) 
              +wx1*wz0*curlb_c(i+1,k,j)+wx1*wz1*curlb_c(i+1,k+1,j);
        }
        //write(*,*) curlbp(1)-wx0*wz0*curlb(i,k,1)-wx0*wz1*curlb(i,k+1,1)-wx1*wz0*curlb(i+1,k,1)-wx1*wz1*curlb(i+1,k+1,1),curlbp(2)-wx0*wz0*curlb(i,k,2)-wx0*wz1*curlb(i,k+1,2)-wx1*wz0*curlb(i+1,k,2)-wx1*wz1*curlb(i+1,k+1,2),curlbp(3)-wx0*wz0*curlb(i,k,3)-wx0*wz1*curlb(i,k+1,3)-wx1*wz0*curlb(i+1,k,3)-wx1*wz1*curlb(i+1,k+1,3)
        //         write(*,*)curlbp(1),curlbp(2),curlbp(3)
        // write(*,*)bdcurlbp
         dbdxp = wx0*wz0*dbdx_c(i,k)+wx0*wz1*dbdx_c(i,k+1) 
                 +wx1*wz0*dbdx_c(i+1,k)+wx1*wz1*dbdx_c(i+1,k+1); 
         dbdzp = wx0*wz0*dbdz_c(i,k)+wx0*wz1*dbdz_c(i,k+1) 
                 +wx1*wz0*dbdz_c(i+1,k)+wx1*wz1*dbdz_c(i+1,k+1);
         bfldp = wx0*wz0*b0_c(i,k)+wx0*wz1*b0_c(i,k+1) 
                 +wx1*wz0*b0_c(i+1,k)+wx1*wz1*b0_c(i+1,k+1); 
         bfldxp = wx0*wz0*b0x_c(i,k)+wx0*wz1*b0x_c(i,k+1) 
                 +wx1*wz0*b0x_c(i+1,k)+wx1*wz1*b0x_c(i+1,k+1); 
         bfldzp = wx0*wz0*b0z_c(i,k)+wx0*wz1*b0z_c(i,k+1) 
                 +wx1*wz0*b0z_c(i+1,k)+wx1*wz1*b0z_c(i+1,k+1); 
         bfldzetap = wx0*wz0*b0zeta_c(i,k)+wx0*wz1*b0zeta_c(i,k+1) 
                 +wx1*wz0*b0zeta_c(i+1,k)+wx1*wz1*b0zeta_c(i+1,k+1); 
         ter = wx0*wz0*t0i_c(i,k)+wx0*wz1*t0i_c(i,k+1) 
                 +wx1*wz0*t0i_c(i+1,k)+wx1*wz1*t0i_c(i+1,k+1); 
         kaptxp = wx0*wz0*captix_c(i,k)+wx0*wz1*captix_c(i,k+1) 
                 +wx1*wz0*captix_c(i+1,k)+wx1*wz1*captix_c(i+1,k+1);

         kapnxp = wx0*wz0*capnix_c(i,k)+wx0*wz1*capnix_c(i,k+1) 
                 +wx1*wz0*capnix_c(i+1,k)+wx1*wz1*capnix_c(i+1,k+1); 

         kaptzp = wx0*wz0*captiz_c(i,k)+wx0*wz1*captiz_c(i,k+1) 
                 +wx1*wz0*captiz_c(i+1,k)+wx1*wz1*captiz_c(i+1,k+1); 
         kapnzp = wx0*wz0*capniz_c(i,k)+wx0*wz1*capniz_c(i,k+1) 
                 +wx1*wz0*capniz_c(i+1,k)+wx1*wz1*capniz_c(i+1,k+1); 

         xnp = wx0*wz0*xn0i_c(i,k)+wx0*wz1*xn0i_c(i,k+1) 
                 +wx1*wz0*xn0i_c(i+1,k)+wx1*wz1*xn0i_c(i+1,k+1); 

         b=1-tor+tor*bfldp;


/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pitch angle collision!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

        ni_temp= wx0*wz0*xn0i_c(i,k)+wx0*wz1*xn0i_c(i,k+1) 
                 +wx1*wz0*xn0i_c(i+1,k)+wx1*wz1*xn0i_c(i+1,k+1); 
         energy0 = (mu_ptr[m]*b+0.5*mims_ptr[0]*pow(u3_ptr[m],2));
         energy =  max(energy0,0.1*T_center);
         if(icollision == 1){
            nudi=pow(nudi0*ni_temp/(energy), 1.5); 
            //  call random_number(rrr)
            //  p_m = int(2*rrr-1)
         }
         //write(*,*) rand_table(globle_integer-1), u2(m),u2(m)*(1-nudi*dt)+rand_table(globle_integer)*sqrt((2*energy0/mims(1)-u2(m)**2)*nudi*dt)!,mu(m),(energy0-0.5*mims(1)*u2(m)**2)/b
         u2_ptr[m]=u2_ptr[m]*(1-nudi*dt)+rand_table_ptr[globle_integer]*sqrt((2*energy0/mims_ptr[0]-pow(u2_ptr[m],2))*nudi*dt);
         globle_integer = (globle_integer+1) % 10007;
         //   write(*,*) rand_table(globle_integer-1), u2(m),u3(m)!,mu(m),(energy0-0.5*mims(1)*u2(m)**2)/b
         //    write(*,*) rand_table(globle_integer-1), mu(m),(energy0-0.5*mims(1)*u2(m)**2)/b
         mu_ptr[m]= (energy0-0.5*mims_ptr[0]*pow(u2_ptr[m],2))/b;

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end of pitch angle collision!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   */



        rhog=sqrt(2*b*mu_ptr[m]*mims_ptr[0])/(q_ptr[0]*b)*iflr;

        rhox[0] = rhog;
        rhoy[0] = 0;
        rhox[1] = -rhox[0];
        rhoy[1] = -rhoy[0];
        rhox[2] = 0;
        rhoy[2] = rhog;
        rhox[3] = 0;
        rhoy[3] = -rhoy[2];
    //calculate avg. e-field...
    //do 1,2,4 point average, where lr is the no. of points...

    exp1=0;
    ezp=0;
    ezetap=0;
    delbxp=0;
    delbzp=0;

//  4 pt. avg. done explicitly for vectorization...
#pragma acc loop seq
        for(l = 0; l < lr_ptr[0]; ++l){

            xt=x2_ptr[m]+rhox[l]; //rwx(1,l)*rhog
            zt=z2_ptr[m]+rhoy[l]; //(rwy(1,l)+sz*rwx(1,l))*rhog;
            //zeta=modulo(zeta2(m),pi2);
     
   //particle can go out of bounds during gyroavg...
            if( (xt<2*dxeq) || (xt>lx-2*dxeq) ) xt=x2_ptr[m];
            if( (zt<2*dzeq) || (zt>lz-2*dzeq) ) zt=z2_ptr[m];
            zeta=zeta2_ptr[m];
            //xt=modulo(xs,xdim)
            //zt=modulo(zt,zdim)
            i=int(xt/dx);
            j=int(zt/dz);
            k=int(zeta/dzeta);


            wx0=float(i+1)-xt/dx;
            wx1=1-wx0;
            wy0=float(j+1)-zt/dz;
            wy1=1-wy0;
            wz0=float(k+1)-zeta/dzeta;
            wz1=1-wz0;

                k_plus_1=k+1;
                if(k==kmx) k_plus_1=0;
            exp1=exp1 + wx0*wy0*wz0*ex_c(i,j,k) + wx1*wy0*wz0*ex_c(i+1,j,k) 
            + wx0*wy1*wz0*ex_c(i,j+1,k) + wx1*wy1*wz0*ex_c(i+1,j+1,k) + 
            wx0*wy0*wz1*ex_c(i,j, k_plus_1) + wx1*wy0*wz1*ex_c(i+1,j, k_plus_1) + 
            wx0*wy1*wz1*ex_c(i,j+1, k_plus_1) + wx1*wy1*wz1*ex_c(i+1,j+1, k_plus_1);

            ezp=ezp + wx0*wy0*wz0*ez_c(i,j,k) + wx1*wy0*wz0*ez_c(i+1,j,k) 
            + wx0*wy1*wz0*ez_c(i,j+1,k) + wx1*wy1*wz0*ez_c(i+1,j+1,k) + 
            wx0*wy0*wz1*ez_c(i,j, k_plus_1) + wx1*wy0*wz1*ez_c(i+1,j, k_plus_1) + 
            wx0*wy1*wz1*ez_c(i,j+1, k_plus_1) + wx1*wy1*wz1*ez_c(i+1,j+1, k_plus_1);

            ezetap =ezetap + wx0*wy0*wz0*ezeta_c(i,j,k) + wx1*wy0*wz0*ezeta_c(i+1,j,k) 
            + wx0*wy1*wz0*ezeta_c(i,j+1,k) + wx1*wy1*wz0*ezeta_c(i+1,j+1,k) + 
            wx0*wy0*wz1*ezeta_c(i,j, k_plus_1) + wx1*wy0*wz1*ezeta_c(i+1,j, k_plus_1) + 
            wx0*wy1*wz1*ezeta_c(i,j+1, k_plus_1) + wx1*wy1*wz1*ezeta_c(i+1,j+1, k_plus_1);

            delbxp =delbxp + wx0*wy0*wz0*delbx_c(i,j,k)  
            + wx1*wy0*wz0*delbx_c(i+1,j,k) 
            + wx0*wy1*wz0*delbx_c(i,j+1,k) 
            + wx1*wy1*wz0*delbx_c(i+1,j+1,k) 
            + wx0*wy0*wz1*delbx_c(i,j, k_plus_1) 
            + wx1*wy0*wz1*delbx_c(i+1,j, k_plus_1) 
            + wx0*wy1*wz1*delbx_c(i,j+1, k_plus_1) 
            + wx1*wy1*wz1*delbx_c(i+1,j+1, k_plus_1);

            delbzp =delbzp + wx0*wy0*wz0*delbz_c(i,j,k) 
            + wx1*wy0*wz0*delbz_c(i+1,j,k) 
            + wx0*wy1*wz0*delbz_c(i,j+1,k) 
            + wx1*wy1*wz0*delbz_c(i+1,j+1,k)  
            + wx0*wy0*wz1*delbz_c(i,j, k_plus_1) 
            + wx1*wy0*wz1*delbz_c(i+1,j, k_plus_1)  
            + wx0*wy1*wz1*delbz_c(i,j+1, k_plus_1)  
            + wx1*wy1*wz1*delbz_c(i+1,j+1, k_plus_1);
        }
        exp1 = exp1/4;
        ezp = ezp/4;
        ezetap = ezetap/4;
        delbxp = delbxp/4;
        delbzp = delbzp/4;

         vfac = 0.5*(mims_ptr[0]*pow(u2_ptr[m],2) + 2*mu_ptr[m]*b);
         kapxp = kapnxp - (1.5-vfac/ter)*kaptxp;
         kapzp = kapnzp - (1.5-vfac/ter)*kaptzp;        

         vpar = u2_ptr[m];
         enerb=(mu_ptr[m]+mims_ptr[0]*vpar*vpar/b)/q_ptr[0]*tor;

         Bstar3[0]=bfldxp +mims_ptr[0]*vpar*curlbp[0]/q_ptr[0]+delbxp;
         Bstar3[1]=bfldzp+mims_ptr[0]*vpar*curlbp[1]/q_ptr[0]+delbzp;
         Bstar3[2]=bfldzetap+mims_ptr[0]*vpar*curlbp[2]/q_ptr[0];


         //bstar=b+mims(1)*vpar*bdcurlbp/q(1)
         bstar=(bfldxp*Bstar3[0]+bfldzp*Bstar3[1]+bfldzetap*Bstar3[2])/bfldp;
         //  write(*,*)bstar-b-mims(1)*vpar*bdcurlbp/q(1), b-bfldp
         //  vcurlbdotE=vpar*(exp1*curlbp(1)+ezp*curlbp(2)+ezetap*curlbp(3))
         
         dum1 = 1;
         // vxdum = (ezp/b+vpar/b*delbxp)*dum1
         //  vxdum = (ezp*bfldzetap-ezetap*bfldzp)/b**2
         // xdot = vxdum*nonlin +vpar*bfldxp/b-enerb/bfldp/bfldp*bfldzetap*dbdzp
         // vzdum = (ezetap*bfldxp-exp1*bfldzetap)/b**2
         //  zdot = (-exp1/b+vpar/b*delbzp)*dum1*nonlin &
         //      +vpar*bfldzp/b+enerb/bfldp/bfldp*bfldzetap*dbdxp
         //  zdot = vzdum*nonlin+vpar*bfldzp/b+enerb/bfldp/bfldp*bfldzetap*dbdxp

         //  vzetadum= (exp1*bfldzp-ezp*bfldxp)/b**2
         //  zetadot = vzetadum/x*nonlin + vpar*bfldzetap/(x*b)+enerb/(x*b*b)*(bfldxp*dbdzp-bfldzp*dbdxp)



//         write(*,*) dbdzetap
         xdot = (vpar*Bstar3[0]+(mu_ptr[m]*(bfldzp*dbdzetap-bfldzetap*dbdzp)/q_ptr[0]+(ezp*bfldzetap-ezetap*bfldzp))/(b))/bstar;
         zdot = (vpar*Bstar3[1]+(mu_ptr[m]*(bfldzetap*dbdxp-bfldxp*dbdzetap)/q_ptr[0]+ (ezetap*bfldxp-exp1*bfldzetap))/(b))/bstar;
         zetadot = (vpar*Bstar3[2]+(mu_ptr[m]*(bfldxp*dbdzp-bfldzp*dbdxp)/q_ptr[0]+(exp1*bfldzp-ezp*bfldxp))/(b))/bstar;


         
          //pzd0 = -mu(m)/mims(1)/b*(bfldxp*dbdxp+bfldzp*dbdzp)
          //write(*,*)
          //write(*,*)(exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)/b*q(1)/mims(1)
          //pzdot = pzd0+(exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)/b*(q(1)/mims(1)+bdcurlbp*vpar/b)*nonlin

          //pzdot = pzd0+((exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)*q(1)/mims(1)+vcurlbdotE)/bstar*nonlin


         pzdot = (Bstar3[0]*(q_ptr[0]*exp1-mu_ptr[m]*dbdxp)+Bstar3[1]*(q_ptr[0]*ezp-mu_ptr[m]*dbdzp)+Bstar3[2]*(q_ptr[0]*ezetap-mu_ptr[m]*dbdzetap))/(mims_ptr[0]*bstar);
          
         
         edot = q_ptr[0]*(xdot*exp1+zdot*ezp+zetadot*ezetap);

         x3_ptr[m] = x2_ptr[m] + 0.5*dt*xdot;
         z3_ptr[m] = z2_ptr[m] + 0.5*dt*zdot;
         zeta3_ptr[m] = zeta2_ptr[m] + 0.5*dt*zetadot;
         u3_ptr[m] = u2_ptr[m] + 0.5*dt*pzdot;

         //dum = 1.0
         //vxdum = (ezp/b+vpar/b*delbxp)*dum1
         //vzdum = (-exp1/b+vpar/b*delbzp)*dum1
         //vxdum = eyp+vpar/b*delbxp
         //w3(m)=w2(m) + 0.5*dt*(vxdum*kapxp + vzdum*kapzp+edot/ter)*dum*xnp

        if( !((x3_ptr[m]>2*dxeq) && (x3_ptr[m]<lx-2*dxeq) && (z3_ptr[m]>2*dzeq) && (z3_ptr[m]<lz-2*dzeq)) ) 
        {
          u3_ptr[m]=u2_ptr[m];
          x3_ptr[m]=x2_ptr[m];
          z3_ptr[m]=z2_ptr[m];
          zeta3_ptr[m]=zeta2_ptr[m];
          w3_ptr[m]=0;
        }
    }
    #pragma acc wait
    end_ppush_tm = MPI_Wtime();
    ppush_tm = ppush_tm + end_ppush_tm - start_ppush_tm;
}

//!-------------- End of subroutine ppush --------------------------------

void cpush_c_(int &n){
        double exp1,ezp,ezetap,delbxp,delbzp,nudi0,nudi=0,ni_temp,energy,rrr,T_center,energy0;
        double wx0,wx1,wy0,wy1,wz0,wz1,dum,vxdum,vzdum,dum1,vzetadum;
        int m,i,j,k,l,k_plus_1=0,p_m;
        double rhog,vfac,kapxp,kapzp,vpar,pidum,kaptxp,kapnxp,kaptzp,kapnzp,xnp,bdcurlbp;
        double b,th,r,enerb,qr,ter,x,z,zeta;
        double xt,xs,zt,xdot,zdot,zetadot,xdt,ydt,pzdot,edot,pzd0,vp0,vcurlbdotE;
        double dbdxp,dbdzp,bfldp,bfldxp,bfldzp,bfldzetap, bstar, dbdzetap=0;
        double rhox[4], rhoy[4], curlbp[3], Bstar3[3];

        //real(8),dimension(3)::curlbp,Bstar3
        start_cpush_tm = MPI_Wtime();
        nudi0 = 1/sqrt(2)*18.4*pow(e,1.5)*(4.7140*pow(10,-8))*(1*pow(10,-6));
        //write(*,*)t0i(200,201)
        T_center = t0i_c(imx/2,jmx/2);
        //write(*,*)T_center*e
        //write(*,*)nudi0
//#pragma parallel loop gang vector private(bstar3,rhoy,rhox) copy(rand_table_ptr)
        for(m = 0; m < mm_ptr[0]; ++m){
            x=x3_ptr[m];
            i = static_cast<int>(x/dxeq);
            i = min(i,nx-1);
            wx0 = (i+1)-x/dxeq;
            wx1 = 1-wx0;

            z = z3_ptr[m];
            k = static_cast<int>(z/dzeq);
            k = min(k,nz-1);
            wz0 = (k+1)-z/dzeq;
            wz1 = 1-wz0;

            ni_temp= wx0*wz0*xn0i_c(i,k)+wx0*wz1*xn0i_c(i,k+1) 
                 +wx1*wz0*xn0i_c(i+1,k)+wx1*wz1*xn0i_c(i+1,k+1); 
           //write(*,*) xn0i(i,k)


           // bdcurlbp =wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
                             //+wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1)
            for(j = 0; j <= 2; ++j){
                curlbp[j]= wx0*wz0*curlb_c(i,k,j)+wx0*wz1*curlb_c(i,k+1,j) 
                +wx1*wz0*curlb_c(i+1,k,j)+wx1*wz1*curlb_c(i+1,k+1,j);
            }

         dbdxp = wx0*wz0*dbdx_c(i,k)+wx0*wz1*dbdx_c(i,k+1) 
                 +wx1*wz0*dbdx_c(i+1,k)+wx1*wz1*dbdx_c(i+1,k+1); 
         dbdzp = wx0*wz0*dbdz_c(i,k)+wx0*wz1*dbdz_c(i,k+1) 
                 +wx1*wz0*dbdz_c(i+1,k)+wx1*wz1*dbdz_c(i+1,k+1);
         bfldp = wx0*wz0*b0_c(i,k)+wx0*wz1*b0_c(i,k+1) 
                 +wx1*wz0*b0_c(i+1,k)+wx1*wz1*b0_c(i+1,k+1); 
         bfldxp = wx0*wz0*b0x_c(i,k)+wx0*wz1*b0x_c(i,k+1) 
                 +wx1*wz0*b0x_c(i+1,k)+wx1*wz1*b0x_c(i+1,k+1); 
         bfldzp = wx0*wz0*b0z_c(i,k)+wx0*wz1*b0z_c(i,k+1) 
                 +wx1*wz0*b0z_c(i+1,k)+wx1*wz1*b0z_c(i+1,k+1); 
         bfldzetap = wx0*wz0*b0zeta_c(i,k)+wx0*wz1*b0zeta_c(i,k+1) 
                 +wx1*wz0*b0zeta_c(i+1,k)+wx1*wz1*b0zeta_c(i+1,k+1);
         ter = wx0*wz0*t0i_c(i,k)+wx0*wz1*t0i_c(i,k+1) 
                 +wx1*wz0*t0i_c(i+1,k)+wx1*wz1*t0i_c(i+1,k+1); 
         kaptxp = wx0*wz0*captix_c(i,k)+wx0*wz1*captix_c(i,k+1) 
                 +wx1*wz0*captix_c(i+1,k)+wx1*wz1*captix_c(i+1,k+1);
         kapnxp = wx0*wz0*capnix_c(i,k)+wx0*wz1*capnix_c(i,k+1)
                 +wx1*wz0*capnix_c(i+1,k)+wx1*wz1*capnix_c(i+1,k+1); 

         kaptzp = wx0*wz0*captiz_c(i,k)+wx0*wz1*captiz_c(i,k+1) 
                 +wx1*wz0*captiz_c(i+1,k)+wx1*wz1*captiz_c(i+1,k+1); 
         kapnzp = wx0*wz0*capniz_c(i,k)+wx0*wz1*capniz_c(i,k+1) 
                 +wx1*wz0*capniz_c(i+1,k)+wx1*wz1*capniz_c(i+1,k+1); 

         xnp = wx0*wz0*xn0i_c(i,k)+wx0*wz1*xn0i_c(i,k+1) 
                 +wx1*wz0*xn0i_c(i+1,k)+wx1*wz1*xn0i_c(i+1,k+1);

         b=1-tor+tor*bfldp;

         rhog=sqrt(2*b*mu_ptr[m]*mims_ptr[0])/(q_ptr[0]*b)*iflr;

         rhox[0] = rhog;
         rhoy[0] = 0;
         rhox[1] = -rhox[0];
         rhoy[1] = -rhoy[0];
         rhox[2] = 0;
         rhoy[2] = rhog;
         rhox[3] = 0;
         rhoy[3] = -rhoy[2];
//    calculate avg. e-field...
//    do 1,2,4 point average, where lr is the no. of points...

         exp1=0;
         ezp=0;
         ezetap=0;
         delbxp=0;
         delbzp=0;
        #pragma loop seq
        for(l = 0; l < lr_ptr[0]; ++l){
//SP            xs=x3(m)+rhox(l) !rwx(1,l)*rhog
            xt=x3_ptr[m]+rhox[l]; //rwx(1,l)*rhog
            zt=z3_ptr[m]+rhoy[l]; //(rwy(1,l)+sz*rwx(1,l))*rhog
//
//   particle can go out of bounds during gyroavg...
            if( (xt<2*dxeq)||(xt>lx-2*dxeq) ) xt=x3_ptr[m];
            if( (zt<2*dzeq)||(zt>lz-2*dzeq) ) zt=z3_ptr[m];
            zeta= fmod(zeta3_ptr[m], 2*M_PI);
            i=int(xt/dx);
            j=int(zt/dz);
            k=int(zeta/dzeta);


            wx0=float(i+1)-xt/dx;
            wx1=1-wx0;
            wy0=float(j+1)-zt/dz;
            wy1=1-wy0;
            wz0=float(k+1)-zeta/dzeta;
            wz1=1-wz0;

                k_plus_1=k+1;
                if(k==kmx) k_plus_1=0;
            exp1=exp1 + wx0*wy0*wz0*ex_c(i,j,k) + wx1*wy0*wz0*ex_c(i+1,j,k) 
            + wx0*wy1*wz0*ex_c(i,j+1,k) + wx1*wy1*wz0*ex_c(i+1,j+1,k) + 
            wx0*wy0*wz1*ex_c(i,j, k_plus_1) + wx1*wy0*wz1*ex_c(i+1,j, k_plus_1) + 
            wx0*wy1*wz1*ex_c(i,j+1, k_plus_1) + wx1*wy1*wz1*ex_c(i+1,j+1, k_plus_1);

            ezp=ezp + wx0*wy0*wz0*ez_c(i,j,k) + wx1*wy0*wz0*ez_c(i+1,j,k) 
            + wx0*wy1*wz0*ez_c(i,j+1,k) + wx1*wy1*wz0*ez_c(i+1,j+1,k) + 
            wx0*wy0*wz1*ez_c(i,j, k_plus_1) + wx1*wy0*wz1*ez_c(i+1,j, k_plus_1) + 
            wx0*wy1*wz1*ez_c(i,j+1, k_plus_1) + wx1*wy1*wz1*ez_c(i+1,j+1, k_plus_1);

            ezetap =ezetap + wx0*wy0*wz0*ezeta_c(i,j,k) + wx1*wy0*wz0*ezeta_c(i+1,j,k) 
            + wx0*wy1*wz0*ezeta_c(i,j+1,k) + wx1*wy1*wz0*ezeta_c(i+1,j+1,k) + 
            wx0*wy0*wz1*ezeta_c(i,j, k_plus_1) + wx1*wy0*wz1*ezeta_c(i+1,j, k_plus_1) + 
            wx0*wy1*wz1*ezeta_c(i,j+1, k_plus_1) + wx1*wy1*wz1*ezeta_c(i+1,j+1, k_plus_1);

            delbxp =delbxp + wx0*wy0*wz0*delbx_c(i,j,k)  
            + wx1*wy0*wz0*delbx_c(i+1,j,k) 
            + wx0*wy1*wz0*delbx_c(i,j+1,k) 
            + wx1*wy1*wz0*delbx_c(i+1,j+1,k) 
            + wx0*wy0*wz1*delbx_c(i,j, k_plus_1) 
            + wx1*wy0*wz1*delbx_c(i+1,j, k_plus_1) 
            + wx0*wy1*wz1*delbx_c(i,j+1, k_plus_1) 
            + wx1*wy1*wz1*delbx_c(i+1,j+1, k_plus_1);

            delbzp =delbzp + wx0*wy0*wz0*delbz_c(i,j,k) 
            + wx1*wy0*wz0*delbz_c(i+1,j,k) 
            + wx0*wy1*wz0*delbz_c(i,j+1,k) 
            + wx1*wy1*wz0*delbz_c(i+1,j+1,k)  
            + wx0*wy0*wz1*delbz_c(i,j, k_plus_1)  
            + wx1*wy0*wz1*delbz_c(i+1,j, k_plus_1)  
            + wx0*wy1*wz1*delbz_c(i,j+1, k_plus_1)  
            + wx1*wy1*wz1*delbz_c(i+1,j+1, k_plus_1);
        }
         exp1 = exp1/4;
         ezp = ezp/4;
         ezetap = ezetap/4;
         delbxp = delbxp/4;
         delbzp = delbzp/4;
       
         vfac = 0.5*(mims_ptr[0]*pow(u2_ptr[m],2) + 2*mu_ptr[m]*b);
         kapxp = kapnxp - (1.5-vfac/ter)*kaptxp;
         kapzp = kapnzp - (1.5-vfac/ter)*kaptzp;        
       
         vpar = u3_ptr[m];
         enerb=(mu_ptr[m]+mims_ptr[0]*vpar*vpar/b)/q_ptr[0]*tor;
       
         Bstar3[0]=bfldxp+mims_ptr[0]*vpar*curlbp[0]/q_ptr[0]+delbxp;
         Bstar3[1]=bfldzp+mims_ptr[0]*vpar*curlbp[0]/q_ptr[0]+delbzp;
         Bstar3[2]=bfldzetap+mims_ptr[0]*vpar*curlbp[2]/q_ptr[0];
       
                
       //bstar=b+mims(1)*vpar*bdcurlbp/q(1)
       
         bstar=(bfldxp*Bstar3[0]+bfldzp*Bstar3[1]+bfldzetap*Bstar3[2])/bfldp;


        //vcurlbdotE=vpar*(exp1*curlbp(1)+ezp*curlbp(2)+ezetap*curlbp(3))


         dum1 = 1;
//         !         vxdum = (ezp/b+vpar/b*delbxp)*dum1
// !         vxdum =(ezp*bfldzetap-ezetap*bfldzp)/b**2
// !         xdot = vxdum*nonlin +vpar*bfldxp/b-enerb/bfldp/bfldp*bfldzetap*dbdzp
//          !         write(*,*) vxdum*nonlin, (-exp1/b+vpar/b*delbzp)*dum1*nonlin
// !         vzdum =(ezetap*bfldxp-exp1*bfldzetap)/b**2
//  !        zdot = (-exp1/b+vpar/b*delbzp)*dum1*nonlin &
//          !            +vpar*bfldzp/b+enerb/bfldp/bfldp*bfldzetap*dbdxp
// !         write(*,*)vzdum,vxdum
// !         zdot = vzdum*nonlin+vpar*bfldzp/b+enerb/bfldp/bfldp*bfldzetap*dbdxp

// !         zetadot =  vpar*bfldzetap/(x*b)+enerb/(x*b*b)*(bfldxp*dbdzp-bfldzp*dbdxp)

// !         vzetadum= (exp1*bfldzp-ezp*bfldxp)/b**2
// !         zetadot = vzetadum/x*nonlin + vpar*bfldzetap/(x*b)+enerb/(x*b*b)*(bfldxp*dbdzp-bfldzp*dbdxp)


// !         write(*,*)dbdzetap
         xdot = (vpar*Bstar3[0]+(mu_ptr[m]*(bfldzp*dbdzetap-bfldzetap*dbdzp)/q_ptr[0]+(ezp*bfldzetap-ezetap*bfldzp))/(b))/bstar;
         zdot = (vpar*Bstar3[1]+(mu_ptr[m]*(bfldzetap*dbdxp-bfldxp*dbdzetap)/q_ptr[0]+ (ezetap*bfldxp-exp1*bfldzetap))/(b))/bstar;
         zetadot = (vpar*Bstar3[2]+(mu_ptr[m]*(bfldxp*dbdzp-bfldzp*dbdxp)/q_ptr[0]+(exp1*bfldzp-ezp*bfldxp))/(b))/bstar;



// !        pzd0 = -mu(m)/mims(1)/b*(bfldxp*dbdxp+bfldzp*dbdzp)
// !         pzdot = pzd0+(exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)/b*q(1)/mims(1)*nonlin
// !         pzdot = pzd0+(exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)/b*(q(1)/mims(1)+bdcurlbp*vpar/b)*nonlin
// !          pzdot = pzd0+((exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)*q(1)/mims(1)+vcurlbdotE)/bstar*nonlin


         pzdot = (Bstar3[0]*(q_ptr[0]*exp1-mu_ptr[m]*dbdxp)+Bstar3[1]*(q_ptr[0]*ezp-mu_ptr[m]*dbdzp)+Bstar3[2]*(q_ptr[0]*ezetap-mu_ptr[m]*dbdzetap))/(mims_ptr[0]*bstar);
 

         edot = q_ptr[0]*(xdot*exp1+zdot*ezp+zetadot*ezetap);

         x3_ptr[m] = x2_ptr[m] + dt*xdot;
         z3_ptr[m] = z2_ptr[m] + dt*zdot;
         zeta3_ptr[m] = zeta2_ptr[m] + dt*zetadot;
         u3_ptr[m] = u2_ptr[m] + dt*pzdot;

// !         dum = 1.0
// !         vxdum = (ezp/b+vpar/b*delbxp)*dum1
// !         vzdum = (-exp1/b+vpar/b*delbzp)*dum1
// !         vxdum = eyp+vpar/b*delbxp
// !         w3(m)=w2(m) + dt*(vxdum*kapxp + vzdum*kapzp+edot/ter)*dum*xnp

         zeta3_ptr[m]= fmod(zeta3_ptr[m],2*M_PI);


//         write(*,*)energy, nudi
if( (x3_ptr[m]>2*dxeq)&&(x3_ptr[m]<lx-2*dxeq)&&(z3_ptr[m]>2*dzeq)&&(z3_ptr[m]<lz-2*dzeq) ){
            //energy0 = (mu(m)*b+0.5*mims(1)*u3(m)**2)
            //energy =  max(energy0,0.1*T_center)
            //if (icollision==1) then
            //      nudi=nudi0*ni_temp/(energy)**1.5
            //           call random_number(rrr)
            //           p_m = int(2*rrr-1)
            //end if
            u2_ptr[m]=u3_ptr[m]; //*(1-nudi*dt)+rand_table(globle_integer)*sqrt((2*energy0/mims(1)-u3(m)**2)*nudi*dt)
            //  !$acc atomic
            //  globle_integer = mod(globle_integer+1,10007)
               
            //write(*,*) rand_table(globle_integer-1), u2(m),u3(m)
            x2_ptr[m]=x3_ptr[m];
            z2_ptr[m]=z3_ptr[m];
            zeta2_ptr[m]=zeta3_ptr[m];
            w2_ptr[m]=w3_ptr[m];
         
            // write(*,*)rand_table(globle_integer-1),mu(m),(energy0-0.5*mims(1)*u3(m)**2)/b
            // mu(m)= (energy0-0.5*mims(1)*u2(m)**2)/b
        } else {
          u3_ptr[m]=u2_ptr[m];
          x3_ptr[m]=x2_ptr[m];
          z3_ptr[m]=z2_ptr[m];
          zeta3_ptr[m]=zeta2_ptr[m];
          w2_ptr[m]=0;
          w3_ptr[m]=0;
        }
        // if(Myid==0)then
        //   open(935, file='flag_debug',status='unknown',position='append')
        //   write(935,*)'before PETSc sovling'
        //   close(935)
        // end if
        //       if(myid==0 .and. m==1)then
        //            open(93, file='test_energy',status='unknown',position='append')
        //            write(93,*) mu(m)*b+0.5*mims(1)*u2(m)**2-q(1)*z2(m)*ez(i,j,k)
        //            close(93)
        //                write(*,*) mu(m)*b+0.5*mims(1)*u2(m)**2-q(1)*z2(m)*ez(i,j,k)
        //          close(19)
        //       end if  
  }
  #pragma acc wait
  end_cpush_tm = MPI_Wtime();
  cpush_tm = cpush_tm + end_cpush_tm - start_cpush_tm;
}