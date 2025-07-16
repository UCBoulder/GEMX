#include "ionPush_c.hpp" 

using namespace std;

#pragma acc routine seq
inline double my_fmod(double x, double y) {
    return x - y * floor(x / y);
}


//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//       Ion pre-push
//
void ppush_c_(const int &n) { 
    double exp1 = 0,ezp = 0,ezetap = 0,delbxp = 0,delbzp  = 0,energy = 0, energy0 = 0,nudi0 = 0,nudi = 0,T_center = 0,ni_temp = 0;
    double wx0 = 0,wx1 = 0,wy0 = 0,wy1 = 0,wz0 = 0,wz1 = 0,dum1 = 0; 
    int m = 0,i = 0,j = 0,k = 0,l = 0,k_plus_1 = 0;
    double rhog = 0,vfac = 0,kapxp = 0,kapzp = 0,vpar = 0,kaptxp = 0,kapnxp = 0,kaptzp = 0,kapnzp = 0,xnp = 0;
    double b = 0,enerb = 0,ter = 0,z = 0,zeta = 0,bstar = 0;
    double x = 0;
    double xt = 0,zt = 0,xdot = 0,zdot = 0,zetadot = 0,pzdot = 0,edot = 0;
    double dbdxp = 0,dbdzp = 0,bfldp = 0,bfldxp = 0,bfldzp = 0,bfldzetap = 0,dbdzetap=0;
    double rhox[4], rhoy[4], BStar3[3], curlbp[3];
//     fill(rhox,  rhox+3, 0.0);
//     fill(rhoy,  rhoy+3, 0.0);
//     fill(curlbp,  curlbp+2, 0.0);
//     fill(BStar3,  BStar3+2, 0.0);
    //real(8),dimension(3)::curlbp,BStar3
    start_ppush_tm = MPI_Wtime();

    nudi0 = 1/sqrt(2.0)*18.4*pow(e,1.5)* 4.7140e-8* 1.e-6;
    //write(*,*)t0i(200,201);
    T_center = t0i(imx/2, jmx/2);

    prepareDeviceData();
    #pragma acc data \
    copyin(x2[0:mmx], z2[0:mmx], zeta2[0:mmx], rand_table[0:10007]) \
    copy(mu[0:mmx], u2[0:mmx], u3[0:mmx], x3[0:mmx], z3[0:mmx], zeta3[0:mmx], w3[0:mmx])
    #pragma acc parallel loop gang vector private(rhoy,BStar3,rhox)
    for(m = 0; m < mm[0]; ++m){
        x = x2[m];
        i = static_cast<int>(x/dxeq);
        i = min(i,nx-1);
        wx0 = (i+1)-x/dxeq;
        wx1 = 1.-wx0;

        z = z2[m];
        k = static_cast<int>(z/dzeq);
        k = min(k,nz-1);
        wz0 = (k+1)-z/dzeq;
        wz1 = 1.-wz0;

        //bdcurlbp =wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
        //                 +wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1)
        #pragma loop seq
        for(j = 0; j <= 2; ++j){
            curlbp[j]= wx0*wz0*curlb(i,k,j)+wx0*wz1*curlb(i,k+1,j) 
              +wx1*wz0*curlb(i+1,k,j)+wx1*wz1*curlb(i+1,k+1,j);
        }
        //write(*,*) curlbp(1)-wx0*wz0*curlb(i,k,1)-wx0*wz1*curlb(i,k+1,1)-wx1*wz0*curlb(i+1,k,1)-wx1*wz1*curlb(i+1,k+1,1),curlbp(2)-wx0*wz0*curlb(i,k,2)-wx0*wz1*curlb(i,k+1,2)-wx1*wz0*curlb(i+1,k,2)-wx1*wz1*curlb(i+1,k+1,2),curlbp(3)-wx0*wz0*curlb(i,k,3)-wx0*wz1*curlb(i,k+1,3)-wx1*wz0*curlb(i+1,k,3)-wx1*wz1*curlb(i+1,k+1,3)
        //         write(*,*)curlbp(1),curlbp(2),curlbp(3)
        // write(*,*)bdcurlbp
         dbdxp = wx0*wz0*dbdx(i,k)+wx0*wz1*dbdx(i,k+1) 
                 +wx1*wz0*dbdx(i+1,k)+wx1*wz1*dbdx(i+1,k+1); 
         dbdzp = wx0*wz0*dbdz(i,k)+wx0*wz1*dbdz(i,k+1) 
                 +wx1*wz0*dbdz(i+1,k)+wx1*wz1*dbdz(i+1,k+1);
         bfldp = wx0*wz0*b0(i,k)+wx0*wz1*b0(i,k+1) 
                 +wx1*wz0*b0(i+1,k)+wx1*wz1*b0(i+1,k+1); 
         bfldxp = wx0*wz0*b0x(i,k)+wx0*wz1*b0x(i,k+1) 
                 +wx1*wz0*b0x(i+1,k)+wx1*wz1*b0x(i+1,k+1); 
         bfldzp = wx0*wz0*b0z(i,k)+wx0*wz1*b0z(i,k+1) 
                 +wx1*wz0*b0z(i+1,k)+wx1*wz1*b0z(i+1,k+1); 
         bfldzetap = wx0*wz0*b0zeta(i,k)+wx0*wz1*b0zeta(i,k+1) 
                 +wx1*wz0*b0zeta(i+1,k)+wx1*wz1*b0zeta(i+1,k+1); 
         ter = wx0*wz0*t0i(i,k)+wx0*wz1*t0i(i,k+1) 
                 +wx1*wz0*t0i(i+1,k)+wx1*wz1*t0i(i+1,k+1); 
         kaptxp = wx0*wz0*captix(i,k)+wx0*wz1*captix(i,k+1) 
                 +wx1*wz0*captix(i+1,k)+wx1*wz1*captix(i+1,k+1);

         kapnxp = wx0*wz0*capnix(i,k)+wx0*wz1*capnix(i,k+1) 
                 +wx1*wz0*capnix(i+1,k)+wx1*wz1*capnix(i+1,k+1); 

         kaptzp = wx0*wz0*captiz(i,k)+wx0*wz1*captiz(i,k+1) 
                 +wx1*wz0*captiz(i+1,k)+wx1*wz1*captiz(i+1,k+1); 
         kapnzp = wx0*wz0*capniz(i,k)+wx0*wz1*capniz(i,k+1) 
                 +wx1*wz0*capniz(i+1,k)+wx1*wz1*capniz(i+1,k+1); 

         xnp = wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) 
                 +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1); 

         b=1.-tor+tor*bfldp;


/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pitch angle collision!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

        // ni_temp= wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) //repeated code here
        //          +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1); 
         energy0 = (mu[m]*b+0.5*mims[0]*(u3[m] * u3[m]));
         energy =  max(energy0,0.1*T_center);
         if(icollision == 1){
            nudi=nudi0*xnp/pow(energy, 1.5);  //whole thing to 1.5 or just energy?
            //  call random_number(rrr)
            //  p_m = int(2*rrr-1)
         }
         //write(*,*) rand_table(globle_integer-1), u2(m),u2(m)*(1-nudi*dt)+rand_table(globle_integer)*sqrt((2*energy0/mims(1)-u2(m)**2)*nudi*dt)!,mu(m),(energy0-0.5*mims(1)*u2(m)**2)/b
       
         u2[m]=u2[m]*(1-nudi*dt)+rand_table[globle_integer]*sqrt((2*energy0/mims[0]-(u2[m]*u2[m]))*nudi*dt);
         globle_integer = (globle_integer+1) % 10007;
         //   write(*,*) rand_table(globle_integer-1), u2(m),u3(m)!,mu(m),(energy0-0.5*mims(1)*u2(m)**2)/b
         //    write(*,*) rand_table(globle_integer-1), mu(m),(energy0-0.5*mims(1)*u2(m)**2)/b

         mu[m]= (energy0-0.5*mims[0]*(u2[m] * u2[m]))/b;

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end of pitch angle collision!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   */



        rhog=sqrt(2.*b*mu[m]*mims[0])/(q[0]*b)*iflr;

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
        for(l = 0; l <= lr[0]; ++l){

            xt=x2[m]+rhox[l]; //rwx(1,l)*rhog
            zt=z2[m]+rhoy[l]; //(rwy(1,l)+sz*rwx(1,l))*rhog;
            //zeta=modulo(zeta2(m),pi2);
     
   //particle can go out of bounds during gyroavg...
            if( (xt<2*dxeq) || (xt>lx-2*dxeq) ) {xt=x2[m];}
            if( (zt<2*dzeq) || (zt>lz-2*dzeq) ) {zt=z2[m];}
            zeta=zeta2[m];
            //xt=modulo(xs,xdim)
            //zt=modulo(zt,zdim)
            i=static_cast<int>(xt/dx);
            j=static_cast<int>(zt/dz);
            k=static_cast<int>(zeta/dzeta);


            wx0=static_cast<float>(i+1)-xt/dx;
            wx1=1.-wx0;
            wy0=static_cast<float>(j+1)-zt/dz;
            wy1=1.-wy0;
            wz0=static_cast<float>(k+1)-zeta/dzeta;
            wz1=1.-wz0;

                k_plus_1=k+1;
                if(k==kmx) k_plus_1=0;
            exp1=exp1 + wx0*wy0*wz0*ex(i,j,k) + wx1*wy0*wz0*ex(i+1,j,k) 
            + wx0*wy1*wz0*ex(i,j+1,k) + wx1*wy1*wz0*ex(i+1,j+1,k) + 
            wx0*wy0*wz1*ex(i,j, k_plus_1) + wx1*wy0*wz1*ex(i+1,j, k_plus_1) + 
            wx0*wy1*wz1*ex(i,j+1, k_plus_1) + wx1*wy1*wz1*ex(i+1,j+1, k_plus_1);

            ezp=ezp + wx0*wy0*wz0*ez(i,j,k) + wx1*wy0*wz0*ez(i+1,j,k) 
            + wx0*wy1*wz0*ez(i,j+1,k) + wx1*wy1*wz0*ez(i+1,j+1,k) + 
            wx0*wy0*wz1*ez(i,j, k_plus_1) + wx1*wy0*wz1*ez(i+1,j, k_plus_1) + 
            wx0*wy1*wz1*ez(i,j+1, k_plus_1) + wx1*wy1*wz1*ez(i+1,j+1, k_plus_1);

            ezetap =ezetap + wx0*wy0*wz0*ezeta(i,j,k) + wx1*wy0*wz0*ezeta(i+1,j,k) 
            + wx0*wy1*wz0*ezeta(i,j+1,k) + wx1*wy1*wz0*ezeta(i+1,j+1,k) + 
            wx0*wy0*wz1*ezeta(i,j, k_plus_1) + wx1*wy0*wz1*ezeta(i+1,j, k_plus_1) + 
            wx0*wy1*wz1*ezeta(i,j+1, k_plus_1) + wx1*wy1*wz1*ezeta(i+1,j+1, k_plus_1);

            delbxp =delbxp + wx0*wy0*wz0*delbx(i,j,k)  
            + wx1*wy0*wz0*delbx(i+1,j,k) 
            + wx0*wy1*wz0*delbx(i,j+1,k) 
            + wx1*wy1*wz0*delbx(i+1,j+1,k) 
            + wx0*wy0*wz1*delbx(i,j, k_plus_1) 
            + wx1*wy0*wz1*delbx(i+1,j, k_plus_1) 
            + wx0*wy1*wz1*delbx(i,j+1, k_plus_1) 
            + wx1*wy1*wz1*delbx(i+1,j+1, k_plus_1);

            delbzp =delbzp + wx0*wy0*wz0*delbz(i,j,k) 
            + wx1*wy0*wz0*delbz(i+1,j,k) 
            + wx0*wy1*wz0*delbz(i,j+1,k) 
            + wx1*wy1*wz0*delbz(i+1,j+1,k)  
            + wx0*wy0*wz1*delbz(i,j, k_plus_1) 
            + wx1*wy0*wz1*delbz(i+1,j, k_plus_1)  
            + wx0*wy1*wz1*delbz(i,j+1, k_plus_1)  
            + wx1*wy1*wz1*delbz(i+1,j+1, k_plus_1);
        }
        exp1 = exp1/4;
        ezp = ezp/4;
        ezetap = ezetap/4;
        delbxp = delbxp/4;
        delbzp = delbzp/4;

         vfac = 0.5*(mims[0] * u2[m] * u2[m] + 2 * mu[m] * b);
         kapxp = kapnxp - (1.5-vfac/ter)*kaptxp;
         kapzp = kapnzp - (1.5-vfac/ter)*kaptzp;        

         vpar = u2[m];
         enerb=(mu[m]+mims[0]*vpar*vpar/b)/q[0]*tor;

         BStar3[0]=bfldxp +mims[0]*vpar*curlbp[0]/q[0]+delbxp;
         BStar3[1]=bfldzp+mims[0]*vpar*curlbp[1]/q[0]+delbzp;
         BStar3[2]=bfldzetap+mims[0]*vpar*curlbp[2]/q[0];


         //bstar=b+mims(1)*vpar*bdcurlbp/q(1)
         bstar=(bfldxp*BStar3[0]+bfldzp*BStar3[1]+bfldzetap*BStar3[2])/bfldp;
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
         xdot = (vpar*BStar3[0]+(mu[m]*(bfldzp*dbdzetap-bfldzetap*dbdzp)/q[0]+(ezp*bfldzetap-ezetap*bfldzp))/(b))/bstar;
         zdot = (vpar*BStar3[1]+(mu[m]*(bfldzetap*dbdxp-bfldxp*dbdzetap)/q[0]+ (ezetap*bfldxp-exp1*bfldzetap))/(b))/bstar;
         zetadot = (vpar*BStar3[2]+(mu[m]*(bfldxp*dbdzp-bfldzp*dbdxp)/q[0]+(exp1*bfldzp-ezp*bfldxp))/(b))/bstar;


         
          //pzd0 = -mu(m)/mims(1)/b*(bfldxp*dbdxp+bfldzp*dbdzp)
          //write(*,*)
          //write(*,*)(exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)/b*q(1)/mims(1)
          //pzdot = pzd0+(exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)/b*(q(1)/mims(1)+bdcurlbp*vpar/b)*nonlin

          //pzdot = pzd0+((exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)*q(1)/mims(1)+vcurlbdotE)/bstar*nonlin


         pzdot = (BStar3[0]*(q[0]*exp1-mu[m]*dbdxp)+BStar3[1]*(q[0]*ezp-mu[m]*dbdzp)+BStar3[2]*(q[0]*ezetap-mu[m]*dbdzetap))/(mims[0]*bstar);
          
         
         edot = q[0]*(xdot*exp1+zdot*ezp+zetadot*ezetap);

         x3[m] = x2[m] + 0.5*dt*xdot;
         z3[m] = z2[m] + 0.5*dt*zdot;
         zeta3[m] = zeta2[m] + 0.5*dt*zetadot;
         if(zeta3[m] < 0) zeta3[m] = zeta3[m] + pi2;

         u3[m] = u2[m] + 0.5*dt*pzdot;
         //dum = 1.0
         //vxdum = (ezp/b+vpar/b*delbxp)*dum1
         //vzdum = (-exp1/b+vpar/b*delbzp)*dum1
         //vxdum = eyp+vpar/b*delbxp
         //w3(m)=w2(m) + 0.5*dt*(vxdum*kapxp + vzdum*kapzp+edot/ter)*dum*xnp

        if( x3[m] <= 2 * dxeq || x3[m] >= lx - 2 * dxeq ||
            z3[m] <= 2 * dzeq || z3[m] >= lz - 2 * dzeq) 
        {
          u3[m]=u2[m];
          x3[m]=x2[m];
          z3[m]=z2[m];
          zeta3[m]=zeta2[m];
          w3[m]=0;
        }
    }
    #pragma acc wait
    freeDeviceData();
    end_ppush_tm = MPI_Wtime();
    ppush_tm = ppush_tm + end_ppush_tm - start_ppush_tm;
}

//!-------------- End of subroutine ppush --------------------------------
void cpush_c_(const int &timestep){  //all warnings from this function are vars used in commented    //declared vars but not used in current version
        double exp1,ezp,ezetap,delbxp,delbzp,nudi0,ni_temp,T_center;                          //nudi=0,energy,rrr,energy0
        double wx0,wx1,wy0,wy1,wz0,wz1,dum1;                                                  //dum,vxdum,vzdum,vzetadum
        int m,i,j,k,l,k_plus_1=0;                                                             //p_m,
        double rhog,vfac,kapxp,kapzp,vpar,kaptxp,kapnxp,kaptzp,kapnzp,xnp;                    //pidum,bdcurlbp
        double b,enerb,ter,x,z,zeta;                                                          //th,r,qr
        double xt,zt,xdot,zdot,zetadot,pzdot,edot;                                            //xs,xdt,ydt,pzd0,vp0,vcurlbdotE
        double dbdxp,dbdzp,bfldp,bfldxp,bfldzp,bfldzetap, bstar, dbdzetap=0;
        double rhox[4], rhoy[4], curlbp[3], BStar3[3];

        //real(8),dimension(3)::curlbp,BStar3
        start_cpush_tm = MPI_Wtime();
        nudi0 = 1/sqrt(2.0)*18.4*pow(e,1.5)*(4.7140e-8)*(1.e-6);
        //write(*,*)t0i(200,201)
        T_center = t0i(imx/2,jmx/2);
        //write(*,*)T_center*e
        //write(*,*)nudi0
        
        prepareDeviceData();
        #pragma acc data \
        copyin(rand_table[0:10007]) \
        copy(mu[0:mmx], u2[0:mmx], u3[0:mmx], x3[0:mmx], z3[0:mmx], zeta3[0:mmx], w3[0:mmx], x2[0:mmx], z2[0:mmx],w2[0:mmx] ,zeta2[0:mmx])
        #pragma acc parallel loop gang vector private(BStar3,rhoy,rhox)
        for(m = 0; m < mm[0]; ++m){
            x=x3[m];
            i = static_cast<int>(x/dxeq);
            i = min(i,nx-1);
            wx0 = (i+1)-x/dxeq;
            wx1 = 1.-wx0;

            z = z3[m];
            k = static_cast<int>(z/dzeq);
            k = min(k,nz-1);
            wz0 = (k+1)-z/dzeq;
            wz1 = 1-wz0;

            ni_temp= wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) 
                 +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1); 
           //write(*,*) xn0i(i,k)


           // bdcurlbp =wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
                             //+wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1)
            #pragma loop seq
            for(j = 0; j <= 2; ++j){
                curlbp[j]= wx0*wz0*curlb(i,k,j)+wx0*wz1*curlb(i,k+1,j) 
                +wx1*wz0*curlb(i+1,k,j)+wx1*wz1*curlb(i+1,k+1,j);
            }

         dbdxp = wx0*wz0*dbdx(i,k)+wx0*wz1*dbdx(i,k+1) 
                 +wx1*wz0*dbdx(i+1,k)+wx1*wz1*dbdx(i+1,k+1); 
         dbdzp = wx0*wz0*dbdz(i,k)+wx0*wz1*dbdz(i,k+1) 
                 +wx1*wz0*dbdz(i+1,k)+wx1*wz1*dbdz(i+1,k+1);
         bfldp = wx0*wz0*b0(i,k)+wx0*wz1*b0(i,k+1) 
                 +wx1*wz0*b0(i+1,k)+wx1*wz1*b0(i+1,k+1); 
         bfldxp = wx0*wz0*b0x(i,k)+wx0*wz1*b0x(i,k+1) 
                 +wx1*wz0*b0x(i+1,k)+wx1*wz1*b0x(i+1,k+1); 
         bfldzp = wx0*wz0*b0z(i,k)+wx0*wz1*b0z(i,k+1) 
                 +wx1*wz0*b0z(i+1,k)+wx1*wz1*b0z(i+1,k+1); 
         bfldzetap = wx0*wz0*b0zeta(i,k)+wx0*wz1*b0zeta(i,k+1) 
                 +wx1*wz0*b0zeta(i+1,k)+wx1*wz1*b0zeta(i+1,k+1);
         ter = wx0*wz0*t0i(i,k)+wx0*wz1*t0i(i,k+1) 
                 +wx1*wz0*t0i(i+1,k)+wx1*wz1*t0i(i+1,k+1); 
         kaptxp = wx0*wz0*captix(i,k)+wx0*wz1*captix(i,k+1) 
                 +wx1*wz0*captix(i+1,k)+wx1*wz1*captix(i+1,k+1);
         kapnxp = wx0*wz0*capnix(i,k)+wx0*wz1*capnix(i,k+1)
                 +wx1*wz0*capnix(i+1,k)+wx1*wz1*capnix(i+1,k+1); 

         kaptzp = wx0*wz0*captiz(i,k)+wx0*wz1*captiz(i,k+1) 
                 +wx1*wz0*captiz(i+1,k)+wx1*wz1*captiz(i+1,k+1); 
         kapnzp = wx0*wz0*capniz(i,k)+wx0*wz1*capniz(i,k+1) 
                 +wx1*wz0*capniz(i+1,k)+wx1*wz1*capniz(i+1,k+1); 

         xnp = wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) 
                 +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1);

         b=1.-tor+tor*bfldp;

         rhog=sqrt(2.*b*mu[m]*mims[0])/(q[0]*b)*iflr;

         rhox[0] = rhog;
         rhoy[0] = 0.;
         rhox[1] = -rhox[0];
         rhoy[1] = -rhoy[0];
         rhox[2] = 0.;
         rhoy[2] = rhog;
         rhox[3] = 0.;
         rhoy[3] = -rhoy[2];
//    calculate avg. e-field...
//    do 1,2,4 point average, where lr is the no. of points...

         exp1=0.;
         ezp=0.;
         ezetap=0.;
         delbxp=0.;
         delbzp=0.;


        #pragma acc loop seq
        for(l = 0; l < lr[0]; ++l){

//SP            xs=x3(m)+rhox(l) !rwx(1,l)*rhog
            xt=x3[m]+rhox[l]; //rwx(1,l)*rhog
            zt=z3[m]+rhoy[l]; //(rwy(1,l)+sz*rwx(1,l))*rhog
//
//   particle can go out of bounds during gyroavg...
            if( (xt<2*dxeq)||(xt>lx-2*dxeq) ) xt=x3[m];
            if( (zt<2*dzeq)||(zt>lz-2*dzeq) ) zt=z3[m];
            zeta = my_fmod(zeta3[m], pi2);
            if(zeta < 0) zeta += pi2;
            i=static_cast<int>(xt/dx);
            j=static_cast<int>(zt/dz);
            k=static_cast<int>(zeta/dzeta);

            wx0=static_cast<float>((i+1)-xt/dx);
            wx1=1.-wx0;
            wy0=static_cast<float>((j+1)-zt/dz);
            wy1=1.-wy0;
            wz0=static_cast<float>((k+1)-zeta/dzeta);
            wz1=1.-wz0;

                k_plus_1=k+1;
              if(k==kmx) {k_plus_1=0;}
            exp1=exp1 + wx0*wy0*wz0*ex(i,j,k) + wx1*wy0*wz0*ex(i+1,j,k)
            + wx0*wy1*wz0*ex(i,j+1,k) + wx1*wy1*wz0*ex(i+1,j+1,k) + 
            wx0*wy0*wz1*ex(i,j, k_plus_1) + wx1*wy0*wz1*ex(i+1,j, k_plus_1) + 
            wx0*wy1*wz1*ex(i,j+1, k_plus_1) + wx1*wy1*wz1*ex(i+1,j+1, k_plus_1);

            ezp=ezp + wx0*wy0*wz0*ez(i,j,k) + wx1*wy0*wz0*ez(i+1,j,k) 
            + wx0*wy1*wz0*ez(i,j+1,k) + wx1*wy1*wz0*ez(i+1,j+1,k) + 
            wx0*wy0*wz1*ez(i,j, k_plus_1) + wx1*wy0*wz1*ez(i+1,j, k_plus_1) + 
            wx0*wy1*wz1*ez(i,j+1, k_plus_1) + wx1*wy1*wz1*ez(i+1,j+1, k_plus_1);

            ezetap =ezetap + wx0*wy0*wz0*ezeta(i,j,k) + wx1*wy0*wz0*ezeta(i+1,j,k) 
            + wx0*wy1*wz0*ezeta(i,j+1,k) + wx1*wy1*wz0*ezeta(i+1,j+1,k) + 
            wx0*wy0*wz1*ezeta(i,j, k_plus_1) + wx1*wy0*wz1*ezeta(i+1,j, k_plus_1) + 
            wx0*wy1*wz1*ezeta(i,j+1, k_plus_1) + wx1*wy1*wz1*ezeta(i+1,j+1, k_plus_1);

            delbxp =delbxp + wx0*wy0*wz0*delbx(i,j,k)  
            + wx1*wy0*wz0*delbx(i+1,j,k) 
            + wx0*wy1*wz0*delbx(i,j+1,k) 
            + wx1*wy1*wz0*delbx(i+1,j+1,k) 
            + wx0*wy0*wz1*delbx(i,j, k_plus_1) 
            + wx1*wy0*wz1*delbx(i+1,j, k_plus_1) 
            + wx0*wy1*wz1*delbx(i,j+1, k_plus_1) 
            + wx1*wy1*wz1*delbx(i+1,j+1, k_plus_1);

            delbzp =delbzp + wx0*wy0*wz0*delbz(i,j,k) 
            + wx1*wy0*wz0*delbz(i+1,j,k) 
            + wx0*wy1*wz0*delbz(i,j+1,k)  
            + wx1*wy1*wz0*delbz(i+1,j+1,k)  
            + wx0*wy0*wz1*delbz(i,j, k_plus_1)  
            + wx1*wy0*wz1*delbz(i+1,j, k_plus_1)  
            + wx0*wy1*wz1*delbz(i,j+1, k_plus_1)  
            + wx1*wy1*wz1*delbz(i+1,j+1, k_plus_1);
        }
         exp1 = exp1*0.25;
         ezp = ezp*0.25;
         ezetap = ezetap*0.25;
         delbxp = delbxp*0.25;
         delbzp = delbzp*0.25;
       
         vfac = 0.5*(mims[0]*pow(u2[m],2) + 2.*mu[m]*b);
         kapxp = kapnxp - (1.5-vfac/ter)*kaptxp;
         kapzp = kapnzp - (1.5-vfac/ter)*kaptzp;        
       
         vpar = u3[m];
         enerb=(mu[m]+mims[0]*vpar*vpar/b)/q[0]*tor;
       
         BStar3[0]=bfldxp+mims[0]*vpar*curlbp[0]/q[0]+delbxp;
         BStar3[1]=bfldzp+mims[0]*vpar*curlbp[1]/q[0]+delbzp; 
         BStar3[2]=bfldzetap+mims[0]*vpar*curlbp[2]/q[0];
                
       //bstar=b+mims(1)*vpar*bdcurlbp/q(1)
       
         bstar=(bfldxp*BStar3[0]+bfldzp*BStar3[1]+bfldzetap*BStar3[2])/bfldp;


        //vcurlbdotE=vpar*(exp1*curlbp(1)+ezp*curlbp(2)+ezetap*curlbp(3))


         dum1 = 1.;
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
         xdot = (vpar*BStar3[0]+(mu[m]*(bfldzp*dbdzetap-bfldzetap*dbdzp)/q[0]+(ezp*bfldzetap-ezetap*bfldzp))/(b))/bstar;
         zdot = (vpar*BStar3[1]+(mu[m]*(bfldzetap*dbdxp-bfldxp*dbdzetap)/q[0]+ (ezetap*bfldxp-exp1*bfldzetap))/(b))/bstar;
         zetadot = (vpar*BStar3[2]+(mu[m]*(bfldxp*dbdzp-bfldzp*dbdxp)/q[0]+(exp1*bfldzp-ezp*bfldxp))/(b))/bstar;

// !        pzd0 = -mu(m)/mims(1)/b*(bfldxp*dbdxp+bfldzp*dbdzp)
// !         pzdot = pzd0+(exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)/b*q(1)/mims(1)*nonlin
// !         pzdot = pzd0+(exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)/b*(q(1)/mims(1)+bdcurlbp*vpar/b)*nonlin
// !          pzdot = pzd0+((exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)*q(1)/mims(1)+vcurlbdotE)/bstar*nonlin


         pzdot = (BStar3[0]*(q[0]*exp1-mu[m]*dbdxp)+BStar3[1]*(q[0]*ezp-mu[m]*dbdzp)+BStar3[2]*(q[0]*ezetap-mu[m]*dbdzetap))/(mims[0]*bstar);
 

         edot = q[0]*(xdot*exp1+zdot*ezp+zetadot*ezetap);

         x3[m] = x2[m] + dt*xdot;
         z3[m] = z2[m] + dt*zdot;
         zeta3[m] = zeta2[m] + dt*zetadot;
         u3[m] = u2[m] + dt*pzdot;

// !         dum = 1.0
// !         vxdum = (ezp/b+vpar/b*delbxp)*dum1
// !         vzdum = (-exp1/b+vpar/b*delbzp)*dum1
// !         vxdum = eyp+vpar/b*delbxp
// !         w3(m)=w2(m) + dt*(vxdum*kapxp + vzdum*kapzp+edot/ter)*dum*xnp

         zeta3[m]= my_fmod(zeta3[m],pi2);
         if(zeta3[m] < 0) zeta3[m] = zeta3[m] + pi2;


//         write(*,*)energy, nudi
        
if( (x3[m]>2*dxeq) && (x3[m]<lx-2*dxeq) && (z3[m]>2*dzeq) && (z3[m]<lz-2*dzeq) ){
            //energy0 = (mu(m)*b+0.5*mims(1)*u3(m)**2)
            //energy =  max(energy0,0.1*T_center)
            //if (icollision==1) then
            //      nudi=nudi0*ni_temp/(energy)**1.5
            //           call random_number(rrr)
            //           p_m = int(2*rrr-1)
            //end if
            u2[m]=u3[m]; //*(1-nudi*dt)+rand_table(globle_integer)*sqrt((2*energy0/mims(1)-u3(m)**2)*nudi*dt)
            //  !$acc atomic
            //  globle_integer = mod(globle_integer+1,10007)
               
            //write(*,*) rand_table(globle_integer-1), u2(m),u3(m)
            x2[m]=x3[m];
            z2[m]=z3[m];
            zeta2[m]=zeta3[m];
            w2[m]=w3[m];
         
            // write(*,*)rand_table(globle_integer-1),mu(m),(energy0-0.5*mims(1)*u3(m)**2)/b
            // mu(m)= (energy0-0.5*mims(1)*u2(m)**2)/b
        } else {
            u3[m]=u2[m];
            x3[m]=x2[m];
            z3[m]=z2[m];
            zeta3[m]=zeta2[m];
            w2[m]=0.;
            w3[m]=0.;
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
  freeDeviceData();
  end_cpush_tm = MPI_Wtime();
  cpush_tm = cpush_tm + end_cpush_tm - start_cpush_tm;
}

inline void prepareDeviceData() {
    curlb.todev();
    ex.todev();  //potentially not needed since sequential loop
    ez.todev();
    ezeta.todev();
    dbdx.todev();
    dbdz.todev();
    b0.todev();
    b0x.todev();
    b0z.todev();
    b0zeta.todev();
    captix.todev();
    captiz.todev();
    capnix.todev();
    capniz.todev();
    xn0i.todev();
    delbx.todev();
    delbz.todev();
    t0i.todev();

//     // Grid-based fields
//     curlb.updatedev();
//     ex.updatedev();
//     ez.updatedev();
//     ezeta.updatedev();
//     dbdx.updatedev();
//     dbdz.updatedev();
//     b0.updatedev();
//     b0x.updatedev();
//     b0z.updatedev();
//     b0zeta.updatedev();
//     captix.updatedev();
//     captiz.updatedev();
//     capnix.updatedev();
//     capniz.updatedev();
//     xn0i.updatedev();
//     delbx.updatedev();
//     delbz.updatedev();
//     t0i.updatedev();
}

inline void freeDeviceData() {
    curlb.fromdev();
    ex.fromdev();
    ez.fromdev();
    ezeta.fromdev();
    dbdx.fromdev();
    dbdz.fromdev();
    b0.fromdev();
    b0x.fromdev();
    b0z.fromdev();
    b0zeta.fromdev();
    captix.fromdev();
    captiz.fromdev();
    capnix.fromdev();
    capniz.fromdev();
    xn0i.fromdev();
    delbx.fromdev();
    delbz.fromdev();
    t0i.fromdev();
}

//as a note: copyin will copy an array to gpu/cpu kernels.
//           Using copy will tell openacc to update the host arrays without this function. 
void updateHost1DArrays() {
    #pragma acc update self(x3[0:mmx])
    #pragma acc update self(z3[0:mmx])
    #pragma acc update self(zeta3[0:mmx])
    #pragma acc update self(u3[0:mmx])
    #pragma acc update self(w3[0:mmx])

    #pragma acc update self(x2[0:mmx])
    #pragma acc update self(z2[0:mmx])
    #pragma acc update self(zeta2[0:mmx])
    #pragma acc update self(u2[0:mmx])
    #pragma acc update self(w2[0:mmx])
}