!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       Ion pre-push
!
      subroutine ppush(n)

      use gemx_com
      use equil
      implicit none
      REAL(8) :: exp1,ezp,ezetap,delbxp,delbzp,energy, energy0,nudi0,nudi,T_center,ni_temp
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,dum,vxdum,vzdum,dum1,vzetadum
      INTEGER :: m,i,j,k,l,n,k_plus_1
      REAL(8) :: rhog,vfac,kapxp,kapzp,vpar,pidum,kaptxp,kapnxp,kaptzp,kapnzp,xnp,bdcurlbp
      REAL(8) :: b,th,r,enerb,qr,ter,x,z,zeta,bstar
      REAL(8) :: xt,zt,xdot,zdot,zetadot,xdt,ydt,pzdot,edot,pzd0,vp0
      REAL(8) :: dbdxp,dbdzp,bfldp,bfldxp,bfldzp,bfldzetap,vcurlbdotE,dbdzetap=0
      REAL(8) :: rhox(4),rhoy(4),psp,pzp,curlbp(3),Bstar3(3)
!      real(8),dimension(3)::curlbp,Bstar3
      start_ppush_tm = MPI_WTIME()


       nudi0 = 1/sqrt(2.0)*18.4*e**1.5*4.7140d-8*1.d-6
!      write(*,*)t0i(200,201)
       T_center = t0i(imx/2,jmx/2)

       
!$acc parallel loop gang vector private(rhoy,bstar3,rhox) copy(rand_table)
      do m=1,mm(1)
         x=x2(m)
         i = int(x/dxeq)
         i = min(i,nx-1)
         wx0 = (i+1)-x/dxeq
         wx1 = 1.-wx0

         z = z2(m)
         k = int(z/dzeq)
         k = min(k,nz-1)
         wz0 = (k+1)-z/dzeq
         wz1 = 1-wz0




 !        bdcurlbp =wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
 !                         +wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1)
         curlbp(:)= wx0*wz0*curlb(i,k,:)+wx0*wz1*curlb(i,k+1,:) &
              +wx1*wz0*curlb(i+1,k,:)+wx1*wz1*curlb(i+1,k+1,:)


         !         write(*,*) curlbp(1)-wx0*wz0*curlb(i,k,1)-wx0*wz1*curlb(i,k+1,1)-wx1*wz0*curlb(i+1,k,1)-wx1*wz1*curlb(i+1,k+1,1),curlbp(2)-wx0*wz0*curlb(i,k,2)-wx0*wz1*curlb(i,k+1,2)-wx1*wz0*curlb(i+1,k,2)-wx1*wz1*curlb(i+1,k+1,2),curlbp(3)-wx0*wz0*curlb(i,k,3)-wx0*wz1*curlb(i,k+1,3)-wx1*wz0*curlb(i+1,k,3)-wx1*wz1*curlb(i+1,k+1,3)
         !         write(*,*)curlbp(1),curlbp(2),curlbp(3)
!         write(*,*)bdcurlbp
         dbdxp = wx0*wz0*dbdx(i,k)+wx0*wz1*dbdx(i,k+1) &
                 +wx1*wz0*dbdx(i+1,k)+wx1*wz1*dbdx(i+1,k+1) 
         dbdzp = wx0*wz0*dbdz(i,k)+wx0*wz1*dbdz(i,k+1) &
                 +wx1*wz0*dbdz(i+1,k)+wx1*wz1*dbdz(i+1,k+1) 
         bfldp = wx0*wz0*b0(i,k)+wx0*wz1*b0(i,k+1) &
                 +wx1*wz0*b0(i+1,k)+wx1*wz1*b0(i+1,k+1) 
         bfldxp = wx0*wz0*b0x(i,k)+wx0*wz1*b0x(i,k+1) &
                 +wx1*wz0*b0x(i+1,k)+wx1*wz1*b0x(i+1,k+1) 
         bfldzp = wx0*wz0*b0z(i,k)+wx0*wz1*b0z(i,k+1) &
                 +wx1*wz0*b0z(i+1,k)+wx1*wz1*b0z(i+1,k+1) 
         bfldzetap = wx0*wz0*b0zeta(i,k)+wx0*wz1*b0zeta(i,k+1) &
                 +wx1*wz0*b0zeta(i+1,k)+wx1*wz1*b0zeta(i+1,k+1) 
         ter = wx0*wz0*t0i(i,k)+wx0*wz1*t0i(i,k+1) &
                 +wx1*wz0*t0i(i+1,k)+wx1*wz1*t0i(i+1,k+1) 
         kaptxp = wx0*wz0*captix(i,k)+wx0*wz1*captix(i,k+1) &
                 +wx1*wz0*captix(i+1,k)+wx1*wz1*captix(i+1,k+1) 
         kapnxp = wx0*wz0*capnix(i,k)+wx0*wz1*capnix(i,k+1) &
                 +wx1*wz0*capnix(i+1,k)+wx1*wz1*capnix(i+1,k+1) 

         kaptzp = wx0*wz0*captiz(i,k)+wx0*wz1*captiz(i,k+1) &
                 +wx1*wz0*captiz(i+1,k)+wx1*wz1*captiz(i+1,k+1) 
         kapnzp = wx0*wz0*capniz(i,k)+wx0*wz1*capniz(i,k+1) &
                 +wx1*wz0*capniz(i+1,k)+wx1*wz1*capniz(i+1,k+1) 

         xnp = wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) &
                 +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1) 

         b=1.-tor+tor*bfldp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pitch angle collision!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ni_temp= wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) &
                 +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1) 
         energy0 = (mu(m)*b+0.5*mims(1)*u3(m)**2)
         energy =  max(energy0,0.1*T_center)
         if (icollision==1) then
                  nudi=nudi0*ni_temp/(energy)**1.5
                        ! call random_number(rrr)
                        ! p_m = int(2*rrr-1)
         end if

!         write(*,*) rand_table(globle_integer-1), u2(m),u2(m)*(1-nudi*dt)+rand_table(globle_integer)*sqrt((2*energy0/mims(1)-u2(m)**2)*nudi*dt)!,mu(m),(energy0-0.5*mims(1)*u2(m)**2)/b

         u2(m)=u2(m)*(1-nudi*dt)+rand_table(globle_integer)*sqrt((2*energy0/mims(1)-u2(m)**2)*nudi*dt)
         globle_integer = mod(globle_integer+1,10007)
!         write(*,*) rand_table(globle_integer-1), u2(m),u3(m)!,mu(m),(energy0-0.5*mims(1)*u2(m)**2)/b
!         write(*,*) rand_table(globle_integer-1), mu(m),(energy0-0.5*mims(1)*u2(m)**2)/b

         mu(m)= (energy0-0.5*mims(1)*u2(m)**2)/b

         


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end of pitch angle collision!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
        



          
         rhog=sqrt(2.*b*mu(m)*mims(1))/(q(1)*b)*iflr

         rhox(1) = rhog
         rhoy(1) = 0.
         rhox(2) = -rhox(1)
         rhoy(2) = -rhoy(1)
         rhox(3) = 0
         rhoy(3) = rhog
         rhox(4) = 0
         rhoy(4) = -rhoy(3)
!    calculate avg. e-field...
!    do 1,2,4 point average, where lr is the no. of points...

         exp1=0.
         ezp=0.
         ezetap=0.
         delbxp=0.
         delbzp=0.

!  4 pt. avg. done explicitly for vectorization...
         !$acc loop seq
         do 200 l=1,lr(1)
!
            xt=x2(m)+rhox(l) !rwx(1,l)*rhog
            zt=z2(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
!            zeta=modulo(zeta2(m),pi2)
!     
!   particle can go out of bounds during gyroavg...
            if( (xt<2*dxeq).or.(xt>lx-2*dxeq) ) xt=x2(m)
            if( (zt<2*dzeq).or.(zt>lz-2*dzeq) ) zt=z2(m)
            zeta=zeta2(m)
!            xt=modulo(xs,xdim)
!            zt=modulo(zt,zdim)
                       i=int(xt/dx)
            j=int(zt/dz)
            k=int(zeta/dzeta)


            wx0=float(i+1)-xt/dx
            wx1=1.-wx0
            wy0=float(j+1)-zt/dz
            wy1=1.-wy0
            wz0=float(k+1)-zeta/dzeta
            wz1=1.-wz0

              k_plus_1=k+1
              if(k==kmx) k_plus_1=0
            exp1=exp1 + wx0*wy0*wz0*ex(i,j,k) + wx1*wy0*wz0*ex(i+1,j,k) &
            + wx0*wy1*wz0*ex(i,j+1,k) + wx1*wy1*wz0*ex(i+1,j+1,k) + &
            wx0*wy0*wz1*ex(i,j, k_plus_1) + wx1*wy0*wz1*ex(i+1,j, k_plus_1) + &
            wx0*wy1*wz1*ex(i,j+1, k_plus_1) + wx1*wy1*wz1*ex(i+1,j+1, k_plus_1)

            ezp=ezp + wx0*wy0*wz0*ez(i,j,k) + wx1*wy0*wz0*ez(i+1,j,k) &
            + wx0*wy1*wz0*ez(i,j+1,k) + wx1*wy1*wz0*ez(i+1,j+1,k) + &
            wx0*wy0*wz1*ez(i,j, k_plus_1) + wx1*wy0*wz1*ez(i+1,j, k_plus_1) + &
            wx0*wy1*wz1*ez(i,j+1, k_plus_1) + wx1*wy1*wz1*ez(i+1,j+1, k_plus_1)

            ezetap =ezetap + wx0*wy0*wz0*ezeta(i,j,k) + wx1*wy0*wz0*ezeta(i+1,j,k) &
            + wx0*wy1*wz0*ezeta(i,j+1,k) + wx1*wy1*wz0*ezeta(i+1,j+1,k) + &
            wx0*wy0*wz1*ezeta(i,j, k_plus_1) + wx1*wy0*wz1*ezeta(i+1,j, k_plus_1) + &
            wx0*wy1*wz1*ezeta(i,j+1, k_plus_1) + wx1*wy1*wz1*ezeta(i+1,j+1, k_plus_1)

            delbxp =delbxp + wx0*wy0*wz0*delbx(i,j,k)  &
            + wx1*wy0*wz0*delbx(i+1,j,k) &
            + wx0*wy1*wz0*delbx(i,j+1,k) &
            + wx1*wy1*wz0*delbx(i+1,j+1,k) &
            + wx0*wy0*wz1*delbx(i,j, k_plus_1) &
            + wx1*wy0*wz1*delbx(i+1,j, k_plus_1) &
            + wx0*wy1*wz1*delbx(i,j+1, k_plus_1) &
            + wx1*wy1*wz1*delbx(i+1,j+1, k_plus_1)

            delbzp =delbzp + wx0*wy0*wz0*delbz(i,j,k) &
            + wx1*wy0*wz0*delbz(i+1,j,k) &
            + wx0*wy1*wz0*delbz(i,j+1,k)  &
            + wx1*wy1*wz0*delbz(i+1,j+1,k)  &
            + wx0*wy0*wz1*delbz(i,j, k_plus_1)  &
            + wx1*wy0*wz1*delbz(i+1,j, k_plus_1)  &
            + wx0*wy1*wz1*delbz(i,j+1, k_plus_1)  &
            + wx1*wy1*wz1*delbz(i+1,j+1, k_plus_1)
 200     continue
         exp1 = exp1/4.
         ezp = ezp/4.
         ezetap = ezetap/4.
         delbxp = delbxp/4.
         delbzp = delbzp/4.
!
         vfac = 0.5*(mims(1)*u2(m)**2 + 2.*mu(m)*b)
         kapxp = kapnxp - (1.5-vfac/ter)*kaptxp
         kapzp = kapnzp - (1.5-vfac/ter)*kaptzp         

         vpar = u2(m)
         enerb=(mu(m)+mims(1)*vpar*vpar/b)/q(1)*tor

         Bstar3(1)=bfldxp +mims(1)*vpar*curlbp(1)/q(1)+delbxp
         Bstar3(2)=bfldzp+mims(1)*vpar*curlbp(2)/q(1)+delbzp
         Bstar3(3)=bfldzetap+mims(1)*vpar*curlbp(3)/q(1)


!         bstar=b+mims(1)*vpar*bdcurlbp/q(1)
         bstar=(bfldxp*Bstar3(1)+bfldzp*Bstar3(2)+bfldzetap*Bstar3(3))/bfldp
!         write(*,*)bstar-b-mims(1)*vpar*bdcurlbp/q(1), b-bfldp
!         vcurlbdotE=vpar*(exp1*curlbp(1)+ezp*curlbp(2)+ezetap*curlbp(3))
         
         dum1 = 1.
         ! vxdum = (ezp/b+vpar/b*delbxp)*dum1
!         vxdum = (ezp*bfldzetap-ezetap*bfldzp)/b**2
!         xdot = vxdum*nonlin +vpar*bfldxp/b-enerb/bfldp/bfldp*bfldzetap*dbdzp
!         vzdum = (ezetap*bfldxp-exp1*bfldzetap)/b**2
 !        zdot = (-exp1/b+vpar/b*delbzp)*dum1*nonlin &
 !            +vpar*bfldzp/b+enerb/bfldp/bfldp*bfldzetap*dbdxp
!         zdot = vzdum*nonlin+vpar*bfldzp/b+enerb/bfldp/bfldp*bfldzetap*dbdxp

!         vzetadum= (exp1*bfldzp-ezp*bfldxp)/b**2
!         zetadot = vzetadum/x*nonlin + vpar*bfldzetap/(x*b)+enerb/(x*b*b)*(bfldxp*dbdzp-bfldzp*dbdxp)



!         write(*,*) dbdzetap
         xdot = (vpar*Bstar3(1)+(mu(m)*(bfldzp*dbdzetap-bfldzetap*dbdzp)/q(1)+(ezp*bfldzetap-ezetap*bfldzp))/(b))/bstar
         zdot = (vpar*Bstar3(2)+(mu(m)*(bfldzetap*dbdxp-bfldxp*dbdzetap)/q(1)+ (ezetap*bfldxp-exp1*bfldzetap))/(b))/bstar
         zetadot = (vpar*Bstar3(3)+(mu(m)*(bfldxp*dbdzp-bfldzp*dbdxp)/q(1)+(exp1*bfldzp-ezp*bfldxp))/(b))/bstar


         
!         pzd0 = -mu(m)/mims(1)/b*(bfldxp*dbdxp+bfldzp*dbdzp)
!         write(*,*)
!         write(*,*)(exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)/b*q(1)/mims(1)
!         pzdot = pzd0+(exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)/b*(q(1)/mims(1)+bdcurlbp*vpar/b)*nonlin

!          pzdot = pzd0+((exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)*q(1)/mims(1)+vcurlbdotE)/bstar*nonlin


         pzdot = (Bstar3(1)*(q(1)*exp1-mu(m)*dbdxp)+Bstar3(2)*(q(1)*ezp-mu(m)*dbdzp)+Bstar3(3)*(q(1)*ezetap-mu(m)*dbdzetap))/(mims(1)*bstar)
          
         
         edot = q(1)*(xdot*exp1+zdot*ezp+zetadot*ezetap)

         x3(m) = x2(m) + 0.5*dt*xdot
         z3(m) = z2(m) + 0.5*dt*zdot
         zeta3(m) = zeta2(m) + 0.5*dt*zetadot
         u3(m) = u2(m) + 0.5*dt*pzdot

!        dum = 1.0
!        vxdum = (ezp/b+vpar/b*delbxp)*dum1
!        vzdum = (-exp1/b+vpar/b*delbzp)*dum1
!         vxdum = eyp+vpar/b*delbxp
!         w3(m)=w2(m) + 0.5*dt*(vxdum*kapxp + vzdum*kapzp+edot/ter)*dum*xnp
         
      if( (x3(m)>2*dxeq).and.(x3(m)<lx-2*dxeq).and.(z3(m)>2*dzeq).and.(z3(m)<lz-2*dzeq) ) then
      else
          u3(m)=u2(m)
          x3(m)=x2(m)
          z3(m)=z2(m)
          zeta3(m)=zeta2(m)
          w3(m)=0.
      endif

   enddo
!!   !$acc wait

      end_ppush_tm = MPI_WTIME()
      ppush_tm = ppush_tm + end_ppush_tm - start_ppush_tm
      return
      end
!
!-------------- End of subroutine ppush --------------------------------


      subroutine cpush(n)
!-----------------------------------------------------------------------
!              Ion corrector push
!-----------------------------------------------------------------------
      use gemx_com
      use equil
      implicit none
      REAL(8) :: exp1,ezp,ezetap,delbxp,delbzp,nudi0,nudi=0,ni_temp,energy,rrr,T_center,energy0
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,dum,vxdum,vzdum,dum1,vzetadum
      INTEGER :: m,i,j,k,l,n,k_plus_1=0,p_m
      REAL(8) :: rhog,vfac,kapxp,kapzp,vpar,pidum,kaptxp,kapnxp,kaptzp,kapnzp,xnp,bdcurlbp
      REAL(8) :: b,th,r,enerb,qr,ter,x,z,zeta
      REAL(8) :: xt,xs,zt,xdot,zdot,zetadot,xdt,ydt,pzdot,edot,pzd0,vp0,vcurlbdotE
      REAL(8) :: dbdxp,dbdzp,bfldp,bfldxp,bfldzp,bfldzetap, bstar, dbdzetap=0
      REAL(8) :: rhox(4),rhoy(4),psp,pzp,curlbp(3),Bstar3(3)
!      real(8),dimension(3)::curlbp,Bstar3
      start_cpush_tm = MPI_WTIME()
      nudi0 = 1/sqrt(2.0)*18.4*e**1.5*4.7140d-8*1.d-6
!      write(*,*)t0i(200,201)
      T_center = t0i(imx/2,jmx/2)
!      write(*,*)T_center*e
      !      write(*,*)nudi0

      
!$acc parallel loop gang vector private(bstar3,rhoy,rhox) copy(rand_table)
      do m=1,mm(1)
         x=x3(m)
         i = int(x/dxeq)
         i = min(i,nx-1)
         wx0 = (i+1)-x/dxeq
         wx1 = 1.-wx0

         z = z3(m)
         k = int(z/dzeq)
         k = min(k,nz-1)
         wz0 = (k+1)-z/dzeq
         wz1 = 1-wz0

         ni_temp= wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) &
                 +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1) 

 !        write(*,*) xn0i(i,k)

         
!         bdcurlbp =wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
!                          +wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1)
         curlbp(:)= wx0*wz0*curlb(i,k,:)+wx0*wz1*curlb(i,k+1,:) &
                 +wx1*wz0*curlb(i+1,k,:)+wx1*wz1*curlb(i+1,k+1,:)

         dbdxp = wx0*wz0*dbdx(i,k)+wx0*wz1*dbdx(i,k+1) &
                 +wx1*wz0*dbdx(i+1,k)+wx1*wz1*dbdx(i+1,k+1) 
         dbdzp = wx0*wz0*dbdz(i,k)+wx0*wz1*dbdz(i,k+1) &
                 +wx1*wz0*dbdz(i+1,k)+wx1*wz1*dbdz(i+1,k+1) 
         bfldp = wx0*wz0*b0(i,k)+wx0*wz1*b0(i,k+1) &
                 +wx1*wz0*b0(i+1,k)+wx1*wz1*b0(i+1,k+1) 
         bfldxp = wx0*wz0*b0x(i,k)+wx0*wz1*b0x(i,k+1) &
                 +wx1*wz0*b0x(i+1,k)+wx1*wz1*b0x(i+1,k+1) 
         bfldzp = wx0*wz0*b0z(i,k)+wx0*wz1*b0z(i,k+1) &
                 +wx1*wz0*b0z(i+1,k)+wx1*wz1*b0z(i+1,k+1) 
         bfldzetap = wx0*wz0*b0zeta(i,k)+wx0*wz1*b0zeta(i,k+1) &
                 +wx1*wz0*b0zeta(i+1,k)+wx1*wz1*b0zeta(i+1,k+1) 
         ter = wx0*wz0*t0i(i,k)+wx0*wz1*t0i(i,k+1) &
                 +wx1*wz0*t0i(i+1,k)+wx1*wz1*t0i(i+1,k+1) 
         kaptxp = wx0*wz0*captix(i,k)+wx0*wz1*captix(i,k+1) &
                 +wx1*wz0*captix(i+1,k)+wx1*wz1*captix(i+1,k+1) 
         kapnxp = wx0*wz0*capnix(i,k)+wx0*wz1*capnix(i,k+1) &
                 +wx1*wz0*capnix(i+1,k)+wx1*wz1*capnix(i+1,k+1) 

         kaptzp = wx0*wz0*captiz(i,k)+wx0*wz1*captiz(i,k+1) &
                 +wx1*wz0*captiz(i+1,k)+wx1*wz1*captiz(i+1,k+1) 
         kapnzp = wx0*wz0*capniz(i,k)+wx0*wz1*capniz(i,k+1) &
                 +wx1*wz0*capniz(i+1,k)+wx1*wz1*capniz(i+1,k+1) 

         xnp = wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) &
                 +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1) 

         b=1.-tor+tor*bfldp

         rhog=sqrt(2.*b*mu(m)*mims(1))/(q(1)*b)*iflr

         rhox(1) = rhog
         rhoy(1) = 0.
         rhox(2) = -rhox(1)
         rhoy(2) = -rhoy(1)
         rhox(3) = 0
         rhoy(3) = rhog
         rhox(4) = 0
         rhoy(4) = -rhoy(3)
!    calculate avg. e-field...
!    do 1,2,4 point average, where lr is the no. of points...

         exp1=0.
         ezp=0.
         ezetap=0.
         delbxp=0.
         delbzp=0.

!  4 pt. avg. done explicitly for vectorization...
         !$acc loop seq
         do 200 l=1,lr(1)
!
!SP            xs=x3(m)+rhox(l) !rwx(1,l)*rhog
            xt=x3(m)+rhox(l) !rwx(1,l)*rhog
            zt=z3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
!
!   particle can go out of bounds during gyroavg...
            if( (xt<2*dxeq).or.(xt>lx-2*dxeq) ) xt=x3(m)
            if( (zt<2*dzeq).or.(zt>lz-2*dzeq) ) zt=z3(m)
            zeta=modulo(zeta3(m),2*pi)
           i=int(xt/dx)
            j=int(zt/dz)
            k=int(zeta/dzeta)


            wx0=float(i+1)-xt/dx
            wx1=1.-wx0
            wy0=float(j+1)-zt/dz
            wy1=1.-wy0
            wz0=float(k+1)-zeta/dzeta
            wz1=1.-wz0

              k_plus_1=k+1
              if(k==kmx) k_plus_1=0
            exp1=exp1 + wx0*wy0*wz0*ex(i,j,k) + wx1*wy0*wz0*ex(i+1,j,k) &
            + wx0*wy1*wz0*ex(i,j+1,k) + wx1*wy1*wz0*ex(i+1,j+1,k) + &
            wx0*wy0*wz1*ex(i,j, k_plus_1) + wx1*wy0*wz1*ex(i+1,j, k_plus_1) + &
            wx0*wy1*wz1*ex(i,j+1, k_plus_1) + wx1*wy1*wz1*ex(i+1,j+1, k_plus_1)

            ezp=ezp + wx0*wy0*wz0*ez(i,j,k) + wx1*wy0*wz0*ez(i+1,j,k) &
            + wx0*wy1*wz0*ez(i,j+1,k) + wx1*wy1*wz0*ez(i+1,j+1,k) + &
            wx0*wy0*wz1*ez(i,j, k_plus_1) + wx1*wy0*wz1*ez(i+1,j, k_plus_1) + &
            wx0*wy1*wz1*ez(i,j+1, k_plus_1) + wx1*wy1*wz1*ez(i+1,j+1, k_plus_1)

            ezetap =ezetap + wx0*wy0*wz0*ezeta(i,j,k) + wx1*wy0*wz0*ezeta(i+1,j,k) &
            + wx0*wy1*wz0*ezeta(i,j+1,k) + wx1*wy1*wz0*ezeta(i+1,j+1,k) + &
            wx0*wy0*wz1*ezeta(i,j, k_plus_1) + wx1*wy0*wz1*ezeta(i+1,j, k_plus_1) + &
            wx0*wy1*wz1*ezeta(i,j+1, k_plus_1) + wx1*wy1*wz1*ezeta(i+1,j+1, k_plus_1)

            delbxp =delbxp + wx0*wy0*wz0*delbx(i,j,k)  &
            + wx1*wy0*wz0*delbx(i+1,j,k) &
            + wx0*wy1*wz0*delbx(i,j+1,k) &
            + wx1*wy1*wz0*delbx(i+1,j+1,k) &
            + wx0*wy0*wz1*delbx(i,j, k_plus_1) &
            + wx1*wy0*wz1*delbx(i+1,j, k_plus_1) &
            + wx0*wy1*wz1*delbx(i,j+1, k_plus_1) &
            + wx1*wy1*wz1*delbx(i+1,j+1, k_plus_1)

            delbzp =delbzp + wx0*wy0*wz0*delbz(i,j,k) &
            + wx1*wy0*wz0*delbz(i+1,j,k) &
            + wx0*wy1*wz0*delbz(i,j+1,k)  &
            + wx1*wy1*wz0*delbz(i+1,j+1,k)  &
            + wx0*wy0*wz1*delbz(i,j, k_plus_1)  &
            + wx1*wy0*wz1*delbz(i+1,j, k_plus_1)  &
            + wx0*wy1*wz1*delbz(i,j+1, k_plus_1)  &
            + wx1*wy1*wz1*delbz(i+1,j+1, k_plus_1)
 200     continue
         exp1 = exp1/4.
         ezp = ezp/4.
         ezetap = ezetap/4.
         delbxp = delbxp/4.
         delbzp = delbzp/4.
!
         vfac = 0.5*(mims(1)*u2(m)**2 + 2.*mu(m)*b)
         kapxp = kapnxp - (1.5-vfac/ter)*kaptxp
         kapzp = kapnzp - (1.5-vfac/ter)*kaptzp         

         vpar = u3(m)
         enerb=(mu(m)+mims(1)*vpar*vpar/b)/q(1)*tor

         Bstar3(1)=bfldxp+mims(1)*vpar*curlbp(1)/q(1)+delbxp
         Bstar3(2)=bfldzp+mims(1)*vpar*curlbp(2)/q(1)+delbzp
         Bstar3(3)=bfldzetap+mims(1)*vpar*curlbp(3)/q(1)

         
!        bstar=b+mims(1)*vpar*bdcurlbp/q(1)

         bstar=(bfldxp*Bstar3(1)+bfldzp*Bstar3(2)+bfldzetap*Bstar3(3))/bfldp
         

!         vcurlbdotE=vpar*(exp1*curlbp(1)+ezp*curlbp(2)+ezetap*curlbp(3))


         dum1 = 1.
         !         vxdum = (ezp/b+vpar/b*delbxp)*dum1
!         vxdum =(ezp*bfldzetap-ezetap*bfldzp)/b**2
!         xdot = vxdum*nonlin +vpar*bfldxp/b-enerb/bfldp/bfldp*bfldzetap*dbdzp
         !         write(*,*) vxdum*nonlin, (-exp1/b+vpar/b*delbzp)*dum1*nonlin
!         vzdum =(ezetap*bfldxp-exp1*bfldzetap)/b**2
 !        zdot = (-exp1/b+vpar/b*delbzp)*dum1*nonlin &
         !            +vpar*bfldzp/b+enerb/bfldp/bfldp*bfldzetap*dbdxp
!         write(*,*)vzdum,vxdum
!         zdot = vzdum*nonlin+vpar*bfldzp/b+enerb/bfldp/bfldp*bfldzetap*dbdxp

!         zetadot =  vpar*bfldzetap/(x*b)+enerb/(x*b*b)*(bfldxp*dbdzp-bfldzp*dbdxp)

!         vzetadum= (exp1*bfldzp-ezp*bfldxp)/b**2
!         zetadot = vzetadum/x*nonlin + vpar*bfldzetap/(x*b)+enerb/(x*b*b)*(bfldxp*dbdzp-bfldzp*dbdxp)


!         write(*,*)dbdzetap
         xdot = (vpar*Bstar3(1)+(mu(m)*(bfldzp*dbdzetap-bfldzetap*dbdzp)/q(1)+(ezp*bfldzetap-ezetap*bfldzp))/(b))/bstar
         zdot = (vpar*Bstar3(2)+(mu(m)*(bfldzetap*dbdxp-bfldxp*dbdzetap)/q(1)+ (ezetap*bfldxp-exp1*bfldzetap))/(b))/bstar
         zetadot = (vpar*Bstar3(3)+(mu(m)*(bfldxp*dbdzp-bfldzp*dbdxp)/q(1)+(exp1*bfldzp-ezp*bfldxp))/(b))/bstar


         
 !        pzd0 = -mu(m)/mims(1)/b*(bfldxp*dbdxp+bfldzp*dbdzp)
         !         pzdot = pzd0+(exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)/b*q(1)/mims(1)*nonlin
!         pzdot = pzd0+(exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)/b*(q(1)/mims(1)+bdcurlbp*vpar/b)*nonlin
!          pzdot = pzd0+((exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)*q(1)/mims(1)+vcurlbdotE)/bstar*nonlin


         pzdot = (Bstar3(1)*(q(1)*exp1-mu(m)*dbdxp)+Bstar3(2)*(q(1)*ezp-mu(m)*dbdzp)+Bstar3(3)*(q(1)*ezetap-mu(m)*dbdzetap))/(mims(1)*bstar)
          

         edot = q(1)*(xdot*exp1+zdot*ezp+zetadot*ezetap)

         x3(m) = x2(m) + dt*xdot
         z3(m) = z2(m) + dt*zdot
         zeta3(m) = zeta2(m) + dt*zetadot
         u3(m) = u2(m) + dt*pzdot

!         dum = 1.0
!         vxdum = (ezp/b+vpar/b*delbxp)*dum1
!         vzdum = (-exp1/b+vpar/b*delbzp)*dum1
!         vxdum = eyp+vpar/b*delbxp
!         w3(m)=w2(m) + dt*(vxdum*kapxp + vzdum*kapzp+edot/ter)*dum*xnp
         
         zeta3(m)=modulo(zeta3(m),pi2)

         
!         write(*,*)energy, nudi
         if( (x3(m)>2*dxeq).and.(x3(m)<lx-2*dxeq).and.(z3(m)>2*dzeq).and.(z3(m)<lz-2*dzeq) ) then
!            energy0 = (mu(m)*b+0.5*mims(1)*u3(m)**2)
!            energy =  max(energy0,0.1*T_center)
!            if (icollision==1) then
!                  nudi=nudi0*ni_temp/(energy)**1.5
                        ! call random_number(rrr)
                        ! p_m = int(2*rrr-1)
!            end if
            u2(m)=u3(m)!*(1-nudi*dt)+rand_table(globle_integer)*sqrt((2*energy0/mims(1)-u3(m)**2)*nudi*dt)
!!            !$acc atomic
!            globle_integer = mod(globle_integer+1,10007)
               
!!          write(*,*) rand_table(globle_integer-1), u2(m),u3(m)
          x2(m)=x3(m)
          z2(m)=z3(m)
          zeta2(m)=zeta3(m)
          w2(m)=w3(m)
          !write(*,*)rand_table(globle_integer-1),mu(m),(energy0-0.5*mims(1)*u3(m)**2)/b
         ! mu(m)= (energy0-0.5*mims(1)*u2(m)**2)/b
          else
          u3(m)=u2(m)
          x3(m)=x2(m)
          z3(m)=z2(m)
          zeta3(m)=zeta2(m)
          w2(m)=0.
          w3(m)=0.
      endif



      !if(Myid==0)then
!  open(935, file='flag_debug',status='unknown',position='append')
!  write(935,*)'before PETSc sovling'
!  close(935)
!end if
      !if(myid==0 .and. m==1)then
      !     open(93, file='test_energy',status='unknown',position='append')
      !     write(93,*) mu(m)*b+0.5*mims(1)*u2(m)**2-q(1)*z2(m)*ez(i,j,k)
      !     close(93)
        !       write(*,*) mu(m)*b+0.5*mims(1)*u2(m)**2-q(1)*z2(m)*ez(i,j,k)
        ! close(19)
      !end if
      


      

   enddo

!!  !$acc wait

      end_cpush_tm = MPI_WTIME()
      cpush_tm = cpush_tm + end_cpush_tm - start_cpush_tm
      return
      end
