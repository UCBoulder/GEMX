!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       Ion pre-push
!
      subroutine ppush(n)

      use gemx_com
      use equil
      implicit none
      REAL(8) :: exp1,ezp,ezetap,delbxp,delbzp,energy, energy0,nudi0,nudi,T_center,ni_temp,phibp,Rbp,Zbp,zetalbp,ebRp,ebZp,ebzetap
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,dum,vxdum,vzdum,dum1,vzetadum,psiparticle,phi_particle,wkk,dldzetabp,ebparp,dzetab
      INTEGER :: m,i,j,k,l,n,k_plus_1,kk,ll_plus_1,ll,nb,nb_temp=10
      REAL(8) :: rhog,vfac,kapxp,kapzp,vpar,pidum,kaptxp,kapnxp,kaptzp,kapnzp,xnp,bdcurlbp
      REAL(8) :: b,th,r,enerb,qr,ter,x,z,zeta,bstar,drtemp,dztemp
      REAL(8) :: xt,zt,xdot,zdot,zetadot,xdt,ydt,pzdot,edot,pzd0,vp0
      REAL(8) :: dbdxp,dbdzp,bfldp,bfldxp,bfldzp,bfldzetap,vcurlbdotE,dbdzetap=0
      REAL(8) :: rhox(4),rhoy(4),psp,pzp,curlbp(3),Bstar3(3),tempb(60)
!      real(8),dimension(3)::curlbp,Bstar3
      start_ppush_tm = MPI_WTIME()


       nudi0 = 1/sqrt(2.0)*18.4*e**1.5*4.7140d-8*1.d-6
!      write(*,*)t0i(200,201)
       T_center = t0i(imx/2,jmx/2)
       dzetab=2*pi/nzeta_blob!zetab_table(2)-zetab_table(1)

       
!$acc parallel loop gang vector private(rhoy,Bstar3,rhox,curlbp,tempb) copy(rand_table)
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
         psiparticle = wx0*wz0*psi_p(i,k)+wx0*wz1*psi_p(i,k+1) &
              +wx1*wz0*psi_p(i+1,k)+wx1*wz1*psi_p(i+1,k+1)
         phi_particle = wx0*wz0*phi(i,k,0)+wx0*wz1*phi(i,k+1,0) &
              +wx1*wz0*phi(i+1,k,0)+wx1*wz1*phi(i+1,k+1,0)

         b=1.-tor+tor*bfldp



!         write(*,*)mims(1),u2(m),mu(m),q(1)
         if (w2(m) > 1.e-100) then
            if( i<xpoint_i_10 .and. k<xpoint_j_10+10)then
               
               w2(m)=0
            else if (a_div*(x2(m)+xctr-xdim/2)+b_div-z2(m)>0 .and. k<xpoint_j_10+10) then
               w2(m)=0
!               write(*,*) psiparticle, u2(m), (mu(m)*b+0.5*mims(1)*u2(m)**2)
               hit_div(m,:)= (/psiparticle,u2(m), mu(m)*b+0.5*mims(1)*u2(m)**2/)
               hit_t(m) = timestep
!               hit_div(m,2)= u2(m)
!               hit_div(m,3)= mu(m)*b+0.5*mims(1)*u2(m)**2
            end if
         end if
         



!         write(*,*) (mu(m)*b+0.5*mims(1)*u3(m)**2+q(1)*phi_particle)
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pitch angle collision!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         ni_temp= wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) &
!                 +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1) 
!         energy0 = (mu(m)*b+0.5*mims(1)*u3(m)**2)
!         energy =  max(energy0,0.1*T_center)
!         if (icollision==1.and.timestep>nstep_E) then
!                  nudi=nudi0*ni_temp/(energy)**1.5

         


!         u2(m)=u2(m)*(1-nudi*dt)+rand_table(globle_integer)*sqrt((2*energy0/mims(1)-u2(m)**2)*nudi*dt)
!         globle_integer = mod(globle_integer+1,10007)


!         mu(m)= (energy0-0.5*mims(1)*u2(m)**2)/b
         !end if
    

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end of pitch angle collision!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
        
   

          
         rhog=sqrt(2.*b*mu(m)*mims(1))/(q(1)*b)*iflr

!         write(*,*)rhog

         rhox(1) = rhog
         rhoy(1) = 0.
         rhox(2) = -rhog!-rhox(1)
         rhoy(2) = 0!-rhoy(1)
         rhox(3) = 0
         rhoy(3) = rhog
         rhox(4) = 0
         rhoy(4) = -rhog!-rhoy(3)
!    calculate avg. e-field...
!    do 1,2,4 point average, where lr is the no. of points...

         exp1=0.
         ezp=0.
         ezetap=0.
         delbxp=0.
         delbzp=0.

!         ebparp=0
!         ebRp=0
         ebZp=0
         !         ebzetap=0
         Rbp=10000
         Zbp=10000

           zeta=zeta2(m)

           do nb=1,nblob
           wkk=(modulo(zeta-zeta0_blob(nb),2*pi)/dzetab)
!           write(*,*)zeta0_blob(nb),zeta,wkk
        
           kk=int(wkk)+1
           wkk=kk-wkk
           do ll=1,ntern_blob
              
!              write(*,*)kk
              dldzetabp=dldzetab_table(kk)*wkk+dldzetab_table(kk+1)*(1-wkk)              
              drtemp=Rb_table(nb,kk)*wkk+Rb_table(nb,kk+1)*(1-wkk)-Rgrid(0)
              dztemp=Zb_table(nb,kk)*wkk+Zb_table(nb,kk+1)*(1-wkk)-Zgrid(0)
              if((drtemp-x)**2+(dztemp-z)**2 < (Rbp-x)**2+(Zbp-z)**2) then
                 Rbp=drtemp
                 Zbp=dztemp
              end if
              kk=kk+nzeta_blob
           end do
          end do
!          write(*,*)m,Rbp+Rgrid(0),Zbp+Zgrid(0)!, x+Rgrid(0),z+Zgrid(0)

        
              

            !            Rbp=
!            tempb=abs((Rb_table(kk::(nzeta_blob/nb_temp))-x-Rgrid(0))**2+(Zb_table(kk::(nzeta_blob/nb_temp))-z-Zgrid(0))**2)
!            write(*,*) tempb(1),tempb(2),tempb(3),tempb(4),tempb(5),tempb(6)
!            ll=(minloc(tempb,1)-1)*nzeta_blob/nb_temp+kk
!            ll_plus_1=ll+1
!            write(*,*)ll,kk
!            write(*,*)m,(Rb_table(ll)),(Zb_table(ll))
!            write(*,*)m,x+Rgrid(0),z+Zgrid(0)

         
!            if (ll==nzeta_blob*ntern_blob+1) then
!               ll_plus_1=ll
!            end if

!            dldzetabp=dldzetab_table(ll)*wkk+dldzetab_table(ll_plus_1)*(1-wkk)
!            Rbp=Rb_table(ll)*wkk+Rb_table(ll_plus_1)*(1-wkk)-Rgrid(0)
!            Zbp=Zb_table(ll)*wkk+Zb_table(ll_plus_1)*(1-wkk)-Zgrid(0)
!            zetalbp=zetalb_table(ll)*wkk+zetalb_table(ll_plus_1)*(1-wkk)
!            write(*,*)zetalb_table(ll),wkk
!            write(*,*) Rb_table(4,kk),Rb_table(5,kk),Rb_table(6,kk)
!            write(*,*)Rbp,Zbp,zetalbp
!            write(*,*)zeta,kk,wkk

!           phibp = -Ablob0/dbZ**2*exp(-(x-Rbp)**2/(2*dbR**2)-(z-Zbp)**2/(2*dbZ**2)-zetalbp**2/(2*dbzeta**2))*(z-Zbp)
!           write(*,*) (mu(m)*b+0.5*mims(1)*u3(m)**2-q(1)*Ablob0/dbZ**2*exp(-(x-Rbp)**2/(2*dbR**2)-(z-Zbp)**2/(2*dbZ**2)-zetalbp**2/(2*dbzeta**2))*(z-Zbp)), -Ablob0/dbZ**2*exp(-(x-Rbp)**2/(2*dbR**2)-(z-Zbp)**2/(2*dbZ**2)-zetalbp**2/(2*dbzeta**2))*(z-Zbp)

!         end do
         

!  4 pt. avg. done explicitly for vectorization...
         !$acc loop seq
         do 200 l=1,lr(1)
!
            xt=x2(m)+rhox(l) !rwx(1,l)*rhog
            zt=z2(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
!            zeta=modulo(zeta2(m),pi2)

!            write(*,*)mu(m)
!   particle can go out of bounds during gyroavg...
            if( (xt<2*dxeq).or.(xt>lx-2*dxeq) ) xt=x2(m)
            if( (zt<2*dzeq).or.(zt>lz-2*dzeq) ) zt=z2(m)
!            xt=modulo(xs,xdim)
            !            zt=modulo(zt,zdim)

!            write(*,*)x2(m),xt,Rbp,z3(m),zt,Zbp,zeta2(m),zetalbp
            phibp=-Ablob0/dbZ**2*exp(-(xt-Rbp)**2/(2*dbR**2)-(zt-Zbp)**2/(2*dbZ**2)-zetalbp**2/(2*dbzeta**2))
!            ebparp=ebparp+(zetalbp)*phibp*(zt-Zbp)/(dldzetabp*dbzeta**2)
!            ebRp=ebRp+(xt-Rbp)*phibp*(zt-Zbp)/dbR**2
            ebZp=ebZp+((zt-Zbp)**2/dbZ**2-1)*phibp
!            ebzetap=ezetap-phibp/(xt+Rgrid(0))*((xt-Rbp)*(zt-Zbp)/dbR**2*(Rb_table(ll_plus_1)-Rb_table(ll))/dzetab+((zt-Zbp)**2/dbZ**2-1)*(Zb_table(ll_plus_1)-Zb_table(ll))/dzetab)!+(zetalbp/dbzeta**2)*(zt-Zbp)*phibp/(xt+Rgrid(0))
!            write(*,*)phibp, xt+Rgrid(0),(xt-Rbp)*(zt-Zbp)/dbR**2*((Rb_table(ll_plus_1)-Rb_table(ll))/dzetab,((zt-Zbp)**2/dbZ**2-1)*(Zb_table(ll_plus_1)-Zb_table(ll)/dzetab)
            i=int(xt/dx)
            j=int(zt/dz)
            k=int(zeta/dzeta)

!           write(*,*)ebparp,ebRp,ebZp,ebzetap,phibp*(zt-Zbp)
            
!            write(*,*)phibp
 
!          write(*,*)i,j,k  


            wx0=float(i+1)-xt/dx
            wx1=1.-wx0
            wy0=float(j+1)-zt/dz
            wy1=1.-wy0
            wz0=float(k+1)-zeta/dzeta
            wz1=1.-wz0
!            write(*,*)wx0,wx1,wy0,wy1,wz0,wz1
!            write(*,*)'ppush,ez=',ez(i,j,k),ez(i,j,k+1),ez(i+1,j,k),ez(i+1,j,k+1)

              k_plus_1=k+1
              if(k==kmx) k_plus_1=0
              !              phib=-phib0*(zt-zb)
              
              
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

!         ebparp=ebparp/4
!         ebRp=ebRp/4
         ebZp=ebZp/4
!         ebzetap=ebzetap/4

!write(*,*)'ppush,Ebz=',ebZp,'Esz=',ezp,'EsR=',exp1
         
!         exp1=exp1+ebRp
!         ezp=ezp+ebZp
!         ezetap=ezetap+ebzetap

!         write(*,*)exp1,ezp,ezetap
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

         !         write(*,*) vpar

         pzdot = (Bstar3(1)*(q(1)*exp1-mu(m)*dbdxp)+Bstar3(2)*(q(1)*ezp-mu(m)*dbdzp)+Bstar3(3)*(q(1)*ezetap-mu(m)*dbdzetap))/(mims(1)*bstar)

         ezp = ezp+ ebZp
         xdot = (vpar*Bstar3(1)+(mu(m)*(bfldzp*dbdzetap-bfldzetap*dbdzp)/q(1)+((ezp)*bfldzetap-(ezetap)*bfldzp))/(b))/bstar
         zdot = (vpar*Bstar3(2)+(mu(m)*(bfldzetap*dbdxp-bfldxp*dbdzetap)/q(1)+ ((ezetap)*bfldxp-(exp1)*bfldzetap))/(b))/bstar
         zetadot = (vpar*Bstar3(3)+(mu(m)*(bfldxp*dbdzp-bfldzp*dbdxp)/q(1)+((exp1)*bfldzp-(ezp)*bfldxp))/(b))/bstar/(x2(m)+xctr-xdim*0.5)

!         write(*,*)bfldxp,bfldzp,bfldzetap
!write(*,*)((ezp)*bfldzetap-(ezetap)*bfldzp)/(b)/bstar, ((ezetap)*bfldxp-(exp1)*bfldzetap)/(b)/bstar,((exp1)*bfldzp-(ezp)*bfldxp)/(b)/bstar/(x2(m)+xctr-xdim*0.5)
         
 !        xdot = (vpar*Bstar3(1)+(mu(m)*(bfldzp*dbdzetap-bfldzetap*dbdzp)/q(1)+((ezp)*bfldzetap-(ezetap)*bfldzp))/(b))/bstar
 !        zdot = (vpar*Bstar3(2)+(mu(m)*(bfldzetap*dbdxp-bfldxp*dbdzetap)/q(1)+ ((ezetap)*bfldxp-(exp1)*bfldzetap))/(b))/bstar
 !        zetadot = (vpar*Bstar3(3)+(mu(m)*(bfldxp*dbdzp-bfldzp*dbdxp)/q(1)+((exp1)*bfldzp-(ezp)*bfldxp))/(b))/bstar/(x2(m)+xctr-xdim*0.5)
         
!         pzd0 = -mu(m)/mims(1)/b*(bfldxp*dbdxp+bfldzp*dbdzp)
!         write(*,*)
!         write(*,*)(exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)/b*q(1)/mims(1)
!         pzdot = pzd0+(exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)/b*(q(1)/mims(1)+bdcurlbp*vpar/b)*nonlin
!         exp1=exp1+ebRp
!         ezp = ezp+ ebZp
!         ezetap=ezetap+ebzetap
         !          pzdot = pzd0+((exp1*bfldxp+ezp*bfldzp+ezetap*bfldzetap)*q(1)/mims(1)+vcurlbdotE)/bstar*nonlin
         
!         pzdot = (Bstar3(1)*(q(1)*exp1-mu(m)*dbdxp)+Bstar3(2)*(q(1)*ezp-mu(m)*dbdzp)+Bstar3(3)*(q(1)*ezetap-mu(m)*dbdzetap)+q(1)*(ebparp*bstar+(Bstar3(1)-bfldxp)*ebRp+(Bstar3(2)-bfldzp)*ebZp+(Bstar3(3)-bfldzetap)*ebzetap))/(mims(1)*bstar)

!          pzdot = (Bstar3(1)*(q(1)*exp1-mu(m)*dbdxp)+Bstar3(2)*(q(1)*ezp-mu(m)*dbdzp)+Bstar3(3)*(q(1)*ezetap-mu(m)*dbdzetap)+q(1)*(ebparp*bstar))/(mims(1)*bstar)
          
!          write(*,*)exp1,ezp,ezetap          
!          write(*,*)(Bstar3(1)*(q(1)*exp1-mu(m)*dbdxp)+Bstar3(2)*(q(1)*ezp-mu(m)*dbdzp)+Bstar3(3)*(q(1)*ezetap-mu(m)*dbdzetap))/(mims(1)*bstar),(ebparp*b+(Bstar3(1)-bfldxp)*ebRp+(Bstar3(2)-bfldzp)*ebZp+(Bstar3(3)-bfldzetap)*ebzetap)*q(1)/(mims(1)*bstar)
         
         edot = q(1)*(xdot*exp1+zdot*ezp+zetadot*ezetap)

         x3(m) = x2(m) + 0.5*dt*xdot
         z3(m) = z2(m) + 0.5*dt*zdot
         zeta3(m) = zeta2(m) + 0.5*dt*zetadot

         u3(m) = u2(m) + 0.5*dt*pzdot

!         write(*,*)x3(m),z3(m)
!        dum = 1.0
!        vxdum = (ezp/b+vpar/b*delbxp)*dum1
!        vzdum = (-exp1/b+vpar/b*delbzp)*dum1
!         vxdum = eyp+vpar/b*delbxp
!         w3(m)=w2(m) + 0.5*dt*(vxdum*kapxp + vzdum*kapzp+edot/ter)*dum*xnp
         w3(m)=w2(m)
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
!   Ion corrector push with electrostatic potential and energy-conserving random walk
!-----------------------------------------------------------------------
      use gemx_com
      use equil
      implicit none

      REAL(8) :: exp1,ezp,ezetap,delbxp,delbzp,energy, energy0,nudi0,nudi,T_center,ni_temp,phibp,Rbp,Zbp,zetalbp,ebRp,ebZp,ebzetap
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,dum,vxdum,vzdum,dum1,vzetadum,psiparticle,phi_particle,wkk,dldzetabp,ebparp,dzetab
      INTEGER :: m,i,j,k,l,n,k_plus_1,kk,ll_plus_1,ll,nb,nb_temp=10
      
      REAL(8) :: rhog,vfac,kapxp,kapzp,vpar,pidum,kaptxp,kapnxp,kaptzp,kapnzp,xnp,bdcurlbp
      REAL(8) :: b,th,r,enerb,qr,ter,x,z,zeta
      REAL(8) :: xt,xs,zt,xdot,zdot,zetadot,xdt,ydt,pzdot,edot,pzd0,vp0,vcurlbdotE
      REAL(8) :: dbdxp,dbdzp,bfldp,bfldxp,bfldzp,bfldzetap, bstar, dbdzetap=0
      REAL(8) :: rhox(4),rhoy(4),psp,pzp,curlbp(3),Bstar3(3),tempb(60),drtemp,dztemp

      ! 用于随机游走、电势插值与能量守恒修正的局部变量
      REAL(8) :: E_total_rw, bfldp_new, E_mag_new, E_kin_par_new, phi_old, phi_new
      REAL(8) :: wx0_new, wx1_new, wy0_new, wy1_new, wz0_new, wz1_new, xt_rw, zt_rw, zeta_rw
      INTEGER :: i_new, j_new, k_new, rand_idx

      start_cpush_tm = MPI_WTIME()
      nudi0 = 1/sqrt(2.0)*18.4*e**1.5*4.7140d-8*1.d-6

      dzetab=2*pi/nzeta_blob

!$acc parallel loop gang vector private(Bstar3,rhoy,rhox,curlbp,tempb, &
!$acc& E_total_rw,bfldp_new,E_mag_new,E_kin_par_new,phi_old,phi_new, &
!$acc& wx0_new,wx1_new,wy0_new,wy1_new,wz0_new,wz1_new,xt_rw,zt_rw,zeta_rw, &
!$acc& i_new,j_new,k_new,rand_idx) copy(rand_table)
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

         exp1=0.
         ezp=0.
         ezetap=0.
         delbxp=0.
         delbzp=0.

         ebRp=0
         ebZp=0

         Rbp=10000
         Zbp=10000

         zeta=modulo(zeta3(m),2*pi)

         do nb=1,nblob
            wkk=(modulo(zeta-zeta0_blob(nb),2*pi)/dzetab)
            kk=int(wkk)+1
            wkk=kk-wkk
            do ll=0,ntern_blob-1
               dldzetabp=dldzetab_table(kk)*wkk+dldzetab_table(kk+1)*(1-wkk)             
               drtemp=Rb_table(nb,kk)*wkk+Rb_table(nb,kk+1)*(1-wkk)-Rgrid(0)
               dztemp=Zb_table(nb,kk)*wkk+Zb_table(nb,kk+1)*(1-wkk)-Zgrid(0)
               if((drtemp-x)**2+(dztemp-z)**2 < (Rbp-x)**2+(Zbp-z)**2) then
                  Rbp=drtemp
                  Zbp=dztemp
               end if
               kk=kk+nzeta_blob
            end do
         end do

         !$acc loop seq
         do 200 l=1,lr(1)
            xt=x3(m)+rhox(l)
            zt=z3(m)+rhoy(l)

            if( (xt<2*dxeq).or.(xt>lx-2*dxeq) ) xt=x3(m)
            if( (zt<2*dzeq).or.(zt>lz-2*dzeq) ) zt=z3(m)

            phibp=-Ablob0/dbZ**2*exp(-(xt-Rbp)**2/(2*dbR**2)-(zt-Zbp)**2/(2*dbZ**2)-zetalbp**2/(2*dbzeta**2))
            ebZp=ebZp+((zt-Zbp)**2/dbZ**2-1)*phibp

            i=int(xt/dx)
            j=int(zt/dz)
            k=int(zeta/dzeta)

            wx0=float(i+1)-xt/dx
            wx1=1.-wx0
            wy0=float(j+1)-zt/dz
            wy1=1.-wy0
            wz0=float(k+1)-zeta/dzeta
            wz1=1.-wz0

            p_psi(m)=wx0*wy0*psi_p(i,j)+wx1*wy0*psi_p(i+1,j)+wx0*wy1*psi_p(i,j+1)+wx1*wy1*psi_p(i+1,j+1)

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
200      continue
         exp1 = exp1/4.
         ezp = ezp/4.
         ezetap = ezetap/4.
         delbxp = delbxp/4.
         delbzp = delbzp/4.

         ebZp=ebZp/4

         vfac = 0.5*(mims(1)*u2(m)**2 + 2.*mu(m)*b)
         kapxp = kapnxp - (1.5-vfac/ter)*kaptxp
         kapzp = kapnzp - (1.5-vfac/ter)*kaptzp         

         vpar = u3(m)
         enerb=(mu(m)+mims(1)*vpar*vpar/b)/q(1)*tor

         Bstar3(1)=bfldxp+mims(1)*vpar*curlbp(1)/q(1)+delbxp
         Bstar3(2)=bfldzp+mims(1)*vpar*curlbp(2)/q(1)+delbzp
         Bstar3(3)=bfldzetap+mims(1)*vpar*curlbp(3)/q(1)

         bstar=(bfldxp*Bstar3(1)+bfldzp*Bstar3(2)+bfldzetap*Bstar3(3))/bfldp

         dum1 = 1.
         pzdot = (Bstar3(1)*(q(1)*exp1-mu(m)*dbdxp)+Bstar3(2)*(q(1)*ezp-mu(m)*dbdzp)+Bstar3(3)*(q(1)*ezetap-mu(m)*dbdzetap))/(mims(1)*bstar)

         ezp = ezp+ ebZp

         xdot = (vpar*Bstar3(1)+(mu(m)*(bfldzp*dbdzetap-bfldzetap*dbdzp)/q(1)+((ezp)*bfldzetap-(ezetap)*bfldzp))/(b))/bstar
         zdot = (vpar*Bstar3(2)+(mu(m)*(bfldzetap*dbdxp-bfldxp*dbdzetap)/q(1)+ ((ezetap)*bfldxp-(exp1)*bfldzetap))/(b))/bstar
         zetadot = (vpar*Bstar3(3)+(mu(m)*(bfldxp*dbdzp-bfldzp*dbdxp)/q(1)+((exp1)*bfldzp-(ezp)*bfldxp))/(b))/bstar/(x2(m)+xctr-xdim*0.5)

         edot = q(1)*(xdot*exp1+zdot*ezp+zetadot*ezetap)

         ! 1. 确定性推进 (Deterministic Push)
         x3(m) = x2(m) + dt*xdot
         z3(m) = z2(m) + dt*zdot
         zeta3(m) = zeta2(m) + dt*zetadot
         u3(m) = u2(m) + dt*pzdot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RANDOM WALK WITH ENERGY CONSERVATION (INC. PHI) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
         if (random_walk == 1) then
             ! A. 插值获取游走前 (x3, z3, zeta3) 处的静场电势 phi_old
             xt_rw = min(max(x3(m), 0.0d0), lx - dx)
             zt_rw = min(max(z3(m), 0.0d0), lz - dz)
             zeta_rw = modulo(zeta3(m), 2.0d0*pi)

             i_new = int(xt_rw / dx)
             i_new = min(i_new, imx - 1)
             wx0_new = float(i_new + 1) - xt_rw / dx
             wx1_new = 1.0d0 - wx0_new

             j_new = int(zt_rw / dz)
             j_new = min(j_new, jmx - 1)
             wy0_new = float(j_new + 1) - zt_rw / dz
             wy1_new = 1.0d0 - wy0_new

             k_new = int(zeta_rw / dzeta)
             k_new = min(k_new, kmx)
             wz0_new = float(k_new + 1) - zeta_rw / dzeta
             wz1_new = 1.0d0 - wz0_new

             k_plus_1 = k_new + 1
             if (k_new == kmx) k_plus_1 = 0

             phi_old = wx0_new*wy0_new*wz0_new*phi(i_new,j_new,k_new) + wx1_new*wy0_new*wz0_new*phi(i_new+1,j_new,k_new) &
                     + wx0_new*wy1_new*wz0_new*phi(i_new,j_new+1,k_new) + wx1_new*wy1_new*wz0_new*phi(i_new+1,j_new+1,k_new) + &
                       wx0_new*wy0_new*wz1_new*phi(i_new,j_new,k_plus_1) + wx1_new*wy0_new*wz1_new*phi(i_new+1,j_new,k_plus_1) + &
                       wx0_new*wy1_new*wz1_new*phi(i_new,j_new+1,k_plus_1) + wx1_new*wy1_new*wz1_new*phi(i_new+1,j_new+1,k_plus_1)

             ! B. 记录游走前的总能量 (平行动能 + 磁能 + 静电势能 q*Phi)
             E_total_rw = 0.5d0 * mims(1) * (u3(m)**2) + mu(m) * bfldp + q(1) * phi_old

             ! C. 取出粒子随机数 (防 OpenACC 数据竞态)
             rand_idx = mod(m + globle_integer, 10007) + 1

             ! D. 垂直磁场方向坐标随机游走
             x3(m) = x3(m) + rand_table(rand_idx) * lambda_n * sqrt(2.0d0 * dt * abs(u2(m)) / (15.0d0 * (bfldxp**2 + bfldzp**2))) * bfldzp
             z3(m) = z3(m) - rand_table(rand_idx) * lambda_n * sqrt(2.0d0 * dt * abs(u2(m)) / (15.0d0 * (bfldxp**2 + bfldzp**2))) * bfldxp

             ! E. 重新插值计算新位置 (x3, z3, zeta3) 处的磁场 B_new 和电势 phi_new
             xt_rw = min(max(x3(m), 0.0d0), lx - dx)
             zt_rw = min(max(z3(m), 0.0d0), lz - dz)
             
             ! 1) 插值新磁场 B_new (使用平衡场网格 dxeq, dzeq)
             i_new = int(xt_rw / dxeq)
             i_new = min(i_new, nx - 1)
             wx0_new = (i_new + 1) - xt_rw / dxeq
             wx1_new = 1.0d0 - wx0_new

             j_new = int(zt_rw / dzeq)
             j_new = min(j_new, nz - 1)
             wy0_new = (j_new + 1) - zt_rw / dzeq
             wy1_new = 1.0d0 - wy0_new

             bfldp_new = wx0_new * wy0_new * b0(i_new, j_new) + wx0_new * wy1_new * b0(i_new, j_new + 1) &
                       + wx1_new * wy0_new * b0(i_new + 1, j_new) + wx1_new * wy1_new * b0(i_new + 1, j_new + 1)

             ! 2) 插值新电势 phi_new (使用场网格 dx, dz, dzeta)
             i_new = int(xt_rw / dx)
             i_new = min(i_new, imx - 1)
             wx0_new = float(i_new + 1) - xt_rw / dx
             wx1_new = 1.0d0 - wx0_new

             j_new = int(zt_rw / dz)
             j_new = min(j_new, jmx - 1)
             wy0_new = float(j_new + 1) - zt_rw / dz
             wy1_new = 1.0d0 - wy0_new

             phi_new = wx0_new*wy0_new*wz0_new*phi(i_new,j_new,k_new) + wx1_new*wy0_new*wz0_new*phi(i_new+1,j_new,k_new) &
                     + wx0_new*wy1_new*wz0_new*phi(i_new,j_new+1,k_new) + wx1_new*wy1_new*wz0_new*phi(i_new+1,j_new+1,k_new) + &
                       wx0_new*wy0_new*wz1_new*phi(i_new,j_new,k_plus_1) + wx1_new*wy0_new*wz1_new*phi(i_new+1,j_new,k_plus_1) + &
                       wx0_new*wy1_new*wz1_new*phi(i_new,j_new+1,k_plus_1) + wx1_new*wy1_new*wz1_new*phi(i_new+1,j_new+1,k_plus_1)

             ! F. 能量重标定 (Energy Re-scaling): 扣除新位置的磁能和电势能后更新平行速度 u3
             E_mag_new = mu(m) * bfldp_new
             E_kin_par_new = E_total_rw - E_mag_new - q(1) * phi_new

             if (E_kin_par_new > 0.0d0) then
                 u3(m) = sign(1.0d0, u3(m)) * sqrt(2.0d0 * E_kin_par_new / mims(1))
             else
                 ! 若游走到强磁场或高电势区导致平行动能耗尽，反转平行速度方向（反弹）
                 u3(m) = -u3(m)
             end if
         end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END OF RANDOM WALK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         zeta3(m)=modulo(zeta3(m),pi2)

         if(w2(m)>1e-100.and. ((x3(m)>2*dxeq).and.(x3(m)<lx-2*dxeq).and.(z3(m)>2*dzeq).and.(z3(m)<lz-2*dzeq)) ) then
            u2(m)=u3(m)
            x2(m)=x3(m)
            z2(m)=z3(m)
            zeta2(m)=zeta3(m)
            w2(m)=w3(m)
         else
            u3(m)=u2(m)
            x3(m)=x2(m)
            z3(m)=z2(m)
            zeta3(m)=zeta2(m)
            w2(m)=0.
            w3(m)=0.
         endif

      enddo

      globle_integer = mod(globle_integer + 1, 10007)

      end_cpush_tm = MPI_WTIME()
      cpush_tm = cpush_tm + end_cpush_tm - start_cpush_tm
      return
      end subroutine cpush
