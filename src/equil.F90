MODULE equil
      IMPLICIT NONE
      real :: mimp=2,chgi=1
      real :: beta,rmaj0,a,q0,r0,q0p,q0abs,shat0
      real :: dR,dth,mu0,e,proton
      integer :: nr=200,nr2=100,ntheta=200,isgnf=1,isgnq=-1,isupae0=0,tor_n
      real :: psi_max=0.31, psi_min=-0.1 ,R_min=1.0, Z_min=-1.5, Z_internal=-1.2, psi_div=0.305,psi_a=0.311647

!     GEM-X
      !integer :: cont=259002 !Calder Edit

      integer :: nx=449,nz=433,nzeta=32
      real :: xdim,zdim,xctr,zctr,dxeq,dzeq
      
      real,dimension(:,:),allocatable :: b0,b0x,b0z,b0zeta,dbdx,dbdz,c2_over_vA2
      real,dimension(:,:),allocatable :: t0i,t0e,xn0i,xn0e,captix,captex,capnix,capnex,captiz,captez,capniz,capnez
      
      real,dimension(:),allocatable :: psi,f,psip,sf, &
                                      vpari,vparip, &
                                      zeff,nue0,phinc,phincp,&
                                      er,upari, Rgrid, Zgrid
      real,dimension(:,:),allocatable :: t0s,xn0s,capts,capns,vpars,vparsp,psi_p,mask,mask2,mask3,mask4
      real:: bu,tu,nu,xu,frequ,vu,eru

!     for including bstar effects
      real,dimension(:),allocatable :: psip2
      real,dimension(:,:),allocatable :: curvbz,srbr,srbz,thbr,thbz,prsrbr,prsrbz,pthsrbr,pthsrbz,bdcrvb

!     for equilibrium current                                                                                                                                         
      real(8),dimension(:,:),allocatable :: upae0,nuob,dnuobdr,dnuobdt

      real,dimension(:,:,:),allocatable ::curlb
      
      
!     for phi average Calder Edit
      integer :: num_lines, line
      real,dimension(:,:),allocatable :: phiavg
      real,dimension(:),allocatable :: psitab,weight00,weight10,weight01,weight11,jacobian,deno !Calder Edits
      integer, dimension(:), allocatable :: gindex,iarray,jarray,priv
       

contains
      subroutine new_equil()
      use gemx_com,only: myid
      implicit none
      real(8) :: pi,pi2,r,th,s

      integer :: i,j,k,m,i1,j1,j2
      real(8) :: dum,x,tempn
      real(8) :: omegau          !,e,proton
 
!     global equilibrium data

      allocate(sf(0:nr), psi(0:nr),&
         psip(0:nr),&
         zeff(0:nr),nue0(0:nr),&
         vpari(0:nr),phinc(0:nr), &
         vparip(0:nr),phincp(0:nr),er(0:nr), &
         upari(0:nr))

      allocate(t0s(1:5,0:nr),xn0s(1:5,0:nr),&
         capts(1:5,0:nr),capns(1:5,0:nr),vpars(1:5,0:nr),&
         vparsp(1:5,0:nr))

      allocate(upae0(0:nr,0:ntheta),nuob(0:nr,0:ntheta),dnuobdr(0:nr,0:ntheta),dnuobdt(0:nr,0:ntheta))      

!     GEM-X
      allocate(phiavg(0:nx,0:nz)) !Calder Edit

      allocate(b0(0:nx,0:nz),b0x(0:nx,0:nz),b0z(0:nx,0:nz),b0zeta(0:nx,0:nz),dbdx(0:nx,0:nz),dbdz(0:nx,0:nz))
      allocate(t0i(0:nx,0:nz),t0e(0:nx,0:nz),xn0i(0:nx,0:nz),xn0e(0:nx,0:nz), &
               capnix(0:nx,0:nz),capnex(0:nx,0:nz),captix(0:nx,0:nz),captex(0:nx,0:nz), &
               capniz(0:nx,0:nz),capnez(0:nx,0:nz),captiz(0:nx,0:nz),captez(0:nx,0:nz))      

      allocate(curlb(0:nx,0:nz,3))
      
      allocate(Rgrid(0:nx),Zgrid(0:nz),bdcrvb(0:nx,0:nz))
      allocate(psi_p(0:nx,0:nz),mask(0:nx,0:nz),c2_over_vA2(0:nx,0:nz),mask2(0:nx,0:nz),mask3(0:nx,0:nz),mask4(0:nx,0:nz))

 
      
                  
               
      open(unit=10, file = 'R.dat',status='old',action='read')
	read(10,*) Rgrid
      close(10)    

      
       open(unit=10, file = 'Z.dat',status='old',action='read')
	read(10,*) Zgrid
        close(10)

        open(unit=10, file='psi_p.dat',status='old',action='read')
        read(10,*) psi_p
        close(10)
        
        !Calder Edits Start
        !character(len=100) :: line_buffer
         ! Open the dataset file
      open(unit=10, file='jacodata.dat', status='old', action='read')
      ! Determine the number of lines in the dataset file
      num_lines = 80817
      ! do
      !     read(10, '(A)')!,end=100)! line_buffer
      !     num_lines = num_lines + 1
      ! enddo
      ! rewind(10)

      ! Allocate arrays based on number of lines
      allocate(gindex(num_lines), psitab(num_lines), iarray(num_lines), jarray(num_lines), &
               weight00(num_lines), weight10(num_lines), weight01(num_lines), weight11(num_lines), &
               jacobian(num_lines), deno(num_lines), priv(num_lines))
  
      ! Read data from the dataset file into arrays
      do line = 1, num_lines
          read(10, *) gindex(line), psitab(line), iarray(line), jarray(line), &
                      weight00(line), weight10(line), weight01(line), weight11(line), &
                       jacobian(line), deno(line), priv(line)
      enddo
      !write(*,*) psitab
      close(10)
        !Calder Edits Finish



!        open(unit=11, file = 'psi_test',status='unknown',action='write')
 !     do i=0,nx
!         write(11,*) psi_p
 !     end do
 !     close(11)     
      
            
        
        xdim=Rgrid(nx)-Rgrid(0)
        zdim=Zgrid(nz)-Zgrid(0)
        dR=(Rgrid(nx)-Rgrid(0))/nx
!        dx=(Rgrid(nx)-Rgrid(0))/nx
!        dZ=(Zgrid(nz)-Zgrid(0))/nz
      
      open(unit=10, file = 'BR.dat',status='old',action='read')
	read(10,*) b0x
      close(10)    

      
!      open(unit=11, file = 'debug.dat',status='unknown',action='write')
!      write(11,*) b0x(5,0), 'dR='dR
!      close(11)
      
      open(10, file = 'Bz.dat',status='old',action='read')
	read(10,*) b0z
      close(10)
      
      open(unit=10, file = 'Bt.dat',status='old',action='read')
	read(10,*) b0zeta
      close(10)

!      B0x=0
!      B0z=0
      
      do i =0,nx
            do j =0,nz
		b0(i,j) = sqrt(b0x(i,j)**2+b0z(i,j)**2+b0zeta(i,j)**2)
    	end do
      end do 




      e = 1.6e-19
      mu0 = 1.25663706212e-6
      proton = 1.67e-27
      Bu = 1.98562

      Bu = 1
      
      Tu = 1000*e
      omegau = e*Bu/proton
      frequ = omegau

      vu = sqrt(Tu/proton)
      vu=1
      xu = proton*vu/(e*Bu)
      xu=1
      nu = 2.5e19
      beta = 4*3.14159*1e-7*nu*Tu/Bu**2


      !     assign T, n profiles 
      do i = 0,nx
         do j = 0,nz
            t0i(i,j) = 1.*Tu
            t0e(i,j) = 1.*Tu
            xn0i(i,j) = 1.*nu
            xn0e(i,j) = 1.*nu
         end do
      end do

      !Calder Edit: realistic profilies for ITG runs
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! open(unit=10, file = 'Profiles/ti0_profile.dat',status='old',action='read')
      ! read(10,*) t0i
      ! close(10)

      ! open(unit=10, file = 'Profiles/te0_profile.dat',status='old',action='read')
      ! read(10,*) t0e
      ! close(10)

      ! open(unit=10, file = 'Profiles/ni0_profile.dat',status='old',action='read')
      ! read(10,*) xn0i
      ! close(10)

      ! open(unit=10, file = 'Profiles/ne0_profile.dat',status='old',action='read')
      ! read(10,*) xn0e
      ! close(10)

      !Manufactured Temperature Profiles: 0.5*(-np.tanh(5*np.array(profiles['psinorm'])-2.5)+1)*np.max(np.array(profiles['ti']))
      open(unit=10, file = 'tif_profile.dat',status='old',action='read')
      read(10,*) t0i
      close(10)

      open(unit=10, file = 'tef_profile.dat',status='old',action='read')
      read(10,*) t0e
      close(10)

      !Put into SI units
      t0i  = t0i*Tu
      t0e  = t0e*Tu
      ! xn0i = xn0i*1.0e21
      ! xn0e = xn0e*1.0e21
      ! xn0i = xn0e
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      open(unit=10, file = 'ne0.dat',status='old',action='read')
	read(10,*) xn0e
      close(10)
 

      open(unit=10, file = 'ti0.dat',status='old',action='read')
	read(10,*) t0i
      close(10)

      xn0e=xn0e*nu
      t0i=t0i*Tu
      xn0i=xn0e


      do i=0,nx
         do j=0,nz
            if (xn0e(i,j)<1e-9) then
               c2_over_vA2(i,j)=mu0*2*proton*0.01*xn0e(nx/2,nz/2)/(b0(i,j)**2)*vu**2
            else
               c2_over_vA2(i,j)=mu0*2*proton*xn0e(i,j)/(b0(i,j)**2)*vu**2!2*Rgrid(i)**2/(Rgrid(0)+Rgrid(nx))**2
            end if
            
           ! if (c2_over_vA2(i,j)<1.0e-19) c2_over_vA2(i,j)=0.01*mu0*2*proton*xn0e(i/2,j/2)/(b0(i/2,j/2)**2)*vu**2
         end do
      end do
      
            
      mask=0
!      c2_over_vA2=1
          do i=2,nx-2
           do j=2,nz-2
             ! c2_over_vA2(i,j)=mu0*2*proton*xn0e(i,j)/(b0(i,j)**2)*vu**2!2*Rgrid(i)**2/(Rgrid(0)+Rgrid(nx))**2
!              write(*,*) c2_over_vA2(i,j)
              if (psi_p(i,j)<psi_max .and.psi_p(i,j)>psi_min .and. (Zgrid(j)>Z_internal .or. (Zgrid(j)>Z_min .and. psi_p(i,j)>psi_div)) .and. Rgrid(i)>R_min)  then
!             if(psi_p(i,j)<0.3 .and. Zgrid(j)>-1.08) then
              mask(i,j)=1
           else
              mask(i,j)=0
           endif
        end do
     end do

     do i=2,nx-2
        do j=2,nz-2
             if(mask(i+1,j)*mask(i-1,j)*mask(i,j+1)*mask(i,j-1)==0) then
                mask2(i,j)=0
             else
              mask2(i,j)=2
           endif
        end do
     end do

     do i=2,nx-2
        do j=2,nz-2
             if(mask2(i+1,j)*mask2(i-1,j)*mask2(i,j+1)*mask2(i,j-1)==0) then
                mask3(i,j)=0
             else
              mask3(i,j)=3
           endif
        end do
     end do

     do i=2,nx-2
        do j=2,nz-2
             if(mask3(i+1,j)*mask3(i-1,j)*mask3(i,j+1)*mask3(i,j-1)==0) then
                mask4(i,j)=0
             else
              mask4(i,j)=4
           endif
        end do
     end do



  if(myid==0)then
    open(unit=11, file = 'test_1_over_vA2',status='unknown',action='write')
      do j=0,nz
         write(11,*) c2_over_vA2(:,j)
      end do
      close(11)
     open(unit=11, file = 'mask',status='unknown',action='write')
      do j=0,nz
         write(11,*) mask(:,j)
      end do
      close(11)

      open(unit=11, file = 'mask2',status='unknown',action='write')
      do j=0,nz
         write(11,*) mask2(:,j)
      end do
      close(11)  
     open(unit=11, file = 'mask3',status='unknown',action='write')
      do j=0,nz
         write(11,*) mask3(:,j)
      end do
      close(11)
      open(unit=11, file = 'mask4',status='unknown',action='write')
      do j=0,nz
         write(11,*) mask4(:,j)
      end do
      close(11)
   end if
   
!     Normalization

      b0 = b0/Bu
      b0x = b0x/Bu
      b0z = b0z/Bu
      b0zeta = b0zeta/Bu


      
!     write(*,*) 'betaU=', beta

      pi = atan(1.0)*4      
      rmaj0 = 1000.
      a = 360.

!      xctr = a*1.5
!      zctr = 0.


      xctr =  0.5*(Rgrid(nx)+Rgrid(0))/xu
      zctr  = 0.5*(Zgrid(nz)+Zgrid(0))/xu
      xdim=xdim/xu
      zdim=zdim/xu
 !     xdim = a*2
 !     zdim = a*3
      dxeq = xdim/nx
      dzeq = zdim/nz



      do i=1,nx-1
         do j=1,nz-1
            curlb(i,j,1)=-((-b0zeta(i,j+1)/b0(i,j+1)+b0zeta(i,j-1)/b0(i,j-1))/(2*dzeq))
            curlb(i,j,2)=-(1/Rgrid(i)*(Rgrid(i+1)*b0zeta(i+1,j)/b0(i+1,j)-Rgrid(i-1)*b0zeta(i-1,j)/b0(i-1,j))/(2*dxeq))
            curlb(i,j,3)=-((b0x(i,j+1)/b0(i,j+1)-b0x(i,j-1)/b0(i,j-1))/(2*dzeq)-(b0z(i+1,j)/b0(i+1,j)-b0z(i-1,j)/b0(i-1,j))/(2*dxeq))
!            write(*,*) curlb(i,j,2)+(b0zeta(i,j)/b0(i,j)/Rgrid(i)*xu+(b0zeta(i+1,j)/b0(i+1,j)-b0zeta(i-1,j)/b0(i-1,j))/(2*dxeq))
         end do
      end do
      
      

      do i=1,nx-1
         do j=1,nz-1
            bdcrvb(i,j)=-1/b0(i,j)*(b0x(i,j)*(-b0zeta(i,j+1)/b0(i,j+1)+b0zeta(i,j-1)/b0(i,j-1))/(2*dzeq)&
                                              +b0z(i,j)/Rgrid(i)*(Rgrid(i+1)*b0zeta(i+1,j)/b0(i+1,j)-Rgrid(i-1)*b0zeta(i-1,j)/b0(i-1,j))/(2*dxeq)&
                                              +b0zeta(i,j)*((b0x(i,j+1)/b0(i,j+1)-b0x(i,j-1)/b0(i,j-1))/(2*dzeq)-(b0z(i+1,j)/b0(i+1,j)-b0z(i-1,j)/b0(i-1,j))/(2*dxeq)))
         end do
      end do
      
      
      
      
!     Gradients of T and n profiles for global simulation. Flux-tube parameters done at the end
      captix = 0.
      captex = 0.
      capnix = 0.
      capnex = 0.
      captiz = 0.
      captez = 0.
      capniz = 0.
      capnez = 0.

      dbdx = 0.
      dbdz = 0.
      do i = 1,nx-1
         do j = 1,nz-1
            dbdx(i,j) = (b0(i+1,j)-b0(i-1,j))/(2*dxeq)
            dbdz(i,j) = (b0(i,j+1)-b0(i,j-1))/(2*dzeq)            
         end do
      end do

    if(myid==0)then
       open(11,file='xpp',status='replace')
       write(11,*)'xu,omegau,vu = ', xu, omegau, vu

10     format(1x,i5,10(1x,e16.9))
       do i = 0,nx
          write(11,10)i,b0(i,nz/2),b0x(i,nz/2),b0z(i,nz/2),b0zeta(i,nz/2),dbdx(i,nz/2),dbdz(i,nz/2)
       end do
       close(11)
      end if


 !     open(unit=11, file = 'debug.dat',status='unknown',action='write')
 !     write(11,*) dR
 !     close(11)

      
  end subroutine new_equil

            
END MODULE equil
