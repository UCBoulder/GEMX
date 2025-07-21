module global
	use inrtype, only: sp
	implicit none

	integer, parameter :: kr=4,kz=4,kp=4
	real(sp), allocatable :: knotr(:),knotz(:)

	real(sp), allocatable :: bscoefscal1(:,:),bscoefscal2(:,:),bscoefscal3(:,:),bscoefscal4(:,:), &
				bscoefpsi(:,:),bscoefbphi(:,:)
	real(sp), allocatable :: bscoefvec1r(:,:),bscoefvec1z(:,:),bscoefvec1p(:,:), &
				bscoefvec2r(:,:),bscoefvec2z(:,:),bscoefvec2p(:,:), &
				bscoefvec3r(:,:),bscoefvec3z(:,:),bscoefvec3p(:,:)
	real(sp), allocatable :: results(:,:)
	! physics constant
	!real(sp), parameter :: mp = 1.6726485d-27, me=9.1094d-31,ee = 1.6021766208d-19, mu0 = 4._dp*pi*1.d-7
end module global

module vars
	use opes
	use inrtype
	implicit none
	integer :: nxefit,nyefit
	real(sp), allocatable :: fold(:,:) ,gridpsi(:),gridr(:),gridz(:)

	integer :: values(1:8), k
	integer, dimension(:), allocatable :: seed

	real(sp) :: R_limiter, Z_limiter
	real(sp) tmpfpol,tmp1,tmp2
	real(sp) xdim,zdim,rcentr,rgrid1,zmid,rmagx,zmagx,simagx,sibdry,bcentr,cpasma,xdum
	real(sp), allocatable :: fpol(:), pres(:),workk1(:),workk2(:), qpsi(:)
	real(sp), allocatable :: knotpsi(:)
	real(sp), allocatable :: cscoeffpol(:,:)
	real(sp), allocatable :: bscoefscalb(:,:)
	type(vector), allocatable :: vec_b(:,:), vec_gradb(:,:), vec_curlb(:,:), vec_bxgradb(:,:)
	real(sp), allocatable :: scal_b(:,:), scal2(:,:),scal3(:,:),scal4(:,:)
	real(sp) x,y(4),dydx(4),yout(4)
	real(sp) rmag,bmag,vmag,tmag

	real(dp), allocatable :: rzbdry(:,:),rzlim(:,:)
	integer nbdry,nlim
	integer(I4B) nstep,maxstep
	character(len=50) gfile
	character(len=50) profiles
	character(len=50) profile
	character(len=100) str
	integer nl,nr,nl2,nr2
	logical alive
	integer, allocatable :: losecone(:,:)
	integer,parameter :: nsep=56
	real(sp) xsep(nsep),ysep(nsep)

	!! variables added by lye
	real(sp) :: Rmax, Rmin, Zmax, Zmin
	integer :: nrmax, nrmin, nzmax, nzmin
	integer :: npsi, ntheta, nctp
	real(sp), allocatable :: vpsi(:)
	real(sp) :: dr, dz
	!lines(:,i,j)记录第j条等高线上第i个点的r,z,deltapsi,b^2值
	!theta(:,i,j)表示第j条等高线上，第i个点的theta,r,z的值。
	real(sp), allocatable ::  lines(:,:,:), theta(:,:,:)
	real(sp) :: d1ofr, d1ofz, d12, lambda  
	integer :: ndata
	real(sp), allocatable :: cscoef(:,:), xdata(:), fdata(:), angle(:), q(:), g(:), &
							alpha(:), iota(:)
	!以下变量forcom(1:5,i,j)记录第j条线第i个点的nu,B,delta1,capitaldelta,capitalq
	!fordelta1(1:5,i,j)表示第j条线，第i个点的nu/psi,z/psi,x/psi,psi/x,psi/z
	real(sp), allocatable :: forcom(:,:,:), temp(:), tempfortheta(:), fordelta1(:,:,:)
end module vars

module miller_eq
	use inrtype, only: sp
	implicit none

	real(sp), allocatable :: br0(:), sr0(:) ! vector small r_0, big R_0
	real(sp), allocatable :: cspsisr(:,:), csbrsr(:,:) ! psi(r0), R_0(r_0)
	real(sp) :: roa, dr0, a0, R0, q0, s0, shift ! reference flux surface for millor model, and two adjecent flux distance (normalized to a0), roa = r0/a0,  dr0 = (r0-rl)/a0
	real(sp) :: RoLne, Rolte, RoLni, RoLti, RoLnc, RoLtc
	real(sp) :: B0, rhoia, teti, betae
	real(sp) :: n_u, X_u, omega_u, t_u, v_u, K_u, J_u, A_u, Phi_u, E_u, Ki_u, Q_u 
	! ki_u energy flux (watt/m^2), Q_u, energy flow (watt)
	real(sp) :: r0r, r0l, r0m, psir, psil, psim, t0e, t0i, t0c, n0e, n0i, n0c, nue, cs ! cs sound speed, nue electron collision frequency
	real(sp), allocatable :: RL(:,:), ZL(:,:), BP(:)
	real(sp) :: asf ! area of the flux surface
	real(sp) :: rhon, rhon2 ! rhon is normalized sqrt(rhot) for comparison with onetwo, use rhon to compute r0m and roa
	integer :: nprf ! number of profile points in radial
	real(sp), allocatable :: vrho(:), ne(:), Te(:), ni(:), Ti(:), nc(:), Tc(:), zeff(:)
	real(sp), allocatable :: vpsit(:), vrrho(:)
	real(sp) :: zi, zc, mimp, mcmp

!jycheng
        real(sp) :: Rpr,Rpz,Zpr,Zpz
        real(sp), allocatable :: r(:),psi(:),grid_r(:),dVr(:),Vr(:),rho(:),rho_r(:) ! 1:npsi
        real(sp), allocatable :: sf(:),ff(:) !,ne1(:),te1(:),ni1(:),ti1(:),capne(:),capte(:),capni(:),capti(:) ! 1:npsi
        real(sp), allocatable :: grid_z(:) ! 1:ntheta
        real(sp), allocatable :: knot_r(:) ! 1:npsi+kr
        real(sp), allocatable :: knot_z(:) ! 1:ntheta+kz
        real(sp), allocatable :: cscoef_z(:,:) ! 4,ntheta-1
        real(sp), allocatable :: cscoef_r(:,:),csrrho(:,:),cscoefqpsi(:,:) ! 4,npsi-1
        real(sp), allocatable :: RL_new(:,:), ZL_new(:,:),bscoef2d_RL(:,:),bscoef2d_ZL(:,:),J_inverse(:,:)

        real(sp), allocatable :: psi_nor(:)
!jycheng

	character(len=50) :: gemfilename
	namelist /miller/ dr0, nprf, rhon, zi, zc, mimp, mcmp

	real(sp) :: dum1, dum2, dum3, dpdr
end module miller_eq

!!At first, it reads a 'g' file neqdsk. Then advance guiding-center orbits.
!!This routine calculate four scals and three vectors.
!!The 4 scales are: scalB, (vecB .dot. CurlB / B^2), (GradB .dot. vecB), (GradB .dot. CurlB).
!!The 3 vectors are: vecB, (CurlB), (vecB .cross. GradB)
Program gcrz
	use inrtype
	use opes
	use splines
	use global
	use vars
	use miller_eq

	implicit none
	interface
		function pnpoly(nvert, vertx, verty, testx, testy)
		use inrtype, only: sp
		integer, intent(in) :: nvert
		real(sp), dimension(:), intent(in) :: vertx, verty
		real(sp), intent(in) :: testx,testy
		logical pnpoly
		end function pnpoly
	end interface
	integer i,j
	namelist /grid/ profiles, npsi, ntheta,nctp
	

	!gfile and particle information
	open(20, file="gcrz.in",status="old",action="read")
		read(20,nml=grid)
		!read(20,nml=particle)
		!read(20,nml=turbulence)
	close(20)
	print *, gfile
	gfile=trim(trim(profiles)//"gfile")
	print *, gfile

	!separatrix is described by polygon
	!open(20, file="separatrix.dat",status="old",action="read")
	!	do i=1,nsep
	!		read(20,*) xsep(i)
	!	end do
	!	do i=1,nsep
	!		read(20,*) ysep(i)
	!	end do
	!close(20)

	inquire(file=trim(gfile),exist=alive)
	if(.not.alive) then
		!open(50,file="gfile.err")
		!write(50,*) "gfile:" // trim(gfile) // " does not exist."
		!close(50)
		print *, "gfile:" // trim(gfile) // " does not exist."
		stop
	end if

	!inquire(file=trim(profiles),exist=alive)
	!if(.not.alive) then
	!	print *, "profiles:" // trim(profiles) // " does not exist."
	!	stop
	!end if

	open(10, file=trim(gfile),status="old",action="read")
		read(10,"(a100)") str !The last two elements should be read into nxefit, nyefit
		nr=len_trim(str)
		nl=index(str(1:nr),' ',.true.)
		nr2=len_trim(str(1:nl))
		nl2=index(str(1:nr2),' ',.true.)
		read(str(nl2:nr),*) nxefit,nyefit

		allocate(fpol(nxefit),pres(nxefit),workk1(nxefit),workk2(nxefit),fold(nxefit,nyefit),qpsi(nxefit))

		read(10,"(5e16.9)") xdim,zdim,rcentr,rgrid1,zmid	!line 2
		read(10,"(5e16.9)") rmagx,zmagx,simagx,sibdry,bcentr	!line 3
		read(10,"(5e16.9)") cpasma,simagx,xdum,rmagx,xdum	!line 4
		read(10,"(5e16.9)") zmagx,xdum,sibdry,xdum,xdum		!line 5

                print *, "sibdry=", sibdry
                print *, "simagx=", simagx
!                sibdry = -simagx
!                simagx = 0.
!                print *, "sibdry=", sibdry
!                print *, "simagx=", simagx
                
		!xdim,zdim	; Size of the domain in meters
		!rcentr,bcentr	; Reference vacuum toroidal field (m, T)
		!rgrid1		; R of left side of domain
		!zmid		; Z at the middle of the domain
		!rmagx,zmagx	; Location of magnetic axis
		!simagx		; Poloidal flux at the axis (Weber / rad)
		!sibdry		; Poloidal flux at plasma boundary (Weber / rad)

		read(10,"(5e16.9)") fpol(1:nxefit)			!line 6 to 31 Poloidal current function on uniform flux grid
		if (fpol(1) < 0.) fpol = - fpol
		read(10,"(5e16.9)") pres(1:nxefit)			!line 32 to 57 Plasma pressure in nt/m^2 on uniform flux grid
                print *, 'pres', pres(1)
                read(10,"(5e16.9)") workk1(1:nxefit)			!line 58 to 83
                print *, 'workk2', workk1(1),workk1(nxefit)                
		read(10,"(5e16.9)") workk2(1:nxefit)			!line 84 to 109
                print *, 'workk2', workk2(1),workk2(nxefit)
                read(10,"(5e16.9)") fold(1:nxefit,1:nyefit)		!line 110 to 3438 Poloidal flux in Weber/rad on grid points
                print *, 'fold', fold(1,1),fold(nxefit,nyefit)
		read(10,"(5e16.9)") qpsi(1:nxefit)			!line 3439 to 3464 q values on uniform flux grid
		!print *, qpsi
		!pause

		read(10,*) nbdry,nlim					!number of plasma boundary points and wall boundary points
		if(nbdry>0) then
			allocate(rzbdry(2,nbdry))
			read(10,"(5e16.9)") rzbdry(1:2,1:nbdry)		!line 3466 to 3502 Plasma boundary
		else
			allocate(rzbdry(2,1))
			rzbdry = 0
		end if
		if(nlim>0) then
			allocate(rzlim(2,nlim))
			read(10,"(5e16.9)") rzlim(1:2,1:nlim)		!line 3503 to 3523 Wall boundary
		else
			allocate(rzlim(2,1))
			rzlim = 0
		end if
	close(10)

	open(30,file="boundary.txt")
		write(30, "(4i5)") nbdry,nlim
		do i=1,nbdry
			write(30,"(2e16.9)") rzbdry(1:2,i)
		end do
		do i=1,nlim
			write(30,"(2e16.9)") rzlim(1:2,i)
		end do
	close(30)

	allocate(gridr(nxefit),gridz(nyefit),gridpsi(nxefit),knotr(nxefit+kr),knotz(nyefit+kz),knotpsi(nxefit+kp))
	allocate(cscoeffpol(4,nxefit-1))
	allocate(bscoefpsi(nxefit,nyefit),bscoefbphi(nxefit,nyefit),bscoefscalb(nxefit,nyefit))
	allocate(vec_b(nxefit,nyefit),vec_gradb(nxefit,nyefit),vec_curlb(nxefit,nyefit),vec_bxgradb(nxefit,nyefit))
	allocate(scal_b(nxefit,nyefit),scal2(nxefit,nyefit),scal3(nxefit,nyefit),scal4(nxefit,nyefit))

	!generate the grid points on R and Z direction
	do i=1,nxefit
		gridr(i) = rgrid1 + xdim*(i-1)/(nxefit-1)
	end do
	do j=1,nyefit
		gridz(j) = (zmid-0.5*zdim) + zdim*(j-1)/(nyefit-1)
	end do

	! set psi = 0 at magnetic axis, by lye
	do i = 1, nxefit
		do j = 1, nyefit
			fold(i,j) = fold(i,j) - simagx
			!if (fold(i,j) < 0.) print *, i, j , fold(i,j)
		end do
	end do
		print *, "after resetting fold"
	
!pause
	!generate the uniform flux grid points
	do i=1,nxefit
		!gridpsi(i) = simagx+(sibdry-simagx)*(i-1)/(nxefit-1)
		gridpsi(i) = (sibdry-simagx)*(i-1)/real(nxefit-1,sp)  ! revised by lye
	end do

	! write for matlab plot
	! write RZ grids
	open(111,file='RZ_plt.dat',status='replace', action='write')
		write(111,101) nxefit, nyefit
		write(111,103) gridr(1)
		write(111,103) gridr(nxefit)
		write(111,103) gridz(1)
		write(111,103) gridz(nyefit)
		do j = 1, nyefit
			do i = 1, nxefit
				write(111,103) gridr(i)
			end do
		end do
		do j = 1, nyefit
			do i = 1, nxefit
				write(111,103) gridz(j)
			end do
		end do
	close(111)

	100 FORMAT (I5)
	101 FORMAT (2I5)
	103 FORMAT (ES25.16)
	104 FORMAT (2ES25.16)
	105 FORMAT (3ES25.16)

	open(111,file='flux2d.dat',status='replace', action='write')
		! wirte flux
		write(111,101) nxefit, nyefit
		do j = 1, nyefit
			do i = 1, nxefit
				write(111,103) fold(i,j) 
			end do
		end do
		! write boundary
		write(111, "(4i5)") nbdry,nlim
		do i=1,nbdry
			write(111,104) rzbdry(1:2,i)
		end do
		do i=1,nlim
			write(111,104) rzlim(1:2,i)
		end do
		! Write magnetic axis
		write(111,104) rmagx,zmagx
	close(111)


        print *, "before knot sequence"
	!generate the knot sequence on R and Z direction
	call cdbbsnak(gridr,kr,knotr)
	call cdbbsnak(gridz,kz,knotz)

	!calculate coef of fpol on uniform flux grid points
	call inrcsnak(gridpsi,fpol,cscoeffpol)
	!calculate coef of psi on (R,Z)
	call cdbbscoef2d(gridr,gridz,fold,knotr,knotz,kr,kz,bscoefpsi)
        print *, "after cdbbscoef2d(gridr,gridz,fold,knotr,knotz,kr,kz,bscoefpsi)"



        
 ! RZ2flux  generate mesh of psi for output
	allocate (vpsi(npsi))
	do i=1, npsi
		vpsi(i) = (sibdry-simagx)*(i-1)/real(npsi-1,sp)  ! revised by lye
	end do

	dr = xdim/(nxefit-1)
	dz = zdim/(nyefit-1)
	Rmax = maxval(rzbdry(1,:))
	Rmin = minval(rzbdry(1,:))
	Zmax = maxval(rzbdry(2,:))
	Zmin = minval(rzbdry(2,:))
	!print *, Rmax, Rmin, Zmax, Zmin
	!pause
	nrmax = inrlocate(gridr,rmax) + 1
	nrmin = inrlocate(gridr,rmin)  
	nzmax = inrlocate(gridz,zmax) + 1
	nzmin = inrlocate(gridz,zmin) 

	allocate(lines(4,ntheta,npsi))
	allocate(theta(3,ntheta,npsi))
	allocate(forcom(5, ntheta, npsi))
        allocate(fordelta1(5,ntheta,npsi))
 
        print *, 'before contourpoints'
	call contourpoints ! revised by lye for line 2D interpolation of contour psi
        print *, 'pass contourpoints'
	!lye 
	allocate(q(npsi), g(npsi), alpha(npsi), iota(npsi))
	do i = 1, npsi
		g(i) = inrcsval(gridpsi,cscoeffpol,vpsi(i),0)
	end do
		
        do j=2,npsi
 !          write(*,*)'lines loop j=', j
		do i=1, ntheta
			!if (abs(lines(1,i,j)) <= 1.0e-3) exit
			!d1ofr = DBS2DR(1,0,lines(1,i,j),lines(2,i,j),5,3,rknot,zknot,GridNumberofR,GridNumberofZ,bscoef)
			!d1ofz = DBS2DR(0,1,lines(1,i,j),lines(2,i,j),5,3,rknot,zknot,GridNumberofR,GridNumberofZ,bscoef)
			d1ofr =	cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,lines(1,i,j),lines(2,i,j),1,0)
			d1ofz =	cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,lines(1,i,j),lines(2,i,j),0,1)
			d12 = d1ofr*d1ofr + d1ofz*d1ofz
			lines(3,i,j) = sqrt(d12)
			!if (j>35 .and. j<37) then
				!print *, lines(1,i,j), lines(2,i,j), d1ofr, d1ofz, lines(3,i,j), 'i=', i, 'j=', j
			!end if
			lines(4,i,j) = (g(j)*g(j) + d12)/(lines(1,i,j)*lines(1,i,j))
		end do
		!!!B_p/r的值。
		!write(*,*) minval(lines(3,:,j),mask=lines(3,:,j)/=0.0),maxval(lines(3,:,j))
		!write(*,*) j,lines(1,1,j)*lines(3,1,j)
	end do


  





 
 !!=========================================================================================
	! for gem local flux tube input, transform to miller equilibrium model, sr, lr, etc...
	! read input roa which is r0 over a0
	open(20, file="gcrz.in",status="old",action="read")
		read(20,nml=miller)
	close(20)
	! compute 3 flux surface for s_delta and s_kappa
	allocate(br0(npsi),sr0(npsi))
	allocate(cspsisr(4,npsi-1),csbrsr(4,npsi-1))
	allocate (cscoef(4,npsi-1))

	sr0(1) = 0.   ! small r0
	br0(1) = rmagx ! large R0

	do i = 2, npsi
		br0(i) = (lines(1,1,i) + lines(1,(ntheta+1)/2,i))/2.  !  vtheta(ntheta+1/2) = pi 
		sr0(i) = (lines(1,1,i) - lines(1,(ntheta+1)/2,i))/2.  ! r at psi grid
		!print *, sr0(i), br0(i), i
	end do
	!print *, ((ntheta+1)/2-1)*twopi/(ntheta-1), pi
	a0 = sr0(npsi)
	call inrcsnak(sr0,vpsi,cspsisr)
	call inrcsnak(sr0,br0,csbrsr)

        write(*,*)'point_1.1'


        

	allocate(r(npsi),psi(npsi),grid_r(npsi),dVr(npsi),Vr(npsi),rho(npsi),rho_r(npsi))
	allocate(sf(npsi),ff(npsi))
	allocate(grid_z(ntheta))
	allocate(knot_r(npsi+kr))
	allocate(knot_z(ntheta+kz))
	allocate(cscoef_z(4,ntheta-1))
	allocate(cscoef_r(4,npsi-1),csrrho(4,npsi-1))
	allocate(cscoefqpsi(4,nxefit-1))
	allocate(RL_new(npsi,ntheta),ZL_new(npsi,ntheta),bscoef2d_RL(npsi,ntheta),bscoef2d_ZL(npsi,ntheta),J_inverse(npsi,ntheta))
	! The RL_new, ZL_new are R(r,\theta), Z(r,\theta)
	! psi(i) is related the equally seperated radius
        !do i=1,npsi
       	!	r(i)=(i-1)*a0/(npsi-1)
	!	psi(i)=inrcsval(sr0,cspsisr,r(i),0)
	!	call ctp(psi(i),RL_new(i,:),ZL_new(i,:))
	!enddo
        psi_nor=vpsi(npsi)*psi_nor
        write(*,*)'point_1.2'
        do i=1,npsi
!                write(*,*)'call ctp loop i=',i
                r(i)=(i-1)*a0/(npsi-1)
                psi(i)=inrcsval(sr0,cspsisr,r(i),0)
                call ctp(psi(i),RL_new(i,:),ZL_new(i,:))
        enddo

	do i=1,npsi
		grid_r(i)=(i-1)*a0/(npsi-1)
	enddo	
	do i=1,ntheta
		grid_z(i)=(i-1)*pi*2.0/(ntheta-1)
	enddo
	call cdbbsnak(grid_r,kr,knot_r)
	call cdbbsnak(grid_z,kz,knot_z)
	call cdbbscoef2d(grid_r,grid_z,RL_new,knot_r,knot_z,kr,kz,bscoef2d_RL)
	call cdbbscoef2d(grid_r,grid_z,ZL_new,knot_r,knot_z,kr,kz,bscoef2d_ZL)
 	do i=1,npsi
!                write(*,*)'J_inverse loop i=',i
                do j=1,ntheta
			Rpr=cdbbsval2d(knot_r,knot_z,kr,kz,bscoef2d_RL,grid_r(i),grid_z(j),1,0)
			Rpz=cdbbsval2d(knot_r,knot_z,kr,kz,bscoef2d_RL,grid_r(i),grid_z(j),0,1)
			Zpr=cdbbsval2d(knot_r,knot_z,kr,kz,bscoef2d_ZL,grid_r(i),grid_z(j),1,0)
			Zpz=cdbbsval2d(knot_r,knot_z,kr,kz,bscoef2d_ZL,grid_r(i),grid_z(j),0,1)
			J_inverse(i,j)=abs(Rpr*Zpz-Zpr*Rpz)
		enddo
	enddo
	write(*,*)'point-0'
	do i=1,npsi
		call inrcsnak(grid_z,J_inverse(i,:),cscoef_z)
		dVr(i)=inrcsitg(grid_z,cscoef_z,grid_z(1),grid_z(ntheta))
	enddo	
	write(*,*)'point-1'
	do i=1,npsi
		call inrcsnak(grid_r,dVr,cscoef_r)
		Vr(i)=inrcsitg(grid_r,cscoef_r,grid_r(1),grid_r(i))
	enddo
	write(*,*)'point-2'
	do i=1,npsi
		rho(i)=sqrt(Vr(i)/Vr(npsi))
	enddo
	write(*,*)'point-3'
	call inrcsnak(rho,grid_r,csrrho)
	write(*,*)'ponit-4',vpsi(npsi)
  	call inrcsnak(gridpsi,qpsi,cscoefqpsi)
	do i=1,npsi
		!rho_r(i)=real((i-1))/real((npsi-1))
                ! rho --> q
                ! q --> psi
                ! psi --> r
		!r(i)=inrcsval(rho,csrrho,rho_r(i),0)
                !r(i)=grid_r(i)
                !rho_r(i)=inrcsval(grid_r,csrrho,r(i),0)
                r(i)=(i-1)*a0/(npsi-1)
                rho_r(i)=r(i)/a0
                psi(i)=inrcsval(sr0,cspsisr,r(i),0)
                ff(i)=inrcsval(gridpsi,cscoeffpol,psi(i),0)
 		sf(i)=inrcsval(gridpsi,cscoefqpsi,psi(i),0)
	enddo
        write(*,*)'a',a0,'R0',rmagx,'r_half',r(npsi/2),'q_half',sf(npsi/2),'B0',bcentr
      	
	open(113,file='profiles-new.dat',status='replace', action='write')
	write(113,"(8A14)") 'r/a', 'r', 'psi', 'ff', 'q'
	do i = 1, npsi
		write(113,204) rho_r(i),r(i),psi(i),ff(i),sf(i)
	enddo
        close(113)

        open(114,file='profiles-1d.dat',status='replace', action='write')
        !write(114,"(8A14)") 'rho_r,', 'r,', 'psi,', 'q,', 'f'
        do i = 1, npsi
                write(114,204) r(i),psi(i),sf(i), ff(i) !, ne1(i), ni1(i), te1(i), ti1(i), capni(i), capne(i), capti(i), capte(i)
        enddo
        close(114)
        
        open(115,file='profiles-2d.dat',status='replace',action='write')
        do i=1,npsi
           write(115,204)(RL_new(i,j),j=1,ntheta),(ZL_new(i,j),j=1,ntheta)
        enddo
        close(115)
        
        open(116,file='rdata.dat',status='replace',action='write') 
        do i=1,npsi
           write(116,204)(RL_new(i,j),j=1,ntheta)
        enddo
        close(116)
      
        open(117,file='zdata.dat',status='replace',action='write')
        do i=1,npsi
           write(117,204)(ZL_new(i,j),j=1,ntheta)
        enddo
        close(117)
  
        open(118,file='qpsi.dat',status='replace',action='write')
        do i=1,npsi
	   write(118,204)sf(i)
	enddo
	close(118)
	call exit()
	! jycheng	
204 FORMAT (5e16.9)

end program gcrz


subroutine ctp(psin,RL,ZL)
	use inrtype
	use global
	use vars
	use splines
	implicit none
	real(sp), intent(in) :: psin
	real(sp), dimension(:), intent(inout) :: RL(ntheta), ZL(ntheta)
	integer :: i, j, kk, nz0, nr0, m(1), n, ii
	real(sp), allocatable :: tmpr(:), tmpz(:), tmppsi(:), tmpcscoef(:,:), tmpl(:)
	real(sp), allocatable :: vr(:), vz(:), cscoef1dr(:,:,:), cscoef1dz(:,:,:)
	real(sp) :: dpsidr, dpsidz, xr, xz

	integer :: nrmagx, nzmagx
	real(sp) :: dtheta, tt, tt1, tt2, tt3, tt4
	real(sp) :: dl, rr, zz, ll

	tt1 = atan((zmax-zmagx)/(rmax-rmagx))
	tt2 = atan((zmax-zmagx)/(rmin-rmagx)) + pi
	tt3 = atan((zmin-zmagx)/(rmin-rmagx)) + pi
	tt4 = atan((zmin-zmagx)/(rmax-rmagx)) + twopi

	nrmagx = inrlocate(gridr, rmagx)
	nzmagx = inrlocate(gridz, zmagx)
	
	nr0 = nrmax - nrmin + 1
	nz0 = nzmax - nzmin + 1

	allocate (vr(nr0), vz(nz0))
	allocate (cscoef1dr(4,nr0-1,nz0), cscoef1dz(4,nz0-1,nr0))
	
	do i = 1, nr0
		vr(i) = gridr(nrmin+i-1)
	end do	
	do i = 1, nz0
		vz(i) = gridz(nzmin+i-1)
	end do

	dtheta = twopi/(ntheta-1)
	
	!print *, rmagx, zmagx
	n = nctp
	allocate(tmppsi(n),tmpcscoef(4,n-1),tmpl(n))
	
	!pause 11
	do j = 1, ntheta
		tt = (j-1)*dtheta

		if ( (tt >= 0 .and. tt< tt1) .or. (tt >= tt4 ) .or. (tt >= tt2 .and. tt< tt3) ) then 
			!print *, 'in 1'
			if (tt >= tt2 .and. tt< tt3) then
				xr = rmin
			else				
				xr = rmax
			end if

			xz = zmagx + tan(tt)*(xr-rmagx)
			dr = (xr - rmagx)/real(n-1,sp)
			dz = (xz - zmagx)/real(n-1,sp)
			dl = sqrt((xr-rmagx)**2+(xz-zmagx)**2)/real(n-1,sp)
	
			tmpl(1) = 0.
			tmppsi(1) = 0. 
		
			!print *, 'in 1.1'
			do i = 2, n
				tmpl(i) = (i-1)*dl
				rr = rmagx + (i-1)*dr
				zz = zmagx + (i-1)*dz
				tmppsi(i) = cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,rr,zz,0,0)
				if (tmppsi(i) < tmppsi(i-1)) then
					print *, 'warning for tt 1', i,j
					do kk= i, n
						tmppsi(kk) = tmppsi(i-1) + (kk-i+1)*(tmppsi(i-1)-tmppsi(i-2))
					end do
					exit
				end if
			end do
			!print *, 'in 1.2'
			call inrcsnak(tmppsi,tmpl,tmpcscoef)
			!print *, 'in 1.3'
			ll = inrcsval(tmppsi,tmpcscoef,psin,0)
                        if(ll<0)write(*,*)'in ctp', psin,j,ll
                        rr = rmagx + ll*cos(tt)
			zz = zmagx + ll*sin(tt)
					!print *, 'in 1.4'
					!print *, rr, zz
					!print *, size(RL), size(ZL)
			RL(j) = rr
					!print *, 'in 1.5'
			ZL(j) = zz
					!print *, 'in 1.6'
		end if

	!pause 12

		if ( (tt >= tt1 .and. tt< tt2) .or. (tt >= tt3 .and. tt< tt4)  ) then

			if  (tt >= tt1 .and. tt< tt2) then
				xz = zmax
			else
				xz = zmin
			end if

			xr = cotan(tt)*(xz-zmagx) + rmagx

			dr = (xr - rmagx)/real(n-1,sp)
			dz = (xz - zmagx)/real(n-1,sp)
			dl = sqrt((xr-rmagx)**2+(xz-zmagx)**2)/real(n-1,sp)
	
			tmpl(1) = 0.
			tmppsi(1) = 0. 		

			do i = 2, n
				tmpl(i) = (i-1)*dl
				rr = rmagx + (i-1)*dr
				zz = zmagx + (i-1)*dz
				tmppsi(i) = cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,rr,zz,0,0)
				if (tmppsi(i) < tmppsi(i-1)) then
					print *, 'warning for tt 1', i,j
					do kk= i, n
						tmppsi(kk) = tmppsi(i-1) + (kk-i+1)*(tmppsi(i-1)-tmppsi(i-2))
					end do
					exit
				end if
			end do
			!pause
			call inrcsnak(tmppsi,tmpl,tmpcscoef)
			
                        ll = inrcsval(tmppsi,tmpcscoef,psin,0)
                        if(ll<0)write(*,*)'in ctp', psin,j,ll
			rr = rmagx + ll*cos(tt)
			zz = zmagx + ll*sin(tt)
			RL(j) = rr
			ZL(j) = zz
		end if

	end do
	return	
end subroutine ctp

subroutine contourpoints!(psi)
	use inrtype
	use global
	use vars
        use splines
	implicit none
	!real(sp), intent(in) :: psi
	integer :: i, j, kk, nz0, nr0, m(1), n, ii
	real(sp), allocatable :: tmpr(:), tmpz(:), tmppsi(:), tmpcscoef(:,:), tmpl(:)
	real(sp), allocatable :: vr(:), vz(:), cscoef1dr(:,:,:), cscoef1dz(:,:,:)
	real(sp) :: dpsidr, dpsidz, xr, xz

	integer :: nrmagx, nzmagx
	real(sp) :: dtheta, tt, tt1, tt2, tt3, tt4
	real(sp) :: dl, rr, zz, ll

	tt1 = atan((zmax-zmagx)/(rmax-rmagx))
	tt2 = atan((zmax-zmagx)/(rmin-rmagx)) + pi
	tt3 = atan((zmin-zmagx)/(rmin-rmagx)) + pi
	tt4 = atan((zmin-zmagx)/(rmax-rmagx)) + twopi

	print *, 'four corner angles=', tt1, tt2, tt3, tt4
	!pause 1234

	nrmagx = inrlocate(gridr, rmagx)
	nzmagx = inrlocate(gridz, zmagx)
	
	nr0 = nrmax - nrmin + 1
	nz0 = nzmax - nzmin + 1

	allocate (vr(nr0), vz(nz0))
	allocate (cscoef1dr(4,nr0-1,nz0), cscoef1dz(4,nz0-1,nr0))
	
	do i = 1, nr0
		vr(i) = gridr(nrmin+i-1)
	end do	
	do i = 1, nz0
		vz(i) = gridz(nzmin+i-1)
	end do

	dtheta = twopi/(ntheta-1)
	
	!print *, rmagx, zmagx
	n = nctp !npsi*2
	allocate(tmppsi(n),tmpcscoef(4,n-1),tmpl(n))

	do j = 1, ntheta
		lines(1,j,1) = rmagx
		lines(2,j,1) = zmagx
	end do
		
	do j = 1, ntheta
		tt = (j-1)*dtheta

		if ( (tt >= 0 .and. tt< tt1) .or. (tt >= tt4 ) .or. (tt >= tt2 .and. tt< tt3) ) then 
			
			if (tt >= tt2 .and. tt< tt3) then
				xr = rmin
			else				
				xr = rmax
			end if

			xz = zmagx + tan(tt)*(xr-rmagx)
			dr = (xr - rmagx)/real(n-1,sp)
			dz = (xz - zmagx)/real(n-1,sp)
			dl = sqrt((xr-rmagx)**2+(xz-zmagx)**2)/real(n-1,sp)
	
			tmpl(1) = 0.
			tmppsi(1) = 0. 
		

			do i = 2, n
				tmpl(i) = (i-1)*dl
				rr = rmagx + (i-1)*dr
				zz = zmagx + (i-1)*dz
				tmppsi(i) = cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,rr,zz,0,0)
				if (tmppsi(i) < tmppsi(i-1)) then
                                  write(*,200)i,j,(rr-rmagx)/zdim,zz/zdim,tmppsi(i)

					do kk= i, n
						tmppsi(kk) = tmppsi(i-1) + (kk-i+1)*(tmppsi(i-1)-tmppsi(i-2))
					end do
					exit
				end if
			end do
			!pause
			call inrcsnak(tmppsi,tmpl,tmpcscoef)
200 format(1x,i4,1x,i4,3(2x,e10.3))			
			do i = 2, npsi
				ll = inrcsval(tmppsi,tmpcscoef,vpsi(i),0)
				rr = rmagx + ll*cos(tt)
				zz = zmagx + ll*sin(tt)
				lines(1,j,i) = rr
				lines(2,j,i) = zz
			end do		
		end if

		if ( (tt >= tt1 .and. tt< tt2) .or. (tt >= tt3 .and. tt< tt4)  ) then

			if  (tt >= tt1 .and. tt< tt2) then
				xz = zmax
			else
				xz = zmin
			end if

			xr = cotan(tt)*(xz-zmagx) + rmagx

			dr = (xr - rmagx)/real(n-1,sp)
			dz = (xz - zmagx)/real(n-1,sp)
			dl = sqrt((xr-rmagx)**2+(xz-zmagx)**2)/real(n-1,sp)
	
			tmpl(1) = 0.
			tmppsi(1) = 0. 		

			do i = 2, n
				tmpl(i) = (i-1)*dl
				rr = rmagx + (i-1)*dr
				zz = zmagx + (i-1)*dz
				tmppsi(i) = cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,rr,zz,0,0)
				if (tmppsi(i) < tmppsi(i-1)) then
					print *, 'warning for tt 2', i,j,rr,zz,tmppsi(i)
					do kk= i, n
						tmppsi(kk) = tmppsi(i-1) + (kk-i+1)*(tmppsi(i-1)-tmppsi(i-2))
					end do
					exit
				end if
			end do
			!pause
			call inrcsnak(tmppsi,tmpl,tmpcscoef)
			
			do i = 2, npsi
				ll = inrcsval(tmppsi,tmpcscoef,vpsi(i),0)
				rr = rmagx + ll*cos(tt)
				zz = zmagx + ll*sin(tt)
				lines(1,j,i) = rr
				lines(2,j,i) = zz
			end do		
		end if

	end do

		

	open(111,file='contourpsi.dat',status='replace', action='write')
		! wirte flux
		write(111,101) npsi, ntheta
		do i = 1, npsi
			write(111,103) vpsi(i)
		end do
		do i = 1, npsi
			do j = 1, ntheta
				write(111,104) lines(1:2,j,i) 
			end do
		end do
	close(111)

	100 FORMAT (I5)
	101 FORMAT (2I5)
	103 FORMAT (ES25.16)
	104 FORMAT (2ES25.16)
	!pause 11
	 return	
end subroutine contourpoints

!http://www.ecse.rpi.edu/~wrf/Research/Short_Notes/pnpoly.html
!Argument	Meaning
!nvert 	Number of vertices in the polygon. Whether to repeat the first vertex at the end is discussed below.
!vertx, verty 	Arrays containing the x- and y-coordinates of the polygon's vertices.
!testx, testy	X- and y-coordinate of the test point. 
function pnpoly(nvert, vertx, verty, testx, testy)
	use inrtype, only: sp
	implicit none
	integer, intent(in) :: nvert
	real(sp), dimension(:), intent(in) :: vertx, verty
	real(sp), intent(in) :: testx,testy
	logical pnpoly
	integer i,j

	pnpoly=.false.
	j=nvert
	do i=1,nvert
		if ( ( (verty(i)>testy) .neqv. (verty(j)>testy) ) .and. &
			(testx < (vertx(j)-vertx(i)) * (testy-verty(i)) / (verty(j)-verty(i)) + vertx(i)) ) then
			pnpoly = .not. pnpoly
		end if
		j = i
	end do
end function pnpoly
