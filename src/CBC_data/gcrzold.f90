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
	integer :: npsi, ntheta
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
        integer :: psi_num_density,psi_num_temp
        real(sp) :: Rpr,Rpz,Zpr,Zpz
        real(sp), allocatable :: r(:),psi(:),grid_r(:),dVr(:),Vr(:),rho(:),rho_r(:) ! 1:npsi
        real(sp), allocatable :: sf(:),ff(:),ne1(:),te1(:),ni1(:),ti1(:),capne(:),capte(:),capni(:),capti(:) ! 1:npsi
        real(sp), allocatable :: grid_z(:) ! 1:ntheta
        real(sp), allocatable :: knot_r(:) ! 1:npsi+kr
        real(sp), allocatable :: knot_z(:) ! 1:ntheta+kz
        real(sp), allocatable :: cscoef_z(:,:) ! 4,ntheta-1
        real(sp), allocatable :: cscoef_r(:,:),csrrho(:,:),cscoefqpsi(:,:) ! 4,npsi-1
        real(sp), allocatable :: RL_new(:,:), ZL_new(:,:),bscoef2d_RL(:,:),bscoef2d_ZL(:,:),J_inverse(:,:)

        real(sp), allocatable :: psi_nor(:),density_new(:),psi_nor_temp(:),temp_new(:),density_new1(:),temp_new1(:)
        real(sp), allocatable :: cspsitemp(:,:),cspsidensity(:,:),csrdensity(:,:),csrtemp(:,:)
        real(sp), allocatable :: density_new_0,temp_new_0
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
	namelist /grid/ profiles, npsi, ntheta
	

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
                print *, "sibdry=", sibdry
                print *, "simagx=", simagx
                
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
        ! jycheng: add the density and temperature profile
        open(40,file='den_gene_case5.txt',status="old",action="read")
        read(40,"(4i5)")psi_num_density
        write(*,*)psi_num_density
        allocate(psi_nor(1:psi_num_density),density_new(1:psi_num_density))
        do i=1,psi_num_density
           read(40,*)psi_nor(i),density_new(i)
        enddo
        write(*,*)psi_nor(1),density_new(1)
        close(40)
        open(50,file='temp_gene_case5_fix.txt',status="old",action="read")
        read(50,"(4i5)")psi_num_temp
        allocate(psi_nor_temp(1:psi_num_temp),temp_new(1:psi_num_temp))
        do i=1,psi_num_temp
           read(50,*)psi_nor_temp(i),temp_new(i)
        enddo
        close(50)
        write(*,*)'after reading'
         !jycheng

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

	!normalise
	!rmag = 2.291807	!?
	!bmag = 2.254007	!?
	!rmag = rgrid1 + 0.5*xdim
	!bmag = abs(bcentr)
	!vmag = sqrt(2.0*E_charge*Ti/Mp)

	!rmag_rho = Mp*vmag/E_charge/bmag
	!epsilon_l=rmag_rho/rmag
	!tmag = rmag/vmag

	!gridr = gridr / rmag
	!gridz = gridz / rmag
	!gridpsi = gridpsi / bmag /rmag /rmag
	!fpol = fpol / bmag /rmag
	!fold = fold / bmag /rmag /rmag
	!R_limiter=R_limiter/rmag
	!Z_limiter=Z_limiter/rmag
	!!dl_tb = length_tb/rmag
	!dt_tb = time_tb/tmag
	!xsep = xsep / rmag
	!ysep = ysep / rmag

        print *, "before knot sequence"
	!generate the knot sequence on R and Z direction
	call cdbbsnak(gridr,kr,knotr)
	call cdbbsnak(gridz,kz,knotz)

	!calculate coef of fpol on uniform flux grid points
	call inrcsnak(gridpsi,fpol,cscoeffpol)
	!calculate coef of psi on (R,Z)
	call cdbbscoef2d(gridr,gridz,fold,knotr,knotz,kr,kz,bscoefpsi)
        print *, "after cdbbscoef2d(gridr,gridz,fold,knotr,knotz,kr,kz,bscoefpsi)"






        goto 998
	do i=1,nxefit
		do j=1,nyefit
			vec_b(i,j)%r = -cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,gridr(i),gridz(j),0,1)/gridr(i)
			vec_b(i,j)%z =  cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,gridr(i),gridz(j),1,0)/gridr(i)
			!vec_b(i,j)%r = -vec_b(i,j)%r !!!change sign of Bp   commentted by lye
			!vec_b(i,j)%z = -vec_b(i,j)%z !!!change sign of Bp
			tmpfpol = inrcsval(gridpsi,cscoeffpol,fold(i,j),0)
			vec_b(i,j)%p = tmpfpol/gridr(i)
			scal_b(i,j) = sqrt(vec_b(i,j) .dot. vec_b(i,j))	!scal_b stands for scale magnetic field on (R,Z)
		end do
	end do

	!calculate coef of bphi on (R,Z)
	call cdbbscoef2d(gridr,gridz,vec_b(:,:)%p,knotr,knotz,kr,kz,bscoefbphi)
	!calculate coef of scal b on (R,Z)
	call cdbbscoef2d(gridr,gridz,scal_b(:,:),knotr,knotz,kr,kz,bscoefscalb)

	do i=1,nxefit
		do j=1,nyefit
			vec_gradb(i,j)%r = cdbbsval2d(knotr,knotz,kr,kz,bscoefscalb,gridr(i),gridz(j),1,0)
			vec_gradb(i,j)%z = cdbbsval2d(knotr,knotz,kr,kz,bscoefscalb,gridr(i),gridz(j),0,1)
			vec_gradb(i,j)%p = 0.0
			vec_bxgradb(i,j) = vec_b(i,j) .cross. vec_gradb(i,j)
		end do
	end do

	do i=1,nxefit
		do j=1,nyefit
         	vec_curlb(i,j)%r = -cdbbsval2d(knotr,knotz,kr,kz,bscoefbphi,gridr(i),gridz(j),0,1)
			tmp1 = cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,gridr(i),gridz(j),2,0)
			tmp2 = cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,gridr(i),gridz(j),0,2)
			vec_curlb(i,j)%z = -(tmp1+tmp2)/gridr(i)
			vec_curlb(i,j)%z = -vec_curlb(i,j)%z !!!change sign of Bp
			vec_curlb(i,j)%p = vec_b(i,j)%p/gridr(i) +  cdbbsval2d(knotr,knotz,kr,kz,bscoefbphi,gridr(i),gridz(j),1,0)
		end do
	end do

	do i=1,nxefit
		do j=1,nyefit
			!(vecB .dot. CurlB / B^2)  
			scal2(i,j) = (vec_b(i,j) .dot. vec_curlb(i,j)) / (scal_b(i,j) * scal_b(i,j))

			!(GradB .dot. vecB)
			scal3(i,j) = vec_gradb(i,j) .dot. vec_b(i,j)

			!(GradB .dot. CurlB)
			scal4(i,j) = vec_gradb(i,j) .dot. vec_curlb(i,j)
		end do
	end do

	allocate(bscoefscal1(nxefit,nyefit),bscoefscal2(nxefit,nyefit),bscoefscal3(nxefit,nyefit),bscoefscal4(nxefit,nyefit))

	!coef of scalb
	call cdbbscoef2d(gridr,gridz,scal_b,knotr,knotz,kr,kz,bscoefscal1)
	!coef of scal2
	call cdbbscoef2d(gridr,gridz,scal2,knotr,knotz,kr,kz,bscoefscal2)
	!coef of scal3
	call cdbbscoef2d(gridr,gridz,scal3,knotr,knotz,kr,kz,bscoefscal3)
	!coef of scal4
	call cdbbscoef2d(gridr,gridz,scal4,knotr,knotz,kr,kz,bscoefscal4)

	allocate(bscoefvec1r(nxefit,nyefit),bscoefvec1z(nxefit,nyefit),bscoefvec1p(nxefit,nyefit))
	allocate(bscoefvec2r(nxefit,nyefit),bscoefvec2z(nxefit,nyefit),bscoefvec2p(nxefit,nyefit))
	allocate(bscoefvec3r(nxefit,nyefit),bscoefvec3z(nxefit,nyefit),bscoefvec3p(nxefit,nyefit))
	!coef of vec_1
	call cdbbscoef2d(gridr,gridz,vec_b(:,:)%r,knotr,knotz,kr,kz,bscoefvec1r)
	call cdbbscoef2d(gridr,gridz,vec_b(:,:)%z,knotr,knotz,kr,kz,bscoefvec1z)
	call cdbbscoef2d(gridr,gridz,vec_b(:,:)%p,knotr,knotz,kr,kz,bscoefvec1p)
	!coef of vec_2
	call cdbbscoef2d(gridr,gridz,vec_curlb(:,:)%r,knotr,knotz,kr,kz,bscoefvec2r)
	call cdbbscoef2d(gridr,gridz,vec_curlb(:,:)%z,knotr,knotz,kr,kz,bscoefvec2z)
	call cdbbscoef2d(gridr,gridz,vec_curlb(:,:)%r,knotr,knotz,kr,kz,bscoefvec2p)
	!coef of vec_3
	call cdbbscoef2d(gridr,gridz,vec_bxgradb(:,:)%r,knotr,knotz,kr,kz,bscoefvec3r)
	call cdbbscoef2d(gridr,gridz,vec_bxgradb(:,:)%z,knotr,knotz,kr,kz,bscoefvec3z)
	call cdbbscoef2d(gridr,gridz,vec_bxgradb(:,:)%r,knotr,knotz,kr,kz,bscoefvec3p)

        print *, 'after coef of vec_3'











998     continue
        
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
           write(*,*)'lines loop j=', j
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

!pause!!!  write(*,*) count(lines(1,:,1)/=0.0),count(lines(1,:,2)/=0.0),count(lines(1,:,3)/=0.0),&
!!!             count(lines(1,:,4)/=0.0),count(lines(1,:,5)/=0.0),count(lines(1,:,10)/=0.0),&
!!!			 count(lines(1,:,20)/=0.0),count(lines(1,:,30)/=0.0),count(lines(1,:,32)/=0.0)
  

	
        goto 999



 !pause
  do j=2,nPsi
     !ndata = count(lines(1,:,j)/=0.0)
!print *, 'j is ', j
	  ndata = ntheta
!!!	 write(*,*) j,lines(:,1,j)-lines(:,ndata,j)
	 allocate(xdata(ndata))
	 allocate(fdata(ndata))
	 allocate(angle(ndata))
	 allocate(temp(ndata))
	 !allocate(break(ndata))
	 allocate(cscoef(4,ndata-1))


     
	 xdata(1) = 0.0
	 do i=2,ndata
	    xdata(i) = xdata(i-1) + sqrt((lines(1,i,j)-lines(1,i-1,j))**2 + &
	 	           (lines(2,i,j)-lines(2,i-1,j))**2)
	 end do
     !计算q
	 do i=1,ndata
	    fdata(i) = 1.0/(lines(1,i,j)*lines(3,i,j))
		!if (j>35 .and. j<37) then
		!print *, fdata(i),lines(1,i,j),lines(3,i,j), 'i=', i
		!end if
		!if (fdata(i) > 100 ) fdata(i) = fdata(i-1)
		!print *, fdata(i)
	 end do
!!!	 write(*,*) j,fdata(1),fdata(ndata)
	 !call DCSPER(ndata,xdata,fdata,break,cscoef)
	call inrcsp(xdata,fdata,cscoef)
	
    !call CSPER(xdata,fdata,break,cscoef)
	 !q(j) = 0.5*g(2,j)*DCSITG(xdata(1),xdata(ndata),ndata-1,break,cscoef)/pi
	!q(j) = 0.5*g(2,j)*CSITG(xdata(1),xdata(ndata),break,cscoef)/pi
	q(j) = 0.5*g(j)*inrcsitg(xdata,cscoef,xdata(1),xdata(ndata))/pi
	!if (j==npsi) q(npsi) = qpsi(nxefit) ! for LFS q is not accurate, use experimental data
	!if (j>35 .and. j<37) then
	!print *, q(j), g(j), j,'================'
	!print *,  fdata
	!end if
	!pause
	!print *, 'start alpha'
     !计算alpha
         do i=1,ndata
	    fdata(i) = lines(4,i,j)*lines(1,i,j)/lines(3,i,j)
	 end do
!!!	 write(*,*) j,fdata(1),fdata(ndata)
	 !call DCSPER(ndata,xdata,fdata,break,cscoef)
	 !call CSPER(xdata,fdata,break,cscoef)
	call inrcsp(xdata,fdata,cscoef)
	 !alpha(j) = 0.5*DCSITG(xdata(1),xdata(ndata),ndata-1,break,cscoef)/pi
     !alpha(j) = 0.5*CSITG(xdata(1),xdata(ndata),break,cscoef)/pi
	alpha(j) = 0.5*inrcsitg(xdata,cscoef,xdata(1),xdata(ndata))/pi
	!print *, q(j), g(j), alpha(j), j,'================'
	 !计算iota
	 iota(j) = alpha(j) - g(j)*q(j)
	!print *, 'start theta'
	 !计算theta
         do i=1,ndata
	    fdata(i) = (lines(1,i,j)*lines(4,i,j))/lines(3,i,j)
	 end do
!!!	 write(*,*) j,fdata(1),fdata(ndata)
	 !call DCSPER(ndata,xdata,fdata,break,cscoef)
	 !call CSPER(xdata,fdata,break,cscoef)
	call inrcsp(xdata,fdata,cscoef)

	 !lambda = 2.0*pi*alpha(j)/DCSITG(xdata(1),xdata(ndata),ndata-1,break,cscoef)
      !lambda = 2.0*pi*alpha(j)/CSITG(xdata(1),xdata(ndata),break,cscoef)
	lambda = 2.0*pi*alpha(j)/inrcsitg(xdata,cscoef,xdata(1),xdata(ndata))
	 !误差修正因子lambda。误差越小，lambda越接近于1。
!!!	 write(*,*) j,lambda
        do i=1,ndata
	    !angle(i) = lambda*DCSITG(xdata(1),xdata(i),ndata-1,break,cscoef)/alpha(j)
        !angle(i) = lambda*CSITG(xdata(1),xdata(i),break,cscoef)/alpha(j)
		angle(i) = lambda*inrcsitg(xdata,cscoef,xdata(1),xdata(i))/alpha(j)
	 end do

	 do i=1,ntheta
	    theta(1,i,j) = 2.0*pi*(i-1)/(ntheta-1)
	 end do
     !计算theta对应的r值
!!!	 write(*,*) j,lines(1,1,j),lines(1,ndata,j),lines(2,1,j),lines(2,ndata,j)
	 !call DCSPER(ndata,angle,lines(1,1:ndata,j),break,cscoef)
	! call CSPER(angle,lines(1,1:ndata,j),break,cscoef)
	call inrcsp(angle, lines(1,1:ndata,j), cscoef)
	 do i=1,ntheta
	    !theta(2,i,j) = DCSVAL(theta(1,i,j),ndata-1,break,cscoef)
            !theta(2,i,j) = CSVAL(theta(1,i,j),break,cscoef)
			theta(2,i,j) = inrcsval(angle,cscoef,theta(1,i,j),0)
	 end do
	 !计算theta对应的z值
	 !call DCSPER(ndata,angle,lines(2,1:ndata,j),break,cscoef)
	 !call CSPER(angle,lines(2,1:ndata,j),break,cscoef)
	call inrcsp(angle, lines(2,1:ndata,j), cscoef)
	 do i=1,ntheta
	    !theta(3,i,j) = DCSVAL(theta(1,i,j),ndata-1,break,cscoef)
            !theta(3,i,j) = CSVAL(theta(1,i,j),break,cscoef)
			theta(3,i,j) = inrcsval(angle,cscoef,theta(1,i,j),0)
	 end do
     !位置修正
     !do i=1, ntheta
     !   if (theta(2,i,1) < R_i .and. theta(2,i,1) > (R_i-1.0e-5)) theta(2,i,1) = R_i
	 !   if (theta(2,i,1) > R_o .and. theta(2,i,1) < (R_o+1.0e-5)) theta(2,i,1) = R_o
	 !   if (theta(3,i,1) < -Z_t .and. theta(3,i,1) > (-Z_t-1.0e-5)) theta(3,i,1) = -Z_t
	 !   if (theta(3,i,1) > Z_t .and. theta(3,i,1) < (Z_t+1.0e-5)) theta(3,i,1) = Z_t
     !end do

	 !计算forcom(1:5,i,j)。第j条线第i个点的nu,B,delta1,capitaldelta,capitalq
     !计算forcom(1,i,j),nu
	 do i=1,ndata
	    fdata(i) = 1.0/(lines(4,i,j)*lines(1,i,j)**2)
	 end do
	 !call DCSPER(ndata,angle,fdata,break,cscoef)
   	 !call CSPER(angle,fdata,break,cscoef)
	call inrcsp(angle,fdata,cscoef)	 

	 do i=1,ndata
	    !temp(i) = g(2,j)*alpha(j)*DCSITG(angle(1),angle(i),ndata-1,break,cscoef)-q(j)*angle(i)
           temp(i) = g(j)*alpha(j)*inrcsitg(angle,cscoef,angle(1),angle(i))-q(j)*angle(i)
         end do
	 !write(*,'(3ES12.3,I3)'), g(j), q(j), alpha(j), j
	 !write(*,*) '======================'
!!!	 write(*,*) j,temp(1),temp(ndata)
     if (abs(temp(1)-temp(ndata))<=1.0e-2) then
	    temp(ndata) = temp(1)
	 else	    
	   write(*,*) "there is something error : nu(0) /= nu(2*pi)."
		!write(*,*) temp(1),temp(ndata),ndata,"n"
	 end if
	 !call DCSPER(ndata,angle,temp,break,cscoef)
  	 !call CSPER(angle,temp,break,cscoef)
	call inrcsp(angle,temp,cscoef)
	 do i=1, ntheta
	    !forcom(1,i,j) = DCSVAL(theta(1,i,j),ndata-1,break,cscoef)
            !forcom(1,i,j) = CSVAL(theta(1,i,j),break,cscoef)
			forcom(1,i,j) = inrcsval(angle,cscoef,theta(1,i,j),0)
	 end do
	 if (abs(forcom(1,1,j)-forcom(1, ntheta,j))<=1.0e-2) then
	    forcom(1,ntheta,j) = forcom(1,1,j)
	 else
	    write(*,*) "there is something error : nu(0) /= nu(2*pi)."
		!write(*,*) forcom(1,1,j),forcom(1, ntheta,j),j,"j"
	 end if
!!!	 write(*,*) j,forcom(1,1,j),forcom(1,33,j)

	 !计算forcom(2,i,j),B
	 !call DCSPER(ndata,angle,lines(4,1:ndata,j),break,cscoef)
!	 call CSPER(angle,lines(4,1:ndata,j),break,cscoef)
	 call inrcsp(angle,lines(4,1:ndata,j),cscoef)
	 do i=1,ntheta
	    !forcom(2,i,j) = sqrt(DCSVAL(theta(1,i,j),ndata-1,break,cscoef))
            forcom(2,i,j) = sqrt(inrcsval(angle,cscoef,theta(1,i,j),0))
	 end do
	 if (abs(forcom(2,1,j)-forcom(2,ntheta,j))<=1.0e-5) then
	    forcom(2,ntheta,j) = forcom(2,1,j)
	 else
	    write(*,*) "there is something error : B(0) /= B(2*pi)."
	 end if
!!!	 write(*,*) j,forcom(2,1,j),forcom(2,33,j)

	 deallocate(xdata)
	 deallocate(fdata)
	 deallocate(angle)
	 deallocate(temp)
!	 deallocate(break)
	 deallocate(cscoef)
  end do


	 !计算nu/psi
	allocate (cscoef(4,npsi-2))
	do i = 1,ntheta-1
	    !call DCSINT(GridNumberofPsi-1,g(1,1:GridNumberofPsi-1),forcom(1,i,1:GridNumberofPsi-1),break2,cscoef2)
	    !call CSINT(g(1,1:GridNumberofPsi-1),forcom(1,i,1:GridNumberofPsi-1),break2,cscoef2)
		call inrcsnak(Vpsi(2:npsi),forcom(1,i,2:npsi),cscoef)
	 	do j=2,nPsi
	 	   !fordelta1(1,i,j) = DCSDER(1,g(1,j),GridNumberofPsi-2,break2,cscoef2)
           !fordelta1(1,i,j) = D_CSDER(1,g(1,j),break2,cscoef2)
			fordelta1(1,i,j) = inrcsval(Vpsi(2:npsi),cscoef,Vpsi(j),1)
                end do
	end do
	fordelta1(1,ntheta,2:nPsi) = fordelta1(1,1,2:npsi)
	 
	 !计算z/psi
	do i=1, ntheta-1
	    !call DCSINT(GridNumberofPsi-1,g(1,1:GridNumberofPsi-1),theta(3,i,1:GridNumberofPsi-1),break2,cscoef2)
	    !call CSINT(g(1,1:GridNumberofPsi-1),theta(3,i,1:GridNumberofPsi-1),break2,cscoef2)
		call inrcsnak(Vpsi(2:npsi), theta(3,i,2:npsi), cscoef)	
		do j=2, npsi
		   !fordelta1(2,i,j) = DCSDER(1,g(1,j),GridNumberofPsi-2,break2,cscoef2)
           !fordelta1(2,i,j) = D_CSDER(1,g(1,j),break2,cscoef2)
			fordelta1(2,i,j) = inrcsval(Vpsi(2:npsi), cscoef, Vpsi(j),1)
		end do
	end do
        fordelta1(2,ntheta,2:npsi) = fordelta1(2,1,2:npsi)
	
	 !计算r/psi
	 do i=1, ntheta-1
	    !call DCSINT(GridNumberofPsi-1,g(1,1:GridNumberofPsi-1),theta(2,i,1:GridNumberofPsi-1),break2,cscoef2)
	    !call CSINT(g(1,2:npsi),theta(2,i,2:npsi),break2,cscoef2)
		call inrcsnak(Vpsi(2:npsi), theta(2,i,2:npsi), cscoef)	
		do j=2, npsi
		   !fordelta1(3,i,j) = DCSDER(1,g(1,j),GridNumberofPsi-2,break2,cscoef2)
           !fordelta1(3,i,j) = D_CSDER(1,g(1,j),break2,cscoef)
			fordelta1(3,i,j) = inrcsval(Vpsi(2:npsi), cscoef, Vpsi(j),1)
		end do
	 end do
         fordelta1(3, ntheta, 2:npsi) = fordelta1(3,1,2:npsi)

	 !计算psi/x,psi/z
	 do j= 2, nPsi
	    do i= 1, ntheta-1 
		   !fordelta1(4,i,j) = DBS2DR(1,0,theta(2,i,j),theta(3,i,j),5,3,rknot,zknot,GridNumberofR,GridNumberofZ,bscoef)
		   !fordelta1(5,i,j) = DBS2DR(0,1,theta(2,i,j),theta(3,i,j),5,3,rknot,zknot,GridNumberofR,GridNumberofZ,bscoef)
		   !fordelta1(4,i,j) = D_BS2DR(1,0,theta(2,i,j),theta(3,i,j),5,3,rknot,zknot,GridNumberofR,GridNumberofZ,bscoef)
		   !fordelta1(5,i,j) = D_BS2DR(0,1,theta(2,i,j),theta(3,i,j),5,3,rknot,zknot,GridNumberofR,GridNumberofZ,bscoef)
			fordelta1(4,i,j) =	cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,theta(2,i,j), theta(3,i,j),1,0)
			fordelta1(5,i,j) =	cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,theta(2,i,j), theta(3,i,j),0,1)
	    end do
	 end do
	 fordelta1(4,ntheta,2:nPsi) = fordelta1(4,1,2:npsi)
	 fordelta1(5,ntheta,2:nPsi) = fordelta1(5,1,2:npsi)

	 !计算delta1
	 do j=2, npsi
	    do i=1, ntheta
	       forcom(3,i,j) = g(j)*fordelta1(1,i,j)+ &
	 	   (fordelta1(2,i,j)*fordelta1(4,i,j)-fordelta1(3,i,j)*fordelta1(5,i,j)) &
	 	   /theta(2,i,j)
	    end do
	 end do

	! at magnetic axis
	q(1) = qpsi(1)
	g(1) = g(1)
	iota(1) = 0.
	alpha(1) = q(1)*g(1)
	do j = 1, ntheta
		theta(2,j,1) = rmagx
		theta(3,j,1) = zmagx
		forcom(1,j,1) = 0. ! nu
		forcom(2,j,1) = g(1)/rmagx ! B
		forcom(3,j,1) = 0.
		forcom(4,j,1) = 0.
		forcom(5,j,1) = 0.
	end do
	 deallocate(cscoef)

	!do j = 1, ntheta
	!print *, fordelta1(:,j,npsi), 'j=', j\

	!end do
	!print *, forcom(3,:,:)

	open(111,file='contourtt.dat',status='replace', action='write')
		! wirte flux
		write(111,101) npsi, ntheta
		do i = 1, npsi
			do j = 1, ntheta
				write(111,104) theta(2:3,j,i) 
			end do
		end do
	close(111)
















999 continue
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

	! jycheng
        allocate(cspsitemp(4,psi_num_temp-1),cspsidensity(4,psi_num_density-1))
        allocate(density_new1(1:npsi),temp_new1(1:npsi))
        allocate(csrdensity(4,npsi-1),csrtemp(4,npsi-1))
        write(*,*)'point_1.1'
        !psi_nor_temp=vpsi(npsi)*psi_nor_temp
        !psi_nor=vpsi(npsi)*psi_nor
        !call inrcsnak(psi_nor_temp,temp_new,cspsitemp) 
        !call inrcsnak(psi_nor,density_new,cspsidensity)
        !write(*,*)'point_1.2'
        !do i=1,npsi
        !        r(i)=(i-1)*a0/(npsi-1)
        !        psi(i)=inrcsval(sr0,cspsisr,r(i),0)
        !        call ctp(psi(i),RL_new(i,:),ZL_new(i,:))
        !        density_new1(i)=inrcsval(psi_nor,cspsidensity,psi(i),0)
        !        temp_new1(i)=inrcsval(psi_nor_temp,cspsitemp,psi(i),0)
        !enddo
        
        !call inrcsnak(r,density_new1,csrdensity)
        !call inrcsnak(r,temp_new1,csrtemp)
        

	allocate(r(npsi),psi(npsi),grid_r(npsi),dVr(npsi),Vr(npsi),rho(npsi),rho_r(npsi))
	allocate(sf(npsi),ff(npsi),ne1(npsi),te1(npsi),ni1(npsi),ti1(npsi),capne(npsi),capte(npsi),capni(npsi),capti(npsi))
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
        psi_nor_temp=vpsi(npsi)*psi_nor_temp
        psi_nor=vpsi(npsi)*psi_nor
        density_new_0=density_new(1)
        temp_new_0=temp_new(1)
        write(*,*)'density in axis',density_new(1),'temp in axis',temp_new(1)
        density_new=density_new/density_new(1)
        temp_new=temp_new/temp_new(1)
        call inrcsnak(psi_nor_temp,temp_new,cspsitemp) 
        call inrcsnak(psi_nor,density_new,cspsidensity)
        write(*,*)'point_1.2'
        do i=1,npsi
                write(*,*)'call ctp loop i=',i
                r(i)=(i-1)*a0/(npsi-1)
                psi(i)=inrcsval(sr0,cspsisr,r(i),0)
                call ctp(psi(i),RL_new(i,:),ZL_new(i,:))
                density_new1(i)=inrcsval(psi_nor,cspsidensity,psi(i),0)
                temp_new1(i)=inrcsval(psi_nor_temp,cspsitemp,psi(i),0)
        enddo

        call inrcsnak(r,density_new1,csrdensity)
        call inrcsnak(r,temp_new1,csrtemp)

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
                write(*,*)'J_inverse loop i=',i
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
                !ff(i)=psi(i)/gridpsi(nxefit)
                ne1(i)=density_new1(i)
                capne(i)=-inrcsval(r,csrdensity,r(i),1)/ne1(i)*rmagx
                te1(i)=temp_new1(i)
                capte(i)=-inrcsval(r,csrtemp,r(i),1)/te1(i)*rmagx
                !call ctp(psi(i),RL_new(i,:),ZL_new(i,:))
                !ne1(i)=inrcsval(psi_nor,cspsidensity,psi(i),0)
		!ne1(i)=(cosh((rho_r(i)-0.5+0.35)/0.02)/cosh((rho_r(i)-0.5-0.35)/0.02))**(-2.22*0.35*0.02/2.)
		!capne(i)=2.22/2.*(tanh((rho_r(i)-0.5+0.35)/0.02)-tanh((rho_r(i)-0.5-0.35)/0.02))
		!te1(i)=(cosh((rho_r(i)-0.5+0.35)/0.02)/cosh((rho_r(i)-0.5-0.35)/0.02))**(-6.91*0.35*0.02/2.)
		!capte(i)=6.91/2.*(tanh((rho_r(i)-0.5+0.35)/0.02)-tanh((rho_r(i)-0.5-0.35)/0.02))
		ni1(i)=ne1(i)
		capni(i)=capne(i)
		ti1(i)=te1(i)
		capti(i)=capte(i)
	enddo
        write(*,*)'a',a0,'R0',rmagx,'r_half',r(npsi/2),'q_half',sf(npsi/2),'n0',density_new_0,'t0',temp_new_0,'B0',bcentr
      	
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
















 ! read equilibrium profiles and compute related parameters
	allocate (vrho(nprf), ne(nprf), Te(nprf), ni(nprf), Ti(nprf), nc(nprf), Tc(nprf), Zeff(nprf) )
	allocate (vpsit(npsi))
	allocate (vrrho(nprf))
	!open(20, file=trim(profiles),status="old",action="read")
	!	read(20,*)
	!	do i = 1, nprf
	!		read(20, *) vrho(i), ne(i), Te(i)
	!		ne(i) = ne(i)*1.d6 ! /cm^3 to /m^3
			!print *, vrho(i), ne(i), Te(i)
	!	end do
	!close(20)
		! vrho is sqrt(psi_t), transform to psi_p and r0
!write(unit=filenamedata,fmt="('../QN_DATA/QNA/','QNA',I2.2,'.bin')") fd_id
	!write(unit=profiles,fmt="(trim(profiles)'/rho')") 

	! load profiles
	! read vrho  ! vrho is sqrt(psi_t), need to transform to psi_p and r0
	!profiles = trim(profiles)
!	profile=trim(trim(profiles)//"rho")
!	open(20, file=trim(profile),status="old",action="read")
!		do i = 1, nprf
!			read(20, *) vrho(i)
!			!print *, vrho(i)
!		end do
!	close(20)
!
!	profile=trim(trim(profiles)//"ne") ! 10^19/m^3
!	open(20, file=trim(profile),status="old",action="read")
!		do i = 1, nprf
!			read(20, *) ne(i)
!			!print *, ne(i)
!		end do
!	close!(!20)
!	
!	!pause
!
!	profile=trim(trim(profiles)//"Te") ! keV
!	open(20, file=trim(profile),status="old",action="read")
!		do i = 1, nprf
!			read(20, *) Te(i)
!			!print *, Te(i)
!		end do
!	close(20)
!
!	profile=trim(trim(profiles)//"z_eff") ! keV
!	open(20, file=trim(profile),status="old",action="read")
!		do i = 1, nprf
!			read(20, *) Zeff(i)
!			!print *, Zeff(i)
!		end do
!	close(20)

	! compute ni and nc by using ne and Zeff
!	do i = 1, nprf
!		ni(i) = ne(i)*(Zeff(i)-Zc)/(Zi*(Zi-Zc))
!		nc(i) = (ne(i)-ni(i)*Zi)/Zc
!		!print *, ni(i)*(3._dp/100._dp)
!		!print *, nc(i)
!	end do
!
!	profile=trim(trim(profiles)//"Ti_1") ! keV
!	open(20, file=trim(profile),status="old",action="read")
!		do i = 1, nprf
!			read(20, *) Ti(i)
!			!print *, (i)
!		end do
!	close(20)
!
!	profile=trim(trim(profiles)//"Ti_3") ! keV
!	open(20, file=trim(profile),status="old",action="read")
!		do i = 1, nprf
!			read(20, *) Tc(i)
!		end do
!	close(20)

	! save profiles to one file for use
	open(111,file='profiles.dat',status='replace', action='write')
		write(111,"(8A14)") 'rho', 'ne', 'Te', 'ni', 'Ti', 'nc', 'Tc', 'zeff'
		do i = 1, nprf
			write(111,204) vrho(i), ne(i), Te(i), ni(i), Ti(i), nc(i), Tc(i), (ni(i)*Zi+nc(i)*Zc**2)/(ni(i)*Zi+nc(i)*Zc)
		end do
	close(111)
	203 FORMAT (8A14)
	!204 FORMAT (8ES14.6)
	204 FORMAT (5e16.9)


	! transform betweent psip and normalized psit
	call inrcsnak(vpsi,q,cscoef)
	!vpsit(1) = 0.
	do i = 1, npsi
		vpsit(i) = inrcsitg(vpsi,cscoef,vpsi(1),vpsi(i))
	end do
	dum1 = vpsit(1)
	dum2 = vpsit(npsi)- vpsit(1)
	do i = 1, npsi
		vpsit(i)=(vpsit(i)-dum1)/dum2
		!vpsit(i) = (vpsit(i)-vpsit(1))/vpsit(npsi) ! normalization of psi_T at psi grid
		!print *, vpsit(i), dum3 
	end do
	!pause
	call inrcsnak(vpsit,sr0,cscoef)
	!call inrcsnak(vpsi,sr0,cspsisr)

	! check dpsi_t/dpsi = q
	!call inrcsnak(vpsi,vpsit,cscoef)
	!do i = 1, npsi
	!	print *, inrcsval(vpsi,cscoef,vpsi(i),1)*dum2,q(i)
	!end do

	! compute corresponding sr at vrho
	do i = 1, nprf
		vrho(i) = vrho(i)**2 ! this is normalized psi_T
	end do
	do i = 1, nprf
		!dum1 = inrcsval(vpsit,cscoef,vrho(i),0) ! corresponding r 
		!dum2 = inrcsval(vpsi,cspsisr,dum1,0) ! corresponding r 
		!vrrho(i) = dum1
		vrrho(i) = inrcsval(vpsit,cscoef,vrho(i),0)
		!print *, vrrho(i), Te(i), a0
	end do
	!pause
	vrrho(1) = 0.
	vrrho(nprf) = a0
	
	rhon2 = rhon**2
	r0m = inrcsval(vpsit,cscoef,rhon2,0)
	roa = r0m/a0
	r0l = (roa - dr0)*a0
	r0r = (roa + dr0)*a0
	!r0m = roa*a0	

	psil = inrcsval(sr0,cspsisr,r0l,0)
	psir = inrcsval(sr0,cspsisr,r0r,0)
	psim = inrcsval(sr0,cspsisr,r0m,0)
	dpdr = inrcsval(sr0,cspsisr,r0m,1)
	! get contours of these 3 flux surfaces
	allocate(RL(ntheta,3), ZL(ntheta,3), BP(ntheta))
	call ctp(psil,RL(:,1),ZL(:,1))
	call ctp(psim,RL(:,2),ZL(:,2))
	call ctp(psir,RL(:,3),ZL(:,3))

	! compute B_p at the reference flux surface
	do j = 1, ntheta
		dum1 = -cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,RL(j,2),ZL(j,2),0,1)/RL(j,2) ! B_Z
		dum2 = cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,RL(j,2),ZL(j,2),1,0)/RL(j,2) !B_R
		BP(j) = sqrt(dum1**2+dum2**2)
	end do

	!  compute miler equilibrium parameters
	R0 = (RL(1,2) + RL(((ntheta+1)/2),2))/2.
	call inrcsnak(vpsi,q,cscoef)
	q0 = inrcsval(vpsi,cscoef,psim,0)
	s0 = inrcsval(vpsi,cscoef,psim,1)*inrcsval(sr0,cspsisr,r0m,1)*r0m/q0
	call inrcsnak(sr0,br0,cscoef)
	shift = inrcsval(sr0,cscoef,r0m,1)
	!print *, s0
	!call inrcsnak(sr0,q,cscoef)  ! another way to check s0
	!s0 = inrcsval(sr0,cscoef,r0m,1)
	!print *, s0*r0m/q0
	!pause 
	! check q profile with initial q in g-file
	!do i = 1, npsi
		!print *, q(i), qpsi(i)
	!end do

	! compute the aera of the reference flux surface
	! s = int( R*dphi * dl), dl = sqrt(del_R^2+delZ^2)
	asf = 0.
	do j = 1, ntheta-1
	    dum1 = (RL(j,2)+RL(j+1,2))/2._dp ! R
	    dum2 = sqrt((RL(j,2)-RL(j+1,2))**2+(ZL(j,2)-ZL(j+1,2))**2) ! dl
	    asf = asf + dum1*dum2
	    !asf = asf + dum2
	end do 	
	asf = asf*twopi

	deallocate (cscoef)
	allocate(cscoef(4,nprf-1))

	call inrcsnak(vrrho,ne,cscoef)	
	dum1 = inrcsval(vrrho,cscoef,r0m,0)
	n0e = dum1  
	dum2 = inrcsval(vrrho,cscoef,r0m,1)
	RoLne = -(R0/(dum1/dum2))

	do i = 1, nprf
		vrho(i) = sqrt(vrho(i)) 
	end do
	call inrcsnak(vrho,ne,cscoef)
	dum1 = inrcsval(vrrho,cscoef,r0m,0)


	call inrcsnak(vrrho,Te,cscoef)
	dum1 = inrcsval(vrrho,cscoef,r0m,0)
	t0e = dum1
	dum2 = inrcsval(vrrho,cscoef,r0m,1)
	RoLte = -(R0/(dum1/dum2))

	!!=======================================
	!do i = 1, nprf
	!	dum1 = inrcsval(vrrho,cscoef,vrrho(i),0)
	!	dum2 = inrcsval(vrrho,cscoef,vrrho(i),1)
	!	!RoLte = -(R0/(dum1/dum2))
	!	dum3 = dum2/dum1*1.8
	!	!print *, vrrho(i), dum3
	!end do

	!! compute ni and nc, Zeff = 3, so that ni = 9 nc. Thus nc = 0.1 ne, ni = 0.9 ne
	!do i = 1, nprf
	!	ni(i) = 0.9*ne(i)
	!	nc(i) = 0.1*ne(i)
	!end do
	
	!! compute Ti, Tc. assume Ti(1) = 0.95keV, Ti has same scale profile with Te
	!Ti(1) = 0.95_sp
	!dum1 = 0.95_sp/Te(1)
	!!print *, dum1,'dum1'
	!do i = 1, nprf
	!	Ti(i) = Te(i)*dum1
	!	Tc(i) = Ti(i)
	!	!print *, Te(i), Ti(i), i
	!end do	

	call inrcsnak(vrrho,ni,cscoef)
	dum1 = inrcsval(vrrho,cscoef,r0m,0)
	n0i = dum1
	dum2 = inrcsval(vrrho,cscoef,r0m,1)
	RoLni = -(R0/(dum1/dum2))

	call inrcsnak(vrrho,Ti,cscoef)
	dum1 = inrcsval(vrrho,cscoef,r0m,0)
	t0i = dum1
	dum2 = inrcsval(vrrho,cscoef,r0m,1)
	RoLti = -(R0/(dum1/dum2))

	call inrcsnak(vrrho,nc,cscoef)
	dum1 = inrcsval(vrrho,cscoef,r0m,0)
	n0c = dum1
	dum2 = inrcsval(vrrho,cscoef,r0m,1)
	dum3 = dum1/dum3
	RoLnc = -(R0/(dum1/dum2))

	call inrcsnak(vrrho,Tc,cscoef)
	dum1 = inrcsval(vrrho,cscoef,r0m,0)
	t0c = dum1
	dum2 = inrcsval(vrrho,cscoef,r0m,1)
	RoLtc = -(R0/(dum1/dum2))	

	! compute profile parameters such Bu, teti, rhoia,etc
	! B_u = B(R_0,Z_0)
	teti = t0e/t0i

	deallocate(cscoef)
	!allocate(cscoef(4, (ntheta+1)/2-1))
	!call inrcsnak(RL(1:(ntheta+1)/2,2), ZL(1:(ntheta+1)/2,2), cscoef)
	!dum3 = inrcsval(RL(1:(ntheta+1)/2,2), cscoef, R0, 0)
	!B0 = cdbbsval2d(knotr,knotz,kr,kz,bscoefscalb,R0,dum3,0,0)
	!print *, B0,'B01'
	B0=inrcsval(gridpsi,cscoeffpol,psim,0)/R0 ! B0 is B_u for GEM
	!print *, B0, 'B02'
	n0e = n0e*1.E19 ! in /m^3
	T0e = T0e*1.E3 ! in eV

	n_u = n0e
	X_u = R0/1000.
	omega_u = E_charge*B0/mp
	t_u = 1._dp/omega_u
	v_u = x_u/t_u
	K_u = mp*V_u**2
	J_u = E_charge*n_u*V_u
	A_u = E_u/(E_charge*V_u)
	Phi_u = E_u/e_charge
	E_u = phi_u/x_u
	Ki_u = n_u*x_u*mp*v_u**2/t_u
	Q_u = n_u*x_u**3*mp*v_u**2/t_u

	cs = sqrt(T0e*E_charge/(2._sp*mp))
	rhoia = mimp*mp*cs/(E_charge*B0)/a0  ! 1d3*E_charge = kev, for deuterrim 
	betae = 2._sp*mu0*n0e*t0e*E_charge/B0**2
	!nue = (n0e*e_charge**4*17.)/(3.*twopi**1.5_dp*eps0**2*sqrt(me)*(t0e*1000._sp*E_charge)**1.5_dp)*a0/cs!
	nue = (n_u*e_charge**4*17.)/(4.*pi*eps0**2*sqrt(me)*(K_u)**1.5_dp)*a0/cs !normalized electron collision frequency
	!print *, nue, omega_u
	!pause
	! save flux data to compute delta' and kappa' 
	write(unit=gemfilename,fmt="(A20,'fs34dk',f3.2,'.dat')") trim(profiles),rhon
	gemfilename = trim(gemfilename)
	open(111,file=gemfilename,status='replace', action='write')
		write(111,101) ntheta
		write(111,104) dpdr, shift
		do j = 1, ntheta
			write(111,104) RL(j,2), ZL(j,2) -zmagx ! force z=0 at the magnetic axis
		end do
		do j = 1, ntheta
			write(111,103) BP(j) ! force z=0 at the magnetic axis
		end do
		!write(111,105) r0l, r0m, r0r
	close(111)

	! save parameters for gem.in
	
	write(unit=gemfilename,fmt="(A20,'g2g_rho',f3.2,'.dat')") trim(profiles),rhon
	gemfilename = trim(gemfilename)
	!print *, trim(gemfilename), rhon
	open(111,file=gemfilename,status='replace', action='write')
		write(111,"(4A11)"),'"mimp','mcmp','chgi','chgc"'
		write(111,"(4I11)"),int(mimp),int(mcmp),int(Zi),int(Zc)
		write(111,*)
		write(111,"(7A11)"), '"R/a','shift','q0','shat0','teti','tcti','rhoia"'
		write(111,"(7G11.4)")R0/a0,shift,q0,s0,teti,t0c/t0i,rhoia
		write(111,*)
		write(111,"(7A11)"),'"R/Lni','R/Lti',   'R/Lne',   'R/Lte', 'R/Lnc','R/Ltc', 'nc/ne"'
		write(111,"(7G11.4)")RoLni, RoLti, RoLne, RoLte, Rolnc, Roltc, n0c*1.E19/n0e
		write(111,*)
		write(111,"(A11,G11.4)"), 'r0a', r0m/a0
		write(111,"(A11,G11.4)"), 'betae', betae
		write(111,"(A11,G11.4)"), 'nue', nue
		write(111,*)
		write(111,"(A11,G11.4)"), 'ki_u*area', ki_u*asf
		write(111,*)
		write(111,303), 'B_unit', B0, 'Tesla B_t(R0,Z0)'
		write(111,303), 'n_unit', n_u, '/m^3'
		write(111,303), 'x_unit', x_u, 'm, R0/1000'
		write(111,303), 'omega_u', omega_u, 'Hz'
		write(111,303), 't_u', t_u, 's'
		write(111,303), 'V_u', v_u, 'm/s'
		write(111,303), 'K_u', K_u, 'J'
		write(111,303), 'J_u', J_u, 'A'
		write(111,303), 'm_u', mp, 'kg'	
		write(111,303), 'Ki_u', Ki_u, 'W/m^2'	
		write(111,303), 'Q_u', Q_u, 'W'	
		write(111,303), 'k_theta*rho_s/n', q0/(r0m/a0)*rhoia
		write(111,303), 'area of surface', asf, 'm^2'
		write(111,303), 'estimate area', twopi*R0*r0m*(twopi+4.*(0.37)), 'm^2'
		write(111,*), 'R0 and a0', R0, a0
		write(111,*)
		write(111,*), 'kappa, s_kappa, delta, s_delta will be produced by numerical fitting'	
	close(111)
		303 FORMAT (A16,G16.4,A16)
	stop

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
	n = npsi*2
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
					 !print *, 'warning for tt 1', i,j
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
					 !print *, 'warning for tt 1', i,j
					do kk= i, n
						tmppsi(kk) = tmppsi(i-1) + (kk-i+1)*(tmppsi(i-1)-tmppsi(i-2))
					end do
					exit
				end if
			end do
			!pause
			call inrcsnak(tmppsi,tmpl,tmpcscoef)
			
			ll = inrcsval(tmppsi,tmpcscoef,psin,0)
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

	!print *, tt1, tt2, tt3, tt4
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
	n = npsi*2
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
					 !print *, 'warning for tt 1', i,j
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
					 !print *, 'warning for tt 1', i,j
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
