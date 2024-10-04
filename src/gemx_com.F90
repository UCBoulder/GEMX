module gemx_com
!common data used for gem
      use mpi
      USE pputil
implicit none

INTERFACE
  real(8) function revers(num,n)
  end function revers

  real(8) function ran2(i)
  end function ran2

  real(8) function en3(s)
      real(8) :: s
  end function en3
END INTERFACE

integer :: imx,jmx,kmx,mmx,nmx,nsmx,nsubd=8,ntube=4,petsc_color,petsc_rank,iBoltzmann,globle_integer=0,eBoltzmann,eAdiabatic,iterations
integer,dimension(0:10006):: rand_table
	 character*70 outname
	 REAL(8) :: endtm,begtm,pstm
	 REAL(8) :: starttm,lasttm,tottm
         REAL(8) :: start_total_tm, end_total_tm, start_integ_tm, end_integ_tm, start_ppush_tm, end_ppush_tm, start_cpush_tm, end_cpush_tm
         REAL(8) :: total_tm = 0.0, integ_tm = 0.0, ppush_tm = 0.0, cpush_tm = 0.0
!          imx,jmx,kmx = max no. of grid pts in x,y,z
!          mmx         = max no. of particles
!          nmx         = max. no. of time steps
!          nsmx        = max. no. of species (including tracer particles

INTEGER,dimension(:),allocatable :: mm,tmm,lr
REAL(8),dimension(:),allocatable :: mims,q
INTEGER :: timestep,iseed,iez

real(8),dimension(:),allocatable :: time
REAL(8) :: dx,dz,dzeta,pi,pi2,dt,totvol,n0,tcurr
REAL(8) :: etaohm
REAL(8) :: lx,lz
INTEGER :: nm,nsm,ncurr,iflr,ifield_solver,ntracer,i3D,icollision
REAL(8) :: cut,amp,tor,amie,emass,qel,rneu
INTEGER :: iput,iget,idg,ision,isham,peritr,iadi
real(8), dimension(:,:,:), allocatable :: phi_k, dphidr, dphi_kdr, d2phidr2, d2phi_kdr2, dphidz, dphi_kdz, d2phidz2, d2phi_kdz2, OPPphi, OPPphik, l_hand, r_hand  !!!!!!!!!! why these 3D? -zhichen


REAL(8) :: vcut
integer :: nonlin,nonline,iflut,ifluid,ipara
COMPLEX(8) :: IU

REAL(8),DIMENSION(:,:,:,:),allocatable :: den
REAL(8),DIMENSION(:,:,:),allocatable :: rho
real(8),dimension(:,:,:),allocatable :: phi!,den_pre,dden

!Calder Edit
!real(8),dimension(:,:),allocatable :: phiavg
!Calder Edit End

REAL(8),DIMENSION(:,:,:),allocatable :: ex
REAL(8),DIMENSION(:,:,:),allocatable :: ez
REAL(8),DIMENSION(:,:,:),allocatable :: ezeta

REAL(8),DIMENSION(:,:,:),allocatable :: delbx,delbz,delby
REAL(8),DIMENSION(:),allocatable :: xg,zg,jac
real(8),dimension(:,:,:),allocatable :: apar,dene
real(8),dimension(:,:,:),allocatable :: upar
real(8),dimension(:,:,:),allocatable :: phis,denes,apars,upars

real(8),dimension(:,:,:),allocatable :: jpar

real(8),dimension(:,:),allocatable :: bmag,gbtor,gbx,gbz,gnuobx,gnuoby,gupae0,xforw,zforw,zbackw,xbackw,den2d1,den2d2,dden2d
real(8),dimension(:,:,:),allocatable :: dnedx,dnedy,dupadx,dupady
real(8),dimension(:,:),allocatable :: gn0i,gt0i,gn0e,gt0e,gcptex,gcptez,gcpnex,gcpnez

!     variables for tracing a grid (i,j) along field lien to the neighboring planes
integer,dimension(:,:),allocatable :: ileft,jleft,iright,jright
      
!          particle array declarations
REAL(8),DIMENSION(:),allocatable :: mu
REAL(8),DIMENSION(:),allocatable :: x2,zeta2,z2,u2
REAL(8),DIMENSION(:),allocatable :: x3,zeta3,z3,u3
REAL(8),DIMENSION(:),allocatable :: w2,w3


!              Various diagnostic arrays and scalars
!    plotting constants

INTEGER :: nplot,xnplt

!    energy diagnostic arrays

REAL(8),DIMENSION(:,:),allocatable :: ke
REAL(8),DIMENSION(:),allocatable :: fe,te
REAL(8),DIMENSION(:),allocatable :: rmsphi,rmsez,rmsapa,avewi
REAL(8),DIMENSION(:,:),allocatable :: nos

!    flux diagnostics
REAL(8),DIMENSION(:),allocatable :: vol
REAL(8),DIMENSION(:,:),allocatable :: efle,pfle
REAL(8),DIMENSION(:,:),allocatable :: pfl,efl

integer,parameter :: Master=0
integer :: numprocs
INTEGER :: MyId,Last,cnt,ierr
INTEGER :: GRID_COMM,TUBE_COMM, PETSC_COMM
INTEGER :: GCLR,TCLR,GLST,TLST
INTEGER :: stat(MPI_STATUS_SIZE)
INTEGER :: lngbr,rngbr,idprv,idnxt

character * (*) directory
parameter(directory='./dump/')

character * (*) outdir
parameter(outdir='./out/')

!real(8) :: ran2,revers
!integer :: mod
!real(8) :: amod
save

contains
subroutine new_gemx_com()

allocate(mm(nsmx),tmm(nsmx),lr(nsmx))
allocate(mims(nsmx),q(nsmx))
allocate(time(0:nmx))

      ALLOCATE( den(2,0:imx,0:jmx,0:kmx))
      
ALLOCATE( rho(0:imx,0:jmx,0:kmx))
allocate( phi(0:imx,0:jmx,0:kmx))!,den_pre(0:imx,0:jmx,0:kmx),dden(0:imx,0:jmx,0:kmx))
!allocate(phiavg(0:imx,0:jmx))!,0:kmx)) !Calder Edit


ALLOCATE( ex(0:imx,0:jmx,0:kmx)) 
ALLOCATE( ez(0:imx,0:jmx,0:kmx)) 
ALLOCATE( ezeta(0:imx,0:jmx,0:kmx))

ALLOCATE( delbx(0:imx,0:jmx,0:kmx),delbz(0:imx,0:jmx,0:kmx),delby(0:imx,0:jmx,0:kmx))
ALLOCATE( xg(0:imx),zg(0:jmx),den2d1(0:imx,0:jmx),den2d2(0:imx,0:jmx),dden2d(0:imx,0:jmx))
allocate( apar(0:imx,0:jmx,0:kmx),dene(0:imx,0:jmx,0:kmx))

allocate( upar(0:imx,0:jmx,0:kmx),jpar(0:imx,0:jmx,0:kmx))
allocate( upars(0:imx,0:jmx,0:kmx),phis(0:imx,0:jmx,0:kmx),&
          denes(0:imx,0:jmx,0:kmx),apars(0:imx,0:jmx,0:kmx))

allocate( jac(0:imx))
allocate( bmag(0:imx,0:jmx),gbtor(0:imx,0:jmx),gbx(0:imx,0:jmx),gbz(0:imx,0:jmx))
allocate( dnedx(0:imx,0:jmx,0:kmx),dnedy(0:imx,0:jmx,0:kmx), &
          dupadx(0:imx,0:jmx,0:kmx),dupady(0:imx,0:jmx,0:kmx))
allocate(gn0i(0:imx,0:jmx),gn0e(0:imx,0:jmx),gt0i(0:imx,0:jmx),gt0e(0:imx,0:jmx),xforw(0:imx,0:jmx),zforw(0:imx,0:jmx),zbackw(0:imx,0:jmx),xbackw(0:imx,0:jmx))
allocate(gcpnex(0:imx,0:jmx),gcpnez(0:imx,0:jmx),gcptex(0:imx,0:jmx),gcptez(0:imx,0:jmx))          
allocate(gnuobx(0:imx,0:jmx),gnuoby(0:imx,0:jmx),gupae0(0:imx,0:jmx)) 
allocate(ileft(0:imx,0:jmx),jleft(0:imx,0:jmx),iright(0:imx,0:jmx),jright(0:imx,0:jmx))

! Boltzmann Electron subroutine arrays for Newton solve
allocate(phi_k(0:imx,0:jmx,0:kmx),dphidr(0:imx,0:jmx,0:kmx),dphi_kdr(0:imx,0:jmx,0:kmx),d2phidr2(0:imx,0:jmx,0:kmx),d2phi_kdr2(0:imx,0:jmx,0:kmx),d2phidz2(0:imx,0:jmx,0:kmx),d2phi_kdz2(0:imx,0:jmx,0:kmx))
allocate(dphidz(0:imx,0:jmx,0:kmx),dphi_kdz(0:imx,0:jmx,0:kmx))
allocate(OPPphi(0:imx,0:jmx,0:kmx),OPPphik(0:imx,0:jmx,0:kmx),l_hand(0:imx,0:jmx,0:kmx),r_hand(0:imx,0:jmx,0:kmx))



!          particle array declarations
allocate( mu(1:mmx))
allocate( x2(1:mmx),zeta2(1:mmx),z2(1:mmx),u2(1:mmx))
allocate( x3(1:mmx),zeta3(1:mmx),z3(1:mmx),u3(1:mmx))
allocate( w2(1:mmx),w3(1:mmx))


ALLOCATE( ke(nsmx,0:nmx),fe(0:nmx),te(0:nmx))
ALLOCATE( rmsphi(0:nmx),rmsez(0:nmx),rmsapa(0:nmx),avewi(0:nmx))
ALLOCATE( nos(nsmx,0:nmx))

!    flux diagnostics
ALLOCATE(vol(1:nsubd),efle(1:nsubd,0:nmx),pfle(1:nsubd,0:nmx), &
         pfl(nsmx+1,0:nmx),efl(nsmx,0:nmx))


end subroutine new_gemx_com

end module gemx_com
