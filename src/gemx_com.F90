module gemx_com
!common data used for gem
      use mpi
      USE pputil
      use iso_c_binding
implicit none

INTERFACE
  real(8) function revers(num,n)
  end function revers

  real(8) function ran2(i)
  end function ran2

  real(8) function en3(s)
      real(8) :: s
  end function en3


! subroutine new_gemx_com_c() bind(c, name='new_gemx_com_c_')
! end subroutine new_gemx_com_c
END INTERFACE


integer, bind(c) :: imx,jmx,kmx,mmx

integer, bind(c) :: nmx,nsmx,nsubd=8,ntube=4,petsc_color,petsc_rank,iBoltzmann,globle_integer=0,eBoltzmann,eAdiabatic,iterations,dbg
integer,dimension(0:10006):: rand_table
	 character*70 outname
	 REAL(8) :: endtm,begtm,pstm
	 REAL(8), bind(c) :: starttm,lasttm,tottm
         REAL(8), bind(c) :: start_total_tm, end_total_tm, start_integ_tm, end_integ_tm, start_ppush_tm, end_ppush_tm, start_cpush_tm, end_cpush_tm
         REAL(8), bind(c) :: total_tm = 0.0, integ_tm = 0.0, ppush_tm = 0.0, cpush_tm = 0.0
!          imx,jmx,kmx = max no. of grid pts in x,y,z
!          mmx         = max no. of particles
!          nmx         = max. no. of time steps
!          nsmx        = max. no. of species (including tracer particles

INTEGER,dimension(:),allocatable :: mm,tmm,lr
REAL(8),dimension(:),allocatable :: mims,q
INTEGER, bind(c) :: timestep,iez
integer, bind(c) ::iseed

real(8),dimension(:),allocatable :: time
REAL(8), bind(c) :: dx,dz,dzeta,pi,pi2,dt,totvol,n0,tcurr
REAL(8) :: etaohm
REAL(8), bind(c) :: lx,lz
INTEGER, bind(c) :: nm,nsm,ncurr,iflr,ifield_solver,ntracer,i3D,icollision
REAL(8), bind(c) :: cut,amp,tor,amie,emass,qel,rneu
INTEGER, bind(c) :: iput,iget,ision,isham,peritr,iadi
integer, bind(c) :: idg
real(8), dimension(:,:,:), allocatable :: phi_k, dphidr, dphi_kdr, d2phidr2, d2phi_kdr2, dphidz, dphi_kdz, d2phidz2, d2phi_kdz2, OPPphi, OPPphik, l_hand, r_hand  !!!!!!!!!! why these 3D? -zhichen


REAL(8), bind(c) :: vcut
integer, bind(c) :: nonlin,nonline,iflut,ifluid,ipara
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

INTEGER, bind(c) :: nplot,xnplt

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
integer, bind(c) :: numprocs
INTEGER, bind(c) :: Last,MyId, cnt , ierr
INTEGER, bind(c) :: GRID_COMM,TUBE_COMM, PETSC_COMM
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

!pointer declarations
type(c_ptr), bind(c) :: oppphi_ptr, oppphik_ptr
type(c_ptr), bind(c) :: tmm_ptr, mm_ptr, zeta2_ptr, x2_ptr, z2_ptr, mims_ptr, u2_ptr, mu_ptr, w2_ptr, x3_ptr, zeta3_ptr, z3_ptr, u3_ptr, w3_ptr
type(c_ptr), bind(c) :: ileft_ptr, xbackw_ptr, zbackw_ptr, jleft_ptr, iright_ptr, xforw_ptr, zforw_ptr, jright_ptr, lr_ptr, jac_ptr, rho_ptr, dene_ptr
type(c_ptr), bind(c) :: phi_ptr, ex_ptr, ez_ptr, ezeta_ptr, phi_k_ptr, dphidr_ptr, dphi_kdr_ptr, dphidz_ptr, dphi_kdz_ptr, d2phidr2_ptr, d2phi_kdr2_ptr, d2phidz2_ptr
type(c_ptr), bind(c) :: d2phi_kdz2_ptr, l_hand_ptr, r_hand_ptr, den2d2_ptr, q_ptr, den_ptr, upar_ptr, apars_ptr, apar_ptr, jpar_ptr, dden2d_ptr, den2d1_ptr
type(c_ptr), bind(c) :: rand_table_ptr, delbx_ptr, delby_ptr, delbz_ptr, phis_ptr, denes_ptr, upars_ptr, gn0e_ptr, gbtor_ptr, bmag_ptr, xg_ptr, gcpnex_ptr, gcpnez_ptr

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

   !1D Arrays
   tmm_ptr = c_loc(tmm(1))
   mm_ptr = c_loc(mm(1))
   zeta2_ptr = c_loc(zeta2(1))
   x2_ptr = c_loc(x2(1))
   z2_ptr = c_loc(z2(1))
   mims_ptr = c_loc(mims(1))
   u2_ptr = c_loc(u2(1))
   mu_ptr = c_loc(mu(1))
   w2_ptr = c_loc(w2(1))
   x3_ptr = c_loc(x3(1))
   zeta3_ptr = c_loc(zeta3(1))
   z3_ptr = c_loc(z3(1))
   u3_ptr = c_loc(u3(1))
   w3_ptr = c_loc(w3(1))
   q_ptr = c_loc(q(1))
   lr_ptr = c_loc(lr(1))
   jac_ptr = c_loc(jac(0))
   rand_table_ptr = c_loc(rand_table(0))
   xg_ptr = c_loc(xg(0))
   !2D Arrays
   ileft_ptr = c_loc(ileft(0,0))
   xbackw_ptr = c_loc(xbackw(0,0))
   zbackw_ptr = c_loc(zbackw(0,0))
   jleft_ptr = c_loc(jleft(0,0))
   iright_ptr = c_loc(iright(0,0))
   xforw_ptr = c_loc(xforw(0,0))
   zforw_ptr = c_loc(zforw(0,0))
   jright_ptr = c_loc(jright(0,0))
   den2d2_ptr = c_loc(den2d2(0,0))
   dden2d_ptr = c_loc(dden2d(0,0))
   den2d1_ptr = c_loc(den2d1(0,0))
   gn0e_ptr = c_loc(gn0e(0,0))
   bmag_ptr = c_loc(bmag(0,0))
   gbtor_ptr = c_loc(gbtor(0,0))
   gcpnex_ptr = c_loc(gcpnex(0,0))
   gcpnez_ptr = c_loc(gcpnez(0,0))
   !3D Arrays
   phi_ptr = c_loc(phi(0,0,0))
   ex_ptr = c_loc(ex(0,0,0))
   ez_ptr = c_loc(ez(0,0,0))
   ezeta_ptr = c_loc(ezeta(0,0,0))
   phi_k_ptr = c_loc(phi_k(0,0,0))
   dphidr_ptr = c_loc(dphidr(0,0,0))
   dphi_kdr_ptr = c_loc(dphi_kdr(0,0,0))
   dphidz_ptr = c_loc(dphidz(0,0,0))
   dphi_kdz_ptr = c_loc(dphi_kdz(0,0,0))
   d2phidr2_ptr = c_loc(d2phidr2(0,0,0))
   d2phi_kdr2_ptr = c_loc(d2phi_kdr2(0,0,0))
   d2phidz2_ptr = c_loc(d2phidz2(0,0,0))
   d2phi_kdz2_ptr = c_loc(d2phi_kdz2(0,0,0))
   oppphi_ptr = c_loc(OPPphi(0,0,0))
   oppphik_ptr = c_loc(OPPphik(0,0,0))
   l_hand_ptr = c_loc(l_hand(0,0,0))
   r_hand_ptr = c_loc(r_hand(0,0,0))
   upar_ptr = c_loc(upar(0,0,0))
   apars_ptr = c_loc(apars(0,0,0))
   apar_ptr = c_loc(apar(0,0,0))
   jpar_ptr = c_loc(jpar(0,0,0))
   rho_ptr = c_loc(rho(0,0,0))
   dene_ptr = c_loc(dene(0,0,0))
   delbx_ptr = c_loc(delbx(0,0,0))
   delby_ptr = c_loc(delby(0,0,0))
   delbz_ptr = c_loc(delbz(0,0,0))
   phis_ptr = c_loc(phis(0,0,0))
   denes_ptr = c_loc(denes(0,0,0))
   upars_ptr = c_loc(upars(0,0,0))
   !4D Arrays
   den_ptr = c_loc(den(1,0,0,0)) !first index uses 1 based index, the rest are like the rest of the 3D Arrays
   
      !call new_gemx_com_c()
end subroutine new_gemx_com

end module gemx_com
