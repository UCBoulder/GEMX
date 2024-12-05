	 subroutine outd(n)
	 use gemx_com
         use equil      
	 implicit none
	 INTEGER :: n,np
	 integer :: i,j,k
         integer m


         if (myid==0) then
            write(*,*) "timestep =", timestep
!         end if
      
      
!SP      open(935, file='./out/tracer.out',status='unknown',position='append')
!SP     with no append file is overwritten
      if (n.eq.ncurr) open(935, file='./out/tracer.out',status='unknown')
      do m=1,ntracer
         write(935,*) timestep,m, (x3(m))*xu+Rgrid(0),(z3(m))*xu+Zgrid(0)!,u3(m)
      end do
      if (n.eq.nm)  close(935)


      open(unit=11, file = './out/testne',status='unknown',action='write')
      do j=0,jmx
         write(11,*) dene(:,j,0)
      enddo
      close(11)

 !     if (k==1) then
     
       open(unit=11, file = './out/testphi',status='unknown',action='write')
       do j=0,jmx
          write(11,*) phi(:,j,0)  !k is assumed to be 0
       enddo
       close(11)
      endif



	 if (mod(n,nplot).eq.0) then

!SEP	   call phixy(phi(:,:,:),'phixy',31,n)
!SEP	   call phixz(phi(:,:,:),'phixz',32,n)

!SEP	   call phixy(apar(:,:,:),'apaxy',36,n)
!SEP       call phixz(apar(:,:,:),'apaxz',37,n)

	 endif

	 return
	 end

!--------------------------------------------------------------

	 subroutine phixy(grd,fl,unt,n)

	 use gemx_com
	 implicit none
	 REAL(8) :: grd(0:imx,0:jmx,0:1)
	 integer :: i,j,k,n
	 INTEGER :: unt,tind,flag,oproc,outk
	 character*5 fl
	 character*70 flnm
!	 save flag

	 oproc=int(kmx/2/kmx*ntube)
	 if (MyId.eq.oproc) then

	 flnm=outdir//fl//'.out'
	 open(unt,file=flnm,form='formatted', &
     	   status='unknown',position='append')

	 if (n.eq.0) then
!      if (flag.ne.-1) then
	   tind=int(nm/nplot)+1
	   if(iget.eq.0)write(unt,110)imx,jmx,tind,nplot,int(dt)
	 endif
!       flag=-1
	 outk=0  !km/2-MyId*kcnt
	 write(unt,99)n
 99	 format('time step= ',I6)
	 do j=0,jmx
	   do i=0,imx
	     write(unt,100)grd(i,j,outk)
	   enddo
	 enddo


 100	 format (e10.3)
 110	 format (5I5)
	 endfile unt
	 close(unt)
	 endif

	 return
	 end

!--------------------------------------------------------------

	 subroutine phixz(grd,fl,unt,n)

	 use gemx_com
	 implicit none
	 integer :: i,j,k,n
	 REAL(8) :: grd(0:imx,0:jmx,0:1)
	 REAL(8) :: sbuf(0:(imx+1)*kmx)
	 REAL(8) :: rbuf(0:(imx+1)*kmx)
	 INTEGER :: unt,tind,flag
	 INTEGER :: lst,m,pkm
	 character*5 fl
	 character*70 flnm
!      save flag

	 flnm=outdir//fl//'.out'

!      if ( (MyId.eq.Last).and.(flag.ne.-1) ) then
	 if ( (MyId.eq.Last).and.(n.eq.0) ) then
	   open(unt,file=flnm,form='formatted', &
     	     status='unknown',position='append')
	   tind=int(nm/nplot)+1
	   if(iget.eq.0)write(unt,110)imx,kmx,tind,nplot,int(dt)
	   endfile unt
	   close(unt)
	 endif
!      flag=-1

	 pkm=int(ntube*kmx/numprocs)
	 do k=0,pkm-1
	   do i=0,imx
	     m=k*(imx+1)+i
	     sbuf(m)=grd(i,jmx/2,k)
	   enddo
	 enddo
	 cnt=pkm*(imx+1)
	 lst=cnt-1
	 call MPI_GATHER(sbuf(0:lst),cnt, &
               MPI_REAL8, &
               rbuf,cnt,MPI_REAL8, &
               glst,tube_comm,ierr)
! 

	 if (MyId.eq.Last) then
	    open(unt,file=flnm,form='formatted', &
                status='unknown',position='append')
	    
	    write(unt,99)n
 99	    format('time step= ',I6)

	    do i=0,((imx+1)*pkm)*numprocs/ntube-1
	       write(unt,100)rbuf(i)
	    enddo
	    do k=pkm,kmx
	       do i=0,imx
		  write(unt,100)grd(i,jmx/2,k)
	       enddo
	    enddo
	    endfile unt
	    close(unt)
	 endif

 100	 format (e10.3)
 110	 format (5I5)

	 return
	 end

!--------------------------------------------------------------
