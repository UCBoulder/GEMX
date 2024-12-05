!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
MODULE pputil
!
      use mpi
      IMPLICIT NONE
  PRIVATE
  PUBLIC :: ppinit
!
  INTEGER, SAVE :: me, nvp,npp,GCLR,TCLR,p_color,p_rank

  INTEGER, SAVE :: TUBE_COMM,GRID_COMM, PETSC_COMM

CONTAINS
!======================================================================
!
 	SUBROUTINE ppinit(idproc,nproc,ntube,kmx,i3D,com1,com2,com_petsc,petsc_color,petsc_rank)
!
     use mpi
     INTEGER, INTENT(IN) :: ntube,kmx,i3D
     INTEGER, INTENT(OUT) :: nproc
	INTEGER, INTENT(OUT) :: idproc,com1,com2, com_petsc,petsc_color,petsc_rank
	INTEGER :: ierr,npp,n_tor_planes
!
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, npp, ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)
	nproc = npp
	IDPROC = me
! 
	GCLR=INT(me/ntube)
	TCLR=MOD(me,ntube)
!
	CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,GCLR,&
     	 &	TCLR,GRID_COMM,ierr)
	CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,TCLR,&
      &	GCLR,TUBE_COMM,ierr)

       if(i3D/=0) then
           n_tor_planes=kmx+1
        else
           n_tor_planes=1
        end if
        
           p_color= int(me*(n_tor_planes)/npp)
           p_rank = mod(me,npp/(kmx+1))
           CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,p_color,p_rank,PETSC_COMM,ierr)
           petsc_color=p_color
           petsc_rank=p_rank
!        else
!           CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,0,nproc,PETSC_COMM,ierr)
!
	com1=TUBE_COMM
        com2=GRID_COMM
        com_petsc=PETSC_COMM
	nvp=npp/ntube
!
	END SUBROUTINE ppinit
END MODULE pputil
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
