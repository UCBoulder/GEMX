      program gemx
#include <petsc/finclude/petscksp.h>
      use gemx_com
      use equil

      use petsc
      use petscdmda
      use petscksp
      ! #include "fftw3.f"
!      use para_com



      implicit none

       integer :: status,mid_i,mid_j
       integer :: n,i,j,k,ip,m,outk,ix=135,jx=68
       integer :: iter, filter_int !Calder Edit

       real::random
       real :: tmp
       PetscInt is,js,iw,jw,idx,n_in_porcs
       PetscInt one,three,vec_start,vec_end
       PetscErrorCode petsc_ierr
       PetscScalar, POINTER ::phi_array(:)
       KSP ksp
       DM dm
       PetscObject  vec
       Vec petsc_phi,mpi_phi
       PetscViewer viewer
       VecScatter   ctx
       external ComputeRHS,ComputeMatrix,ComputeInitialGuess


         
!       call init
       call initialize
      


       
  outk=0!(kmx+1)/2
   
       one = 1
       three = 3


      if (eBoltzmann == 0) then !Calder Edit

       PETSC_COMM_WORLD =  PETSC_COMM
      

       
       PetscCallA(PetscInitialize(petsc_ierr))



       PetscCallA(KSPCreate(PETSC_COMM_WORLD,ksp,petsc_ierr))
       PetscCallA(DMDACreate2D(PETSC_COMM_WORLD, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,imx+1,jmx+1,PETSC_DECIDE,PETSC_DECIDE,one,one, PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, dm, petsc_ierr))
       PetscCallA(DMSetFromOptions(dm,petsc_ierr))
       PetscCallA(DMSetUp(dm,petsc_ierr))
       PetscCallA(KSPSetDM(ksp,dm,petsc_ierr))
       PetscCallA(KSPSetComputeInitialGuess(ksp,ComputeInitialGuess,0,petsc_ierr))
       PetscCallA(KSPSetComputeOperators(ksp,ComputeMatrix,0,petsc_ierr))      	
       PetscCallA(DMDAGetCorners(dm,is,js,PETSC_NULL_INTEGER,iw,jw,PETSC_NULL_INTEGER,petsc_ierr))
       PetscCallA(KSPSetFromOptions(ksp,petsc_ierr))
       PetscCallA(KSPSetUp(ksp,petsc_ierr))
      end if !Calder Edit

!  include "Initialize_petsc.h"
     

       if(iget.eq.0)call loadi
       call integ(2)

               if(myid==0)then
                open(unit=11, file = 'testden',status='unknown',action='write')
                do j=0,jmx                 
                  write(11,*) den2d2(:,j)
                enddo
                  close(11)
               end if
               if(i3D==0)then
                  do k=1,kmx
                     xn0i=den2d2
                  end do
               end if
               
       
        starttm=MPI_WTIME()
        upar=0

      !   mid_i=imx/2
      !   mid_j=jmx/2
      !   mid_i=257
      !   mid_j=257
      mid_i = (imx+1)/2
      mid_j = (jmx+1)/2
        tor_n=1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!initialize perturbation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (checkpoint == 0) then
        do k=0,kmx
              do i=0,imx
                 do j=0,jmx
!                    call random_number(random)
!                    dene(i,j,k)=mask(i,j)*cos(2*pi*k/(kmx+1))*2*exp(-((i-mid_i)**2+(j-mid_j)**2)/(0.09*min(mid_i,mid_j))**2)




                    !                    apar(i,j,k)=mask(i,j)*cos(tor_n*2*pi*k/(kmx+1))*2*exp(-((i-mid_i)**2+(j-mid_j)**2)/(0.09*min(mid_i,mid_j))**2)
                    ! apar(i,j,k)=mask2(i,j)*2*exp(-((i-ix)**2+(j-jx)**2)/(0.04*257)**2)!*cos(tor_n*2*pi*k/(kmx+1))
                    !                    apar(i,j,k)=mask2(i,j)*(exp(-(psi_p(i,j)-0.13)**2/0.01**2)-exp(-(psi_p(i,j)-0.16)**2/0.01**2))*cos(tor_n*2*pi*k/(kmx+1))
                   ! if (i==mid_i .and. j==mid_j) then
                   !    apar(i,j,k)=0
                   ! else
                      ! apar(i,j,k)=mask2(i,j)*(exp(-(psi_p(i,j)-0.1905)**2/0.01**2))*cos(tor_n*2*pi*k/(kmx+1))*(2*(j-mid_j)**2*dz**2/((j-mid_j)**2*dz**2+(i-mid_i)**2*dx**2)-1)!cos(2*pi*2*ATAN((j-mid_j)/(i-mid_i)))
                   ! endif
                    
                    
!                    apar(i,j,k)=mask(i,j)*(ran2(iseed)-0.5)
!                    apar(i,j,k)=0
                    apars(i,j,k)=0
!                    apar(i,j,k)=0
                     jpar(i,j,k)=0
                     dene(i,j,k)=0!ran(-0.5)
                     phi(i,j,k)=0!-j*0.01
                     ez(i,j,k)=0!0.01/dz
                     
                     ex(i,j,k)=0 !Calder Edits 1/9/25
                     ezeta(i,j,k)=0 !Calder Edits 1/9/25
                 enddo
              enddo
           enddo
         else
            open(unit=10, file = 'out/checkpoint_apar',status='old',action='read')
            read(10,*) apars
            close(10)
            
            open(unit=10, file = 'out/checkpoint_jpar',status='old',action='read')
            read(10,*) jpar
            close(10)

            open(unit=10, file = 'out/checkpoint_ne',status='old',action='read')
            read(10,*) dene
            close(10)

            open(unit=10, file = 'out/checkpoint_phi',status='old',action='read')
            read(10,*) phi
            close(10)

            open(unit=10, file = 'out/checkpoint_ex',status='old',action='read')
            read(10,*) ex
            close(10)

            open(unit=10, file = 'out/checkpoint_ez',status='old',action='read')
            read(10,*) ez
            close(10)

            open(unit=10, file = 'out/checkpoint_ezeta',status='old',action='read')
            read(10,*) ezeta
            close(10)
         end if


           phiavg=0 !Calder Edit

           call get_jpar(apar)
           call get_ne(0)

           if(i3d==0)then
              apar=0
              dene=0
              call integ(2)
           end if
           

           if (MyId==0) then

           open(unit=11, file = 'testj0',status='unknown',action='write')
               do j=0,jmx
                  
                  write(11,*) jpar(:,j,outk)
                  enddo
               close(11)

               open(unit=11, file = 'testne0',status='unknown',action='write')
               do j=0,jmx
                  
                  write(11,*) dene(:,j,outk)
                  enddo
               close(11)

               open(unit=11, file = 'testapar0',status='unknown',action='write')
               do j=0,jmx
                  
                  write(11,*) apar(:,j,outk)
                  enddo
               close(11)
               
               open(unit=11, file = 'testne0_zeta',status='unknown',action='write')
                 do k=0,kmx
                   write(11,*) dene(mid_i,mid_j,k)
                 enddo
               close(11)

              open(unit=11, file = 'testjpar0_zeta',status='unknown',action='write')
                 do k=0,kmx
                   write(11,*) jpar(mid_i,mid_j,k)
                 enddo
               close(11)
            end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end of init perturbation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
         ! do i=0,imx
            ! do j=0,jmx
            ! write(*,*)   -(xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_i(i,j)**2+c2_over_vA2(i,j)
            ! end do
         ! end do
 if (ifield_solver .eq. 1) then         
                 ncurr = 1                      
 end if

  start_total_tm = MPI_WTIME()
  do  timestep=ncurr,nm
     do i=0,10006
          if (ran2(iseed)-0.5>0) then
             rand_table(i)=1
          else
             rand_table(i)=-1
      end if
          
      end do
           tcurr = tcurr+dt

!	   call accumulate(timestep-1,0)
!	   call ezamp
!	   call gkps
!     call field(timestep-1,0)
   

   open(unit=11, file = 'testweights',status='unknown',position='append')
      weight_diag = 0
      do m=1,mm(1)         
         weight_diag =+ w2(m)/mm(1)
      enddo
      write(11,*) timestep, weight_diag
   close(11)


      if(ifield_solver .eq. 1 ) then

       phi = 0.0 
       phiavg = 0.0
       denes=dene

      if (eBoltzmann == 1) then
         k = 0
         do iter=0, iterations
            ! do iter=0, iterations
               ! call fluxavg(phi,phiavg)
            ! end do
            call boltzsolve(phi)
         end do
         ! call  MPI_Allreduce(MPI_IN_PLACE, phi, (imx+1)*(jmx+1)*(kmx+1),MPI_Real8, MPI_SUM, MPI_COMM_WORLD,ierr)


      else if (i3D /= 0) then
         phi=0
      do iter=0, iterations 
       do k=MyId*(kmx+1)/(numprocs),(MyId+1)*(kmx+1)/(numprocs)-1
         ! do iter=0, iterations
            call fluxavg(phi,phiavg) !Uses previous time step phi
               ! phi = 0 

         !  PetscCallA(KSPSetComputeRHS(ksp,ComputeRHS,k,petsc_ierr))
         PetscCallA(KSPSetComputeRHS(ksp,ComputeRHS,k,petsc_ierr))
         PetscCallA(KSPSolve(ksp,PETSC_NULL_VEC,PETSC_NULL_VEC,petsc_ierr))
         PetscCallA(KSPGetSolution(ksp,petsc_phi,petsc_ierr))
         PetscCall(VecGetArrayReadF90(petsc_phi, phi_array, petsc_ierr))
         PetscCallA(VecGetOwnershipRange(petsc_phi,vec_start,vec_end,petsc_ierr))

         
         do idx=1, vec_end-vec_start
            i=mod(idx-1,(iw))+is
            j=(idx-1)/(iw)+js
            phi(i,j,k)=phi_array(idx)!*mask(i,j)
         enddo
          PetscCall(VecRestoreArrayReadF90(petsc_phi,phi_array,petsc_ierr))
         enddo
         call  MPI_Allreduce(MPI_IN_PLACE, phi, (imx+1)*(jmx+1)*(kmx+1),MPI_Real8, MPI_SUM, MPI_COMM_WORLD,ierr)
         ! call fluxavg(phi,phiavg) !Turned off for CBC
       enddo
            
   else   
       k=0
       phi = 0
       do iter=0, iterations
         call fluxavg(phi,phiavg) !Turned off for CBC
         ! phi = 0

         PetscCallA(KSPSetComputeRHS(ksp,ComputeRHS,k,petsc_ierr))
         PetscCallA(KSPSolve(ksp,PETSC_NULL_VEC,PETSC_NULL_VEC,petsc_ierr))
         PetscCallA(KSPGetSolution(ksp,petsc_phi,petsc_ierr))
         PetscCallA(VecGetOwnershipRange(petsc_phi,vec_start,vec_end,petsc_ierr))
         PetscCall(VecGetArrayReadF90(petsc_phi, phi_array, petsc_ierr))

         do idx=1, vec_end-vec_start
            i=mod(idx-1,(iw))+is
            j=(idx-1)/(iw)+js
            phi(i,j,k)=phi_array(idx)!*mask(i,j)
         enddo

        PetscCall(VecRestoreArrayReadF90(petsc_phi,phi_array,petsc_ierr))
       
         call  MPI_Allreduce(MPI_IN_PLACE, phi, (imx+1)*(jmx+1)*(kmx+1),MPI_Real8, MPI_SUM, MPI_COMM_WORLD,ierr)

      end do !Calder Edit
      ! call  MPI_Allreduce(MPI_IN_PLACE, phi, (imx+1)*(jmx+1)*(kmx+1),MPI_Real8, MPI_SUM, MPI_COMM_WORLD,ierr)   
         do k=1,kmx !Put outside Calder Loops
            phi(:,:,k)=phi(:,:,0)
         end do
      end if
      
      ! do iter=0, iterations 
         ! call fluxavg(phi,phiavg)
      ! end do

      !PING Function
      if (timestep <= 10 .and. nonlin /= 1) then
         do i=0, imx
            do j=0, jmx
               do k=0, kmx
                  ! phi(i,j,k) = 100
                  phi(i,j,k) = 1e-8*cos(modes*((pi2*k)/(kmx+1)-1.3*atan2(Zgrid(j)-Zgrid(jmx/2),Rgrid(i)-Rgrid(imx/2))))* &
                  exp(-(sqrt((Zgrid(j)-Zgrid(jmx/2))**2+(Rgrid(i)-Rgrid(imx/2))**2)-0.25)**2/(2*0.05**2))!*cos(atan2(Zgrid(j)-Zgrid(jmx/2),Rgrid(i)-Rgrid(imx/2))/2)
                  ! phi(i,j,k) = 1e-5*cos(modes*(-1.3*atan2(Zgrid(j)+Zgrid(0)-Zgrid(jmx/2),Rgrid(i)+Rgrid(0)-Rgrid(imx/2)))) * &
                  ! exp(-((sqrt((Rgrid(i)+Rgrid(0)-Rgrid(imx/2))**2+(Zgrid(j)+Zgrid(0)-Zgrid(jmx/2))**2)-0.25)**2)/(2*(0.15)**2))
                  if (mask(i,j) < 0.99) then
                     phi(i,j,k) = 0
                  end if
               end do
            end do
         end do
      end if

      if (modes /= 0) then
         call fourier_modes(phi,modes)
      end if

      ! call binomial_filter(phi)
      
      if (filtering_iterations /= 0) then
         do filter_int=1, filtering_iterations !Edit 6/3/2025
            call poloidal_filter_methods(phi)
         end do
      end if

      !EFIELD TESTING
      if (cold_start == 1) then
         if (timestep >= 300) then
            call efieldcalc(phi)
         end if
      else
         call efieldcalc(phi)
      end if

      call get_apar(-1)
!         call  smooth(apars,2)
      call get_jpar(apars)
!           call smooth(jpar,3)
      call get_ne(-1)

!           if(myid==0)then
!               open(unit=11, file = 'testapars',status='unknown',action='write')
!               do j=0,jmx
                  
!                  write(11,*) apars(:,j,outk)-apars(:,j,outk+1)
!                  enddo
!                  close(11)

!               open(unit=11, file = 'testjpars',status='unknown',action='write')
!               do j=0,jmx
                  
!                  write(11,*) jpar(:,j,outk)-jpar(:,j,outk+1)
!                  enddo
!               close(11)   


 !              open(unit=11, file = 'testnes',status='unknown',action='write')
 !              do j=0,jmx
                  
 !                 write(11,*) denes(:,j,outk)-denes(:,j,outk+1)
 !                 enddo
 !                 close(11)

 !              open(unit=11, file = 'testBR',status='unknown',action='write')
 !              do j=0,jmx
                  
 !                 write(11,*) b0x(:,j)
 !                 enddo
 !                 close(11)
 !              end if


               if(ision==1)call ppush(timestep)
               if(ifluid==1)call integ(1)
               

          else
             if(ision==1)call ppush(timestep)
             !             if(ifluid==1)call pintef
             if(ifluid==1)call integ(1)

!             if(myid==0)then
!                open(unit=11, file = 'testden',status='unknown',action='write')
!                do j=0,jmx                 
!                  write(11,*) den2d2(:,j)
!                enddo
!                  close(11)
!              end if
                              
               
            endif
            
 ! write(*,*)'dx=', dx, 'dz=',dz


!	   call accumulate(timestep,1)
!	   call ezamp
!	   call gkps
!	   call field(timestep,1)

     if (ifield_solver .eq. 1) then
        phi=0.0

        if (eBoltzmann == 1) then
         k=0
         do iter=0, iterations
            ! do iter=0, iterations
            !    call fluxavg(phi,phiavg)
            ! end do
            call boltzsolve(phi)
         end do

         ! call  MPI_Allreduce(MPI_IN_PLACE, phi, (imx+1)*(jmx+1)*(kmx+1),MPI_Real8, MPI_SUM, MPI_COMM_WORLD,ierr)

       else if(i3D /= 0) then
         phi=0
      do iter=0,iterations
         ! call fluxavg(phi,phiavg)
         ! phi = 0 !Should already be zero 
         do k=MyId*(kmx+1)/(numprocs),(MyId+1)*(kmx+1)/(numprocs)-1  
       
         ! do iter=0,iterations
            call fluxavg(phi,phiavg)    

         PetscCallA(KSPSetComputeRHS(ksp,ComputeRHS,k,petsc_ierr))
         PetscCallA(KSPSolve(ksp,PETSC_NULL_VEC,PETSC_NULL_VEC,petsc_ierr))
         PetscCallA(KSPGetSolution(ksp,petsc_phi,petsc_ierr))
         PetscCallA(VecGetOwnershipRange(petsc_phi,vec_start,vec_end,petsc_ierr))
         PetscCall(VecGetArrayReadF90(petsc_phi, phi_array, petsc_ierr))

         do idx=1, vec_end-vec_start
            i=mod(idx-1,(iw))+is
            j=(idx-1)/(iw)+js
            phi(i,j,k)=phi_array(idx)!*mask(i,j)
         enddo
         PetscCall(VecRestoreArrayReadF90(petsc_phi,phi_array,petsc_ierr))
	      enddo
         call  MPI_Allreduce(MPI_IN_PLACE, phi, (imx+1)*(jmx+1)*(kmx+1),MPI_Real8, MPI_SUM, MPI_COMM_WORLD,ierr)
         ! call fluxavg(phi,phiavg) !Turned off for CBC
      enddo

      
      ! call  MPI_Allreduce(MPI_IN_PLACE, phi, (imx+1)*(jmx+1)*(kmx+1),MPI_Real8, MPI_SUM, MPI_COMM_WORLD,ierr)
           
   else

      k=0
      phi = 0
         do iter=0, iterations
            call fluxavg(phi,phiavg) !Turned off for CBC
               ! phi = 0
         
         PetscCallA(KSPSetComputeRHS(ksp,ComputeRHS,k,petsc_ierr))
         PetscCallA(KSPSolve(ksp,PETSC_NULL_VEC,PETSC_NULL_VEC,petsc_ierr))
         PetscCallA(KSPGetSolution(ksp,petsc_phi,petsc_ierr))
         PetscCall(VecGetArrayReadF90(petsc_phi, phi_array, petsc_ierr))

         do idx=1, vec_end-vec_start
            i=mod(idx-1,(iw))+is
            j=(idx-1)/(iw)+js
            phi(i,j,k)=phi_array(idx)!*mask(i,j)
         enddo
         
         
            

            PetscCall(VecRestoreArrayReadF90(petsc_phi,phi_array,petsc_ierr))

         call  MPI_Allreduce(MPI_IN_PLACE, phi, (imx+1)*(jmx+1)*(kmx+1),MPI_Real8, MPI_SUM, MPI_COMM_WORLD,ierr)

     end do !Calder Edit
   !   call  MPI_Allreduce(MPI_IN_PLACE, phi, (imx+1)*(jmx+1)*(kmx+1),MPI_Real8, MPI_SUM, MPI_COMM_WORLD,ierr)
         do k=1,kmx
            phi(:,:,k)=phi(:,:,0)
         end do
   end if

   ! do iter=0, iterations 
      ! call fluxavg(phi,phiavg)
   ! end do


   !PING Function
   if (timestep <= 10 .and. nonlin == 0) then
      do i=0, imx
         do j=0, jmx
            do k=0, kmx
               phi(i,j,k) = 1e-8*cos(modes*((pi2*k)/(kmx+1)-1.3*atan2(Zgrid(j)-Zgrid(jmx/2),Rgrid(i)-Rgrid(imx/2))))* &
                  exp(-(sqrt((Zgrid(j)-Zgrid(jmx/2))**2+(Rgrid(i)-Rgrid(imx/2))**2)-0.25)**2/(2*0.05**2))!*cos(atan2(Zgrid(j)-Zgrid(jmx/2),Rgrid(i)-Rgrid(imx/2))/2)
               ! phi(i,j,k) = 1e-5*cos(modes*(pi2*k-1.3*atan2(Zgrid(j)+Zgrid(0)-Zgrid(jmx/2),Rgrid(i)+Rgrid(0)-Rgrid(imx/2)))) * &
               ! exp(-((sqrt((Rgrid(i)+Rgrid(0)-Rgrid(imx/2))**2+(Zgrid(j)+Zgrid(0)-Zgrid(jmx/2))**2)-0.25)**2)/(2*(0.15)**2))
               if (mask(i,j) < 0.99) then
                  phi(i,j,k) = 0
               end if
            end do
         end do
      end do
   end if
   
   if (modes /= 0) then
      call fourier_modes(phi,modes)
   end if

   !!!!!!!!!!!!!!!!!!!!!!!
   ! call ftcamp(phi,timestep)
   !!!!!!!!!!!!!!!!!!!!!!!
   ! call binomial_filter(phi)

   if (filtering_iterations /= 0) then
      do filter_int=1, filtering_iterations !Edit
         call poloidal_filter_methods(phi)
      end do
   end if


   !EFIELD TESTING
   if (cold_start == 1) then
      if (timestep >= 300) then
         call efieldcalc(phi)
      end if
   else
      call efieldcalc(phi)
   end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PHI DIAGNOSTIC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   phi_diag = 0.0
   phi_diag_freq = 0.0

   do i = 0, imx
      do j = 0, jmx
         do k = 0, kmx
            phi_diag = phi_diag + (abs(phi(i,j,k))**2)/(imx*jmx*kmx)
            phi_diag_freq = phi_diag_freq + phi(i,j,k)/(imx*jmx*kmx)
         end do
      end do
   end do

   open(unit=11, file = 'testPhiDiag',status='unknown',position='append')
   write(11,*) timestep, phi_diag
   close(11)

   open(unit=11, file='testPhiFreq',status='unknown',position='append')
   write(11,*) timestep, phi_diag_freq
   close(11)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (i3D == 1) then
      ! call growthdiag(phi)
   end if

               if(MyId==0 .and. mod(timestep,10)==0)then
                  write(*,*)'outk=',outk
                  ! write(*,*) phi(:,:,outk)

               open(unit=11, file = 'testphi',status='unknown',action='write')
               do j=0,jmx
                  
                  write(11,*) phi(:,j,outk)
               enddo
               
               close(11)

               end if
            

    call get_apar(1)
!    call smooth(apar,2)
      call get_jpar(apar)
!      call smooth(jpar,3)
      call get_ne(1)


      if(MyId==0 .and. mod(timestep,1)==0)then

      open(unit=11, file = 'testphiavg',status='unknown',action='write')
      do j=0,jmx
         write(11,*) phiavg(:,j)
      enddo
      close(11)

      open(unit=11, file = 'testER', status='unknown',action='write')
      do j=0, jmx
         write(11,*) ex(:,j,0)
      end do
      close(11)

      open(unit=11,file='testEZ', status='unknown',action='write')
      do j=0,jmx
         write(11,*) ez(:,j,0)
      end do
      close(11)
      end if

      

       if(ision==1)call cpush(timestep)
        !        if(ifluid==1)call cintef(timestep)
       if(ifluid==1)call integ(2)
    else
        if(ision==1)call cpush(timestep)
        !        if(ifluid==1)call cintef(timestep)
        if(ifluid==1)call integ(2)
        !        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     end if
     

        if(myid==0 .and. mod(timestep,10)==0)then
           open(unit=11, file = 'testden2',status='unknown',action='write')
           do j=0,jmx
              write(11,*) den2d2(:,j)
           enddo
           close(11)
           open(unit=11, file = 'testdiffden',status='unknown',action='write')
           do j=0,jmx
              write(11,*) dden2d(:,j)
           enddo
           open(unit=11, file = 'testupar',status='unknown',action='write')
           do j=0,jmx
              write(11,*) upar(:,j,0)
           enddo
           close(11)
        end if
        

    

     call outd(timestep)

         if(MyId==0 .and. ifield_solver==1 .and. mod(timestep,10)==0) then    !5/30/25
               open(unit=11, file = 'testapar',status='unknown',action='write')
               do j=0,jmx
                  
                  write(11,*) apar(:,j,outk)
                  enddo
                  close(11)

               open(unit=11, file = 'testjpar',status='unknown',action='write')
               do j=0,jmx
                  
                  write(11,*) jpar(:,j,outk)
                  enddo
               close(11)   


               open(unit=11, file = 'testne',status='unknown',action='write')
               do j=0,jmx
                  
                  write(11,*) dene(:,j,outk)
                  enddo
                  close(11)

               open(unit=11, file = 'testphi_r_phi',status='unknown',action='write')
               do k=0,kmx
                  
                  write(11,*) phi(:,mid_j,k)
                  enddo
               close(11)
            end if
            




           if(myid.eq.master .and. mod(timestep,xnplt)==0)then
              open(9,file='plot',status='unknown',position='append')
              m = 2
              i = 4
              write(9,10)timestep,(x2(m+i*j),z2(m+i*j),j=1,7)
 10           format(1x,i6,16(1x,e10.3))
              close(9)
           end if

           if(myid.eq.master .and. ifield_solver.eq.1)then
            !   open(unit=11, file = 'testAtPhitrhotjt',status='unknown',position='append')                
            ! !   write(11,*)apar(387,256,outk),phi(387,256,outk),dene(387,256,outk),jpar(387,256,outk)
            !   close(11)

              if(mod(timestep,100)==0)then
                 open(unit=11, file = 'testPhit',status='unknown',position='append')
                 write(11,*)phi(:,:,outk)
                 
               !   open(unit=11, file = 'testPhi_major_diagnostic',status='unknown',position='append')
               !    write(11,*)phi(:,:,k)
              end if
              
                 
          
           write(*,*)'time_step=', timestep
           write(*,*)'dx=', dx, 'dz=',dz,'dzeta=',dzeta,'omega_A0=', tor_n/(Rgrid(mid_i)/xu*sqrt(c2_over_vA2(mid_i,mid_j)))
           write(*,*)'v_A=', 1/sqrt(c2_over_vA2(mid_i,mid_j)), 'Omega_i=', q(1)*b0(mid_i,mid_j)/mims(1)
        end if
        

 end do
end_total_tm = MPI_WTIME()
total_tm = total_tm + end_total_tm - start_total_tm
        call MPI_reduce(ppush_tm, tmp, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if(myid==0)ppush_tm = tmp/real(numprocs)
        call MPI_reduce(cpush_tm, tmp, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if(myid==0)cpush_tm = tmp/real(numprocs)
        call MPI_reduce(integ_tm, tmp, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if(myid==0)integ_tm = tmp/real(numprocs)
        call MPI_reduce(total_tm, tmp, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if(myid==0)total_tm = tmp/real(numprocs)
        if(myid==0)open(123, file = "gemx_timing.txt", status = "replace", action="write")
        if(myid==0)write(123,*)'ppush time', ppush_tm, 'cpush time', cpush_tm, 'integ time', integ_tm, 'total time', total_tm, 'other (including field solver)', total_tm - ppush_tm - cpush_tm- integ_tm
        if(myid==0)call flush(123)
        if(myid==0)close(123)
               
	 lasttm=MPI_WTIME()
  tottm=lasttm-starttm


      if (eBoltzmann == 0) then !Calder Edit
         PetscCallA(PetscFinalize(petsc_ierr))
      end if !Calder Edit

      !Set Checkpoints
      write(*,*) "Last timestep reached, writing checkpoints."

      open(unit=11, file = 'out/checkpoint_apar',status='unknown',action='write')
      do i=0,imx
         do j=0,jmx
            do k=0,kmx
               write(11,*) apar(i,j,k)
            end do
         end do
      end do
      close(11)

      open(unit=11, file = 'out/checkpoint_jpar',status='unknown',action='write')
      do i=0,imx
         do j=0,jmx
            do k=0,kmx
               write(11,*) jpar(i,j,k)
            end do
         end do
      end do
      close(11)   

      open(unit=11, file = 'out/checkpoint_ne',status='unknown',action='write')
      do i=0,imx
         do j=0,jmx
            do k=0,kmx
               write(11,*) dene(i,j,k)
            end do
         end do
      end do
      close(11)

      open(unit=11, file = 'out/checkpoint_phi',status='unknown',action='write')
      do i=0,imx
         do j=0,jmx
            do k=0,kmx
               write(11,*) phi(i,j,k)
            end do
         end do
      end do
      close(11)

      open(unit=11, file = 'out/checkpoint_ex',status='unknown',action='write')
      do i=0,imx
         do j=0,jmx
            do k=0,kmx
               write(11,*) ex(i,j,k)
            end do
         end do
      end do
      close(11)

      open(unit=11, file = 'out/checkpoint_ez',status='unknown',action='write')
      do i=0,imx
         do j=0,jmx
            do k=0,kmx
               write(11,*) ez(i,j,k)
            end do
         end do
      end do
      close(11)

      open(unit=11, file = 'out/checkpoint_ezeta',status='unknown',action='write')
      do i=0,imx
         do j=0,jmx
            do k=0,kmx
               write(11,*) ezeta(i,j,k)
            end do
         end do
      end do
      close(11)

      open(unit=11, file = 'out/checkpoint_x3',status='unknown',action='write')
      do m=1,mm(1)
         write(11,*) x3(m)
      end do
      close(11)

      open(unit=11, file = 'out/checkpoint_z3',status='unknown',action='write')
      do m=1,mm(1)
         write(11,*) z3(m)
      end do
      close(11)

      open(unit=11, file = 'out/checkpoint_zeta3',status='unknown',action='write')
      do m=1,mm(1)
         write(11,*) zeta3(m)
      end do
      close(11)

      open(unit=11, file = 'out/checkpoint_u3',status='unknown',action='write')
      do m=1,mm(1)
         write(11,*) u3(m)
      end do
      close(11)

      open(unit=11, file = 'out/checkpoint_w3',status='unknown',action='write')
      do m=1,mm(1)
         write(11,*) w3(m)
      end do
      close(11)

 100     call MPI_FINALIZE(ierr)
         end program gemx
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine init
      
      use gemx_com
      use equil
      implicit none
      character*(62) dumchar
      INTEGER :: i,j,k,n,ns,idum,i1,k1,m,j1
      INTEGER :: mm1,lr1
      REAL(8) :: x,z,dum,zdum
      REAL(8) :: dbdrp,dbdtp,bfldp,upae0p,dnuobdrp,dnuobdtp,btorp,bxp,bzp
      REAL(8) :: gn0ip,gn0ep,gt0ip,gt0ep,capnxp,capnzp
      REAL(8) :: wx0,wx1,wz0,wz1,b

      IU=cmplx(0.,1.)
      pi=4.0*atan(1.0)
      pi2 = pi*2.
      open(115,file='gemx.in')
      read(115,*) dumchar
      read(115,*) imx,jmx,kmx,mmx,nmx,nsmx,ntube
      read(115,*) dumchar
      read(115,*) dt,nm,nsm,iez
      read(115,*) dumchar
      read(115,*) iput,iget,ision,peritr
      read(115,*) dumchar
      read(115,*) nplot,xnplt
      read(115,*) dumchar
      read(115,*) cut,amp,tor
      read(115,*) dumchar
      read(115,*) etaohm
      read(115,*) dumchar
      read(115,*) ifluid,amie,rneu
      read(115,*) dumchar
      read(115,*) beta,nonlin,vcut
      read(115,*) dumchar
      read(115,*) ntracer,ifield_solver,i3D,iBoltzmann,eAdiabatic,iterations,icollision
      read(115,*) dumchar    !Calder Edit: eBoltzmann
      read(115,*) eBoltzmann !Calder Edit: eBoltzmann
      read(115,*) dumchar
      read(115,*) iflr,PADE
      read(115,*) dumchar
      read(115,*) CST !Load CST
      read(115,*) dumchar
      read(115,*) weightscheme !Option for full-f or delta-f
      read(115,*) dumchar
      read(115,*) modes,filtering_iterations !input as integer
      read(115,*) dumchar
      read(115,*) cold_start
      read(115,*) dumchar
      read(115,*) psi_max,psi_min,R_min,Z_min, Z_internal, psi_div,psi_a
      read(115,*) dumchar
      read(115,*) checkpoint
      close(115)
      
      nsm=1

      call new_gemx_com()
      ns = 1
      tmm(ns)=mmx!ntracer
!      mm(ns)=int(ntracer/numprocs)
      mm(ns)=mmx
      mims(ns)=2.0*1.67e-27
      q(ns)=1.0*1.6e-19
      lr(ns)=4

      emass = 1./amie
      qel = -1

      call new_equil()
      lx = xdim
      lz = zdim

      if(myid.eq.master)then
         open(9,file='plot',status='unknown',position='append')
         write(9,*)'a,rmaj0,lx,lz= ',a,rmaj0,lx,lz
         write(9,*)'xctr,xdim=',xctr,xdim         
         close(9)
      end if

      iadi = 0

      if(iget.eq.1) amp=0.

      dx=lx/real(imx)
      dz=lz/real(jmx)
      dzeta=pi2/(kmx+1)
!     
      do 10 i=0,imx
         xg(i)=i*dx 
 10   continue
      do 14 k=0,jmx
         zg(k)=k*dz
 14   continue

      do i1 = 0,imx
         x = i1*dx+xctr-xdim/2
         i = int(x/dxeq)
         i = min(i,nx-1)
         wx0 = ((i+1)*dxeq-x)/dxeq
         wx1 = 1.-wx0

         do k1 = 0,jmx
            z = k1*dz
            k = int(z/dzeq)
            k = min(k,nz-1)            
            wz0 = ((k+1)*dzeq-z)/dzeq
            wz1 = 1-wz0

            bfldp = wx0*wz0*b0(i,k)+wx0*wz1*b0(i,k+1) &
                 +wx1*wz0*b0(i+1,k)+wx1*wz1*b0(i+1,k+1) 
            btorp = wx0*wz0*b0zeta(i,k)+wx0*wz1*b0zeta(i,k+1) &
                 +wx1*wz0*b0zeta(i+1,k)+wx1*wz1*b0zeta(i+1,k+1) 
            bxp = wx0*wz0*b0x(i,k)+wx0*wz1*b0x(i,k+1) &
                 +wx1*wz0*b0x(i+1,k)+wx1*wz1*b0x(i+1,k+1) 
            bzp = wx0*wz0*b0z(i,k)+wx0*wz1*b0z(i,k+1) &
                 +wx1*wz0*b0z(i+1,k)+wx1*wz1*b0z(i+1,k+1) 
            gt0ip = wx0*wz0*t0i(i,k)+wx0*wz1*t0i(i,k+1) &
                 +wx1*wz0*t0i(i+1,k)+wx1*wz1*t0i(i+1,k+1) 
            gt0ep = wx0*wz0*t0e(i,k)+wx0*wz1*t0e(i,k+1) &
                 +wx1*wz0*t0e(i+1,k)+wx1*wz1*t0e(i+1,k+1) 
            gn0ip = wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) &
                 +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1) 
            gn0ep = wx0*wz0*xn0e(i,k)+wx0*wz1*xn0e(i,k+1) &
                 +wx1*wz0*xn0e(i+1,k)+wx1*wz1*xn0e(i+1,k+1) 
            capnxp = wx0*wz0*capnex(i,k)+wx0*wz1*capnex(i,k+1) &
                 +wx1*wz0*capnex(i+1,k)+wx1*wz1*capnex(i+1,k+1) 
            capnzp = wx0*wz0*capnez(i,k)+wx0*wz1*capnez(i,k+1) &
                 +wx1*wz0*capnez(i+1,k)+wx1*wz1*capnez(i+1,k+1) 

            b=1.-tor+tor*bfldp
            bmag(i1,k1) = b
            gbtor(i1,k1) = btorp
            gbx(i1,k1) = bxp
            gbz(i1,k1) = bzp                        

            gt0i(i1,k1) = gt0ip
            gt0e(i1,k1) = gt0ep
            gn0e(i1,k1) = gn0ep
            gn0i(i1,k1) = gn0ip
            gcpnex(i1,k1) = capnxp
            gcpnez(i1,k1) =  capnzp           

            gupae0(i1,k1) = upae0p
!            gnuoby(i1,k1) = (-dydrp*dnuobdtp+r0/q0*qhatp*dnuobdrp)*fp/radiusp*grcgtp
!            gnuobx(i1,k1) = dnuobdtp*fp/radiusp*grcgtp
         end do
      end do


      iseed = -(1777+myid*13)
      idum = ran2(iseed)
      phi = 0.
      apar = 0.
      dene = 0.
      upar = 0.


      do i = 0,imx
         do j = 0,jmx
            do k = 0,kmx
               phi(i,j,k) = amp*(ran2(idum)-0.5)*ifluid*1e-8  !amp*sin(nzcrt*xg(i)*pi/lx) !
               dene(i,j,k) = amp*(ran2(idum)-0.5)*ifluid *1e-8
               apar(i,j,k) = amp*(ran2(idum)-0.5)*ifluid *1e-10 
            end do
         end do
      end do

      if(myid.eq.master)then
         open(9,file='plot',status='unknown',position='append')
         write(9,*)'inner,outer = ',xctr-xdim/2,xctr+xdim/2
         write(9,*)'mi,qi=',mims(1),q(1)
         write(9,*)'mm(1)=',mm(1)         
         close(9)
      end if

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grad(ip)
  
!  currently set up for periodic in x,y,z

      use gemx_com
      use equil
      implicit none
      INTEGER :: i,j,k,ip
      real(8) :: ux(0:imx,0:jmx,0:1),uy(0:imx,0:jmx,0:1)
      real(8) :: tmp(0:imx,0:jmx,0:1),uoverb(0:imx,0:jmx,0:1)
      real(8) :: v(0:imx-1),dum,dum1

      call gradu(phi(:,:,:),ux,uy)
      ex(:,:,:) = -ux(:,:,:)
      ez(:,:,:) = -uy(:,:,:)

      delbx = 0.
      delby = 0.
      if(ifluid.eq.1)then
         call gradu(apar(:,:,:),ux,uy)
         delbx(:,:,:) = uy(:,:,:)
         delby(:,:,:) = -ux(:,:,:)
      end if

      call gradu(tmp(:,:,:),ux,uy)
      dnedx(:,:,:) = ux(:,:,:)
      dnedy(:,:,:) = uy(:,:,:)
      do i = 0,imx
         do j = 0,jmx
            do k = 0,1
               uoverb(i,j,k) = upar(i,j,k)  !/bfld(i,k)
            end do
         end do
      end do
      call gradu(uoverb(:,:,:),ux,uy)
      dupadx(:,:,:) = ux(:,:,:)
      dupady(:,:,:) = uy(:,:,:)

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid1(ip,n)

!    source quantities are are calculated: n_i
!    right now only ion quantitities are calculated...

      use gemx_com
      use equil
      implicit none
      REAL(8) :: x,z
      INTEGER :: m,n,i,j,k,l,ns,ip
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,ter
      REAL(8) :: wght,r,b,bfldp,dv
      REAL(8) :: xt,zt,rhog,pidum,vpar,xs,vfac
      real(8) :: myden(0:imx,0:jmx,0:kmx),myjpar(0:imx,0:jmx,0:kmx)
      REAL(8) :: rhox(4),rhoy(4)

      ns=1
      rho=0.
      den=0.
      jpar = 0.
      myden = 0.
      myjpar = 0.

      do m=1,mm(1)
         dv=float(lr(1))*(dx*dz*dzeta)

         x=x3(m)
         i = int(x/dxeq)
         wx0 = ((i+1)*dxeq-x)/dxeq
         wx1 = 1.-wx0

         z = z3(m)
         k = int(z/dzeq)
         wz0 = ((k+1)*dzeq-z)/dzeq
         wz1 = 1-wz0


         bfldp = wx0*wz0*b0(i,k)+wx0*wz1*b0(i,k+1) &
                 +wx1*wz0*b0(i+1,k)+wx1*wz1*b0(i+1,k+1) 
         ter = wx0*wz0*t0i(i,k)+wx0*wz1*t0i(i,k+1) &
                 +wx1*wz0*t0i(i+1,k)+wx1*wz1*t0i(i+1,k+1) 

         b=1.-tor+tor*bfldp

         rhog=sqrt(2.*b*mu(m)*mims(1))/(q(1)*b) * iflr

         rhox(1) = rhog
         rhoy(1) = 0.
         rhox(2) = -rhox(1)
         rhoy(2) = -rhoy(1)
         rhox(3) = 0
         rhoy(3) = rhog
         rhox(4) = 0
         rhoy(4) = -rhoy(3)

         vfac=0.5*(mims(1)*u3(m)**2 + 2.*mu(m)*b )
         wght=w3(m)/dv
         if(vfac/ter > vcut)wght=0.
         vpar = u3(m)

!    now do 1,2,4 point average, where lr is the no. of points...
         do 100 l=1,lr(1)
            xs=x3(m)+rhox(l) !rwx(1,l)*rhog
            zt=z3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=modulo(xs,xdim)
            zt=modulo(zt,zdim)

            include "gridli.h"
 100     continue
      enddo
if(idg.eq.1)write(*,*)myid,'pass ion grid1'
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!   enforce periodicity

      do 110 i=0,imx
         do 120 j=0,jmx
            do 130 k=0,kmx
               den(1,i,j,k)=q(ns)*myden(i,j,k)/n0/jac(i)
               jpar(i,j,k) = q(ns)*myjpar(i,j,k)/n0/jac(i)*ifluid
 130        continue
 120     continue
 110  continue


      do 150 i=0,imx
         do 160 j=0,jmx
            do 170 k=0,kmx
               rho(i,j,k)=rho(i,j,k)+den(1,i,j,k)
 170        continue
 160     continue
 150  continue

 499  continue
      do i = 0,imx
         do j = 0,jmx
            do k = 0,kmx
               rho(i,j,k) = ision*rho(i,j,k) + dene(i,j,k)*qel/ntube
            enddo
         enddo
      enddo      

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!        Normal distribution random no. generator, stand. dev. = 1.
!        Version 2 does it Willy's way...

         subroutine parperp(vpar,vperp2,m,pi,cnt,MyId)

         REAL(8) :: vpar,vperp2,r1,r2,t,pi
         INTEGER :: m,iflag,cnt,MyId
         REAL(8) :: c0,c1,c2
         REAL(8) :: d1,d2,d3
         data c0,c1,c2/2.515517,0.802853,0.010328/
         data d1,d2,d3/1.432788,0.189269,0.001308/


          r1=revers(m+MyId*cnt,7)
          r2=revers(m+MyId*cnt,11)


!.....quiet start---see denavit pf '71(?) & abramowitz hand book
!.....fibonacci start---see denavit comm. pla. phy. & con. fus. '81
! warning: we have g1=1 in the x-direction. This surpresses all odd
!          modes in the x-direction!!!

         iflag=1
         if(r1.le.0.5) go to 110
         r1=1.-r1
         iflag=-1
  110    continue
         if(r1.ge.1.e-6) then
           t=sqrt(log(1.0/(r1*r1)))
         else
           t=5.0
           write(*,*)'parperp2 warning  m= ',m
         endif

         vpar=t-(c0+c1*t+c2*t**2)/(1.+d1*t+d2*t**2+d3*t**3)
         vpar=vpar*iflag

          vperp2=-2.0*dlog(r2)

        return
        end

!---------------------------------------------------------------------- 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gkps()
      use gemx_com
      use equil
      use petsc
      use petscdmda
      use petscksp
      

      return
      end subroutine gkps
      !real,dimension(nx,nz,nzeta)::nepredict,aparpredic
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_apar(flagnumber)
         use gemx_com
         use equil

         IMPLICIT NONE

         integer::flagnumber, i, j, k
!         real:: E_par

         
         if (flagnumber == -1) then
            
            apars=apar-0.5*dt*(gradpar(phi))!+0.04*(jpar+q(1)*mu0*upar)) !+Epar)
            else
               apar = apar -dt*(gradpar(phi))!+0.04*(jpar+q(1)*mu0*upar)) !+Epar)    
            endif
          

          CONTAINS

            function  gradparz(matrix)
              real,dimension(0:imx,0:jmx,0:kmx)::matrix,gradparz
              call gradz(matrix,gradparz)
              do k=0,kmx
                 do i=0,imx
                    do j=0,jmx
                       gradparz(i,j,k)=0.5*gradparz(i,j,k)*mask2(i,j)
                    end do
                 end do
              end do
              
            end function gradparz
            



      function  gradpar(matrix)
        real,dimension(0:imx,0:jmx,0:kmx)::matrix,gradpar

        gradpar=0
      do k=1,(kmx-1)
         do i=2,(imx-2)
            do j=2,(jmx-2)
               if (mask2(i,j)<1.99)then
               else
               gradpar(i,j,k)=(b0x(i,j)/b0(i,j)*(matrix(i+1,j,k)-matrix(i-1,j,k))*0.5/dx  &
               +b0z(i,j)/b0(i,j)*(matrix(i,j+1,k)-matrix(i,j-1,k))*0.5/dz                &
               +b0zeta(i,j)/b0(i,j)*(matrix(i,j,k+1)-matrix(i,j,k-1))/((Rgrid(i)/xu)*2*dzeta))
               endif
            enddo
         enddo
      enddo

         do i=2,(imx-2)
            do j=2,(jmx-2)
               if (mask2(i,j)<1.99)then
               else
               gradpar(i,j,0)=(b0x(i,j)/b0(i,j)*(matrix(i+1,j,0)-matrix(i-1,j,0))*0.5/dx  &
               +b0z(i,j)/b0(i,j)*(matrix(i,j+1,0)-matrix(i,j-1,0))*0.5/dz                &
               +b0zeta(i,j)/b0(i,j)*(matrix(i,j,1)-matrix(i,j,kmx))/((Rgrid(i)/xu)*2*dzeta))
               endif
            enddo
         enddo


         do i=2,(imx-2)
            do j=2,(jmx-2)
               if (mask2(i,j)<1.99)then
               else
               gradpar(i,j,kmx)=(b0x(i,j)/b0(i,j)*(matrix(i+1,j,kmx)-matrix(i-1,j,kmx))*0.5/dx  &
               +b0z(i,j)/b0(i,j)*(matrix(i,j+1,kmx)-matrix(i,j-1,kmx))*0.5/dz                        &
               +b0zeta(i,j)/b0(i,j)*(matrix(i,j,0)-matrix(i,j,kmx-1))/((Rgrid(i)/xu)*2*dzeta))
               endif
            enddo
         enddo
      
      return
      end function gradpar
      end subroutine get_apar


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine get_jpar(matrix)
         use gemx_com
         use equil
         IMPLICIT NONE
         real, dimension(0:imx,0:jmx,0:kmx):: matrix
         integer::i,j,k,k0

         k0=0
         
        do k=0,kmx
            do i=2,imx-2
               do j=2,jmx-2
                  if (mask3(i,j)<2.99)then
                  else                     
                      jpar(i,j,k)=(-(matrix(i+1,j,k)+matrix(i-1,j,k)-2*matrix(i,j,k))/dx**2       &
                                 -(matrix(i,j+1,k)+matrix(i,j-1,k)-2*matrix(i,j,k))/dz**2   &
                                 -(matrix(i+1,j,k)-matrix(i-1,j,k))*0.5/(dx*Rgrid(i)/xu))  &
                                 -q(1)*mu0*upar(i,j,k)
                  endif
               enddo
            enddo
         enddo

         
         end subroutine get_jpar

 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         subroutine get_ne(flagnumber)
            use gemx_com
            use equil

            IMPLICIT NONE

            integer::flagnumber, i, j, k
            
            
            if (flagnumber == -1) then
               denes=dene+0.5*dt*(gradpar(jpar)) 
               elseif (flagnumber ==1) then
                  dene = dene + dt*(gradpar(jpar))
               elseif (flagnumber ==0) then
                  do i=2,imx-2
                     do j=2,jmx-2
                        do k=0,kmx
                           if(mask4(i,j)==4)then
                              dene(i,j,k)=jpar(i,j,k)*sqrt(c2_over_vA2(i,j))
                           end if
                           
                        enddo
                     enddo               
                  enddo      
               endif
               if (i3D==0) then
                  do k=1,kmx
                     if (flagnumber ==-1)then
                        denes(:,:,0)=denes(:,:,0)+denes(:,:,k)
                     else
                        dene(:,:,0)=dene(:,:,0)+dene(:,:,k)
                     end if
                  end do
                  if (flagnumber ==-1)then
                     denes=denes/(kmx+1)
                  else
                     dene=dene/(kmx+1)
                  end if
               end if
               
               
  
   
         CONTAINS
           function  gradparz(matrix)
           real,dimension(0:imx,0:jmx,0:kmx)::matrix,gradparz
           call gradz(matrix,gradparz)

              do k=0,kmx
                 do i=0,imx
                    do j=0,jmx
                       gradparz(i,j,k)=0.25*gradparz(i,j,k)*mask4(i,j)
                    end do
                 end do
              end do
         end function gradparz
         
         

         function  gradpar(matrix)
           real,dimension(0:imx,0:jmx,0:kmx)::matrix,gradpar

           gradpar=0

         do k=1,(kmx-1)
            do i=2,(imx-2)
               do j=2,(jmx-2)
                  if (mask4(i,j)<3.99)then
                  else                
                  gradpar(i,j,k)=(b0x(i,j)/b0(i,j)*(matrix(i+1,j,k)-matrix(i-1,j,k))*0.5/dx  &
                  +b0z(i,j)/b0(i,j)*(matrix(i,j+1,k)-matrix(i,j-1,k))*0.5/dz                &
                  +b0zeta(i,j)/b0(i,j)*(matrix(i,j,k+1)-matrix(i,j,k-1))/(Rgrid(i)/xu*2*dzeta))
                  endif
               enddo
            enddo
         enddo
   
            do i=2,(imx-2)
               do j=2,(jmx-2)
                  if (mask4(i,j)<3.99)then
                  else                     
                  gradpar(i,j,0)=(b0x(i,j)/b0(i,j)*(matrix(i+1,j,0)-matrix(i-1,j,0))*0.5/dx  &
                  +b0z(i,j)/b0(i,j)*(matrix(i,j+1,0)-matrix(i,j-1,0))*0.5/dz                &
                  +b0zeta(i,j)/b0(i,j)*(matrix(i,j,1)-matrix(i,j,kmx))/(Rgrid(i)/xu*2*dzeta))!*mask(i,j)
                  endif
               enddo
            enddo
   


         do i=2,(imx-2)
            do j=2,(jmx-2)
               if (mask4(i,j)<3.99)then
               else                  
               gradpar(i,j,kmx)=(b0x(i,j)/b0(i,j)*(matrix(i+1,j,kmx)-matrix(i-1,j,kmx))*0.5/dx  &
               +b0z(i,j)/b0(i,j)*(matrix(i,j+1,kmx)-matrix(i,j-1,kmx))*0.5/dz                        &
               +b0zeta(i,j)/b0(i,j)*(matrix(i,j,0)-matrix(i,j,kmx-1))/(Rgrid(i)/xu*2*dzeta))!*mask(i,j)
               endif
            enddo
         enddo


           
    !     return
         end function gradpar
       end subroutine get_ne

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      subroutine eqmo(ip)
      use gemx_com
      use equil
      implicit none
      integer :: i,j,k,ip
      real(8) :: eta

      ez(:,:,:) = 0.
      if(iez==0)return

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine spec(n)
      use gemx_com
      use equil
      implicit none
      integer :: i,j,k,l,m,n

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ezamp

      use gemx_com
      use equil

      implicit none

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(8) function ran2(idum)
      parameter( IM1=2147483563,  &
                IM2=2147483399, &
                AM=1.0/IM1,&
                IMM1=IM1-1,&
                IA1=40014,&
                IA2=40692,&
                IQ1=53668,&
                IQ2=52774,&
                IR1=12211,&
                IR2=3791,&
                NTAB=32,&
                NDIV=1+IMM1/NTAB,&
                EPS=1.2e-7,&
                RNMX=1.0-EPS &
               )
      integer :: j,k,idum2=123456789,iy=0,iv(0:NTAB-1)
      real(8) :: temp

      save idum2, iy,iv
      if(idum.le.0)then
         if(-idum.lt.1)then
            idum=1
         else
            idum = -idum
         end if
         idum2 = idum
         do j = NTAB+7,0,-1
            k = idum/IQ1
            idum = IA1*(idum-k*IQ1)-k*IR1
            if(idum.lt.0)idum = idum+IM1
            if(j.lt.NTAB)iv(j) = idum
         end do
         iy = iv(0)
      end if

      k = idum/IQ1
      idum = IA1*(idum-k*IQ1)-k*IR1
      if(idum.lt.0)idum = idum+IM1
      k = idum2/IQ2
      idum2 = IA2*(idum2-k*IQ2)-k*IR2
      if(idum2.lt.0)idum2 = idum2+IM2
      j = iy/NDIV
      iy = iv(j)-idum2
      iv(j) = idum
      if(iy<1)iy = iy+IMM1
      temp = AM*iy
      if(temp>RNMX)then
         ran2 = RNMX
      else
         ran2 = temp
      end if
      return
    end function ran2
    


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine loadi

      use gemx_com
      use equil
      implicit none
      INTEGER :: i,j,k,m,idum,ns,m1
      REAL(8) :: vpar,vperp2,r,x,z,b,ter,bfldp
      REAL(8) :: avgv,myavgv,avgw,myavgw
      real(8) :: dumx,dumy,dumz,jacp,rand(4)
      REAL(8) :: wx0,wx1,wz0,wz1,avex=0

      cnt=int(tmm(1)/numprocs)
      cnt=mmx

      myavgv=0.
      avgv=0.
      avgw = 0.
      myavgw = 0.

      if (CST /= 0) then
         open(unit=10,file='rdata.dat',status='old',action='read')
         read(10,*) x2
         close(10)
   
         open(unit=10,file='zdata.dat',status='old',action='read')
         read(10,*) z2
         close(10)
      end if 

      m=1
      do while(m<=mm(1))

!     load a slab of ions...

!         dumx=xdim*(ran2(iseed)+0.01)*0.9
!         dumy=zdim*(ran2(iseed)+0.01)*0.9

         !revers(MyId*cnt+j,2) !ran2(iseed)
         dumx=2*dxeq+(xdim-4*dxeq)*ran2(iseed)  !revers(MyId*cnt+j,2) !ran2(iseed)
         dumy=2*dzeq+(zdim-4*dzeq)*ran2(iseed) !revers(MyId*cnt+j,3) !ran2(iseed)
         dumz=pi2*ran2(iseed) !revers(MyId*cnt+j,5) !ran2(iseed)


!         dumx=dxeq+(xdim-2*dxeq)*m/((mm(1)))

         if (CST /= 0) then
            dumx = x2(m)-xctr+xdim/2
            dumy = z2(m)-zctr+zdim/2
         end if
         
         r = xctr-xdim/2+dumx
         jacp = r/(xctr+xdim/2)
!         if(ran2(iseed)<jacp)then
!            x2(m)=min(dumx,xdim-dxeq)
!            z2(m)=min(dumy,zdim-dzeq)
!            x2(m)=max(dumx,dxeq)
!            z2(m)=max(dumz,dzeq)
            zeta2(m)=dumz
            x2(m)=dumx
            z2(m)=dumy
            call parperp(vpar,vperp2,m,pi,cnt,MyId)

            x=x2(m)
            i = int(x/dxeq)
            wx0 = ((i+1)*dxeq-x)/dxeq
            wx1 = 1.-wx0

            z = z2(m)
            k = int(z/dzeq)
            wz0 = ((k+1)*dzeq-z)/dzeq
            wz1 = 1-wz0

            bfldp = wx0*wz0*b0(i,k)+wx0*wz1*b0(i,k+1) &
                   +wx1*wz0*b0(i+1,k)+wx1*wz1*b0(i+1,k+1) 
            ter = wx0*wz0*t0i(i,k)+wx0*wz1*t0i(i,k+1) &
                   +wx1*wz0*t0i(i+1,k)+wx1*wz1*t0i(i+1,k+1) 

            u2(m)=vpar/sqrt(mims(1)/ter)
            mu(m)=0.5*vperp2/bfldp*ter

            myavgv=myavgv+u2(m)

!    LINEAR: perturb w(m) to get linear growth...
!            w2(m)=2.*amp*ran2(iseed)
            
               w2(m)= (wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) &
                 +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1))*r/xctr*((imx-3)*(jmx-3)*(kmx+1))/(numprocs*mmx)!*xctr/(x+xctr-xdim/2.)
               gw(m) = 1
!            w2(m) = r/xctr*((imx-1)*(jmx-1)*(kmx+1))/(numprocs*mmx)
            if (weightscheme == 1) then
               if (nonlin == 1) then 
               w2(m) = 0
                  ! w2(m) = 1e-12
                  ! w2(m) = 1e-12*cos(modes*(zeta2(m)-1.3*atan2(dumy+Zgrid(0)-Zgrid(jmx/2),dumx+Rgrid(0)-Rgrid(imx/2)))) * &
                  ! exp(-((sqrt((dumx+Rgrid(0)-Rgrid(imx/2))**2+(dumy+Zgrid(0)-Zgrid(jmx/2))**2)-0.25)**2)/(2*(0.05)**2))
               else
                  ! w2(m) = 1e-12 
               !Initilizing weight as w_0 = A_i \cos{n(\zeta - q \arctan{z/r})} \exp{-\frac{(\sqrt{r^2+z^2}-a/2)^2}{2\Delta r^2}}
               w2(m) = 1e-12*cos(modes*(zeta2(m)-1.3*atan2(dumy+Zgrid(0)-Zgrid(jmx/2),dumx+Rgrid(0)-Rgrid(imx/2)))) * &
               exp(-((sqrt((dumx+Rgrid(0)-Rgrid(imx/2))**2+(dumy+Zgrid(0)-Zgrid(jmx/2))**2)-0.25)**2)/(2*(0.05)**2))
               end if
               gw(m) =  (wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) + wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1))*r/xctr*((imx-3)*(jmx-3)*(kmx+1))/(numprocs*mmx)!*xctr/(x+xctr-xdim/2.)
            end if
 
            myavgw=myavgw+w2(m)
            m = m+1            
         end do

         if (CST /= 0) then
            open(unit=12, file = 'testinitialv_par',status='unknown',action='write')
            write(12,*) u2
            close(12)
   
            open(unit=14, file = 'testenergeez',status='unknown',action='write')
            write(14,*) mu*bfldp + 0.5*mims(1)*u2**2
            close(14)
         end if 
!             do i=1,mmx
!                avex=avex+x2(i)
!             end do
!         write(*,*)avex/mmx
             if (MyId==0) then

         ! open(unit=11, file = 'testdepo_posi',status='unknown',action='write')
         !       do j=mmx-10000,mmx
                  
         !          write(11,*) x2(j),z2(j),zeta2(j)
         !       enddo
         !       close(11)
             end if

      
      myavgw = myavgw/mm(1)

      call MPI_ALLREDUCE(myavgv,avgv,1, &
          MPI_REAL8, &
          MPI_SUM,MPI_COMM_WORLD,ierr)
      if(idg.eq.1)write(*,*)'all reduce'
      avgv=avgv/float(tmm(1))
      do 180 m=1,mm(1)
         u2(m)=u2(m)-avgv
         x3(m)=x2(m)
         z3(m)=z2(m)
         zeta3(m)=zeta2(m)
         u3(m)=u2(m)
!         w2(m) = w2(m)-myavgw
         w3(m)=w2(m)
 180  continue

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gradu(u,ux,uz)
      use gemx_com
      use equil
      implicit none
      real(8) :: u(0:imx,0:jmx,0:1)
      real(8) :: ux(0:imx,0:jmx,0:kmx),uz(0:imx,0:jmx,0:kmx)
      integer :: i,j,k,l,m,n,jj,ju,jl
      real(8) :: ydum,wy1,ul

      do j=0,jmx-1
         ju = j+1
         jl = j-1
         if(j.eq.0)jl = jmx-1
         do i=0,imx-1
            do k=0,kmx
               uz(i,j,k)=(u(i,ju,k)-u(i,jl,k))/(2.*dz)
            enddo
         enddo
      enddo

      do i=1,imx-1
         do j=0,jmx-1
            do k=0,kmx
               ux(i,j,k)=(u(i+1,j,k)-u(i-1,j,k))/(2.*dx)
            enddo
         enddo
      enddo

! do boundary i=0
      do j=0,jmx-1
         do k=0,kmx
            ul=u(imx-1,j,k)
            ux(0,j,k)=(u(1,j,k)-ul)/(2.*dx)
         enddo
      enddo

      return
    end subroutine gradu
    
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gradz(u,uz)
      use gemx_com
      use equil
      implicit none
      real(8) :: u(0:imx,0:jmx,0:kmx)
      real(8) :: uz(0:imx,0:jmx,0:kmx)
      integer :: i,j,k,kleft,kright
      real(8) :: wx0,wx1,wz0,wz1,uleft,uright

      uz = 0.
      do k = 0,kmx
         kleft = k-1
         if(k==0)kleft=kmx
         kright = k+1
         if(k==kmx)kright = 0
         do i = 1,imx-1
            do j=1,jmx-1
               wx0 = ((ileft(i,j)+1)*dx-xbackw(i,j))/dx
               wx1 = 1.0-wx0
               wz0 = ((jleft(i,j)+1)*dz-zbackw(i,j))/dz
               wz1 = 1.0-wz0
               uleft = wx0*wz0*u(ileft(i,j),jleft(i,j),kleft) &
                      +wx1*wz0*u(ileft(i,j)+1,jleft(i,j),kleft) &
                      +wx0*wz1*u(ileft(i,j),jleft(i,j)+1,kleft) &
                      +wx1*wz1*u(ileft(i,j)+1,jleft(i,j)+1,kleft)
               wx0 = ((iright(i,j)+1)*dx-xforw(i,j))/dx
               wx1 = 1.0-wx0
               wz0 = ((jright(i,j)+1)*dz-zforw(i,j))/dz
               wz1 = 1.0-wz0
               uright = wx0*wz0*u(iright(i,j),jright(i,j),kright) &
                      +wx1*wz0*u(iright(i,j)+1,jright(i,j),kright) &
                      +wx0*wz1*u(iright(i,j),jright(i,j)+1,kright) &
                      +wx1*wz1*u(iright(i,j)+1,jright(i,j)+1,kright)

               uz(i,j,k)=(uright-uleft)/(2.*b0(i,j)/b0zeta(i,j)*dzeta*Rgrid(i)/xu)
            enddo
         enddo
      enddo
 !     uz(:,:,kmx) = uz(:,:,0)

      return
    end subroutine gradz
    






      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine initialize
         use gemx_com
      use equil

	implicit none
        real(8) :: dum,dum1,dum2,jacp,xndum,r,wx0,wx1
!        complex(8),dimension(0:1) :: x,y
        real(8),dimension(0:1) :: x,y
	integer :: n,i,j,k,ip

        call ppinit(MyId,numprocs,ntube,kmx,i3D,TUBE_COMM,GRID_COMM, PETSC_COMM,petsc_color,petsc_rank)
         
!     reset timestep counter.
         Last=numprocs-1
         timestep=0
         tcurr = 0.

         do i=0,Last
            if (MyId.eq.i) call init
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         enddo
      
         dum = 0.
         do i = 0,imx-1
            dum = dum+(jac(i)+jac(i+1))/2
         end do
         call MPI_ALLREDUCE(dum,jacp,1,  &
             MPI_REAL8,MPI_SUM,           &
             tube_comm,ierr)
         totvol = lx*lz*pi2*xctr    
         n0=float(tmm(1))/totvol
         

!         do k=0,kmx
!            den_pre(:,:,k)=xn0i
!         end do
         
         
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         ncurr = 1

	end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine accumulate(n,ip)
         use gemx_com
         use equil
	implicit none

	integer :: n,i,j,k,ip
	call grid1(ip,n)
	if(idg.eq.1)write(*,*)'pass grid1'
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine accumulate
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine field(n,ip)
         use gemx_com
         use equil
	implicit none
        integer :: n,i,j,k,ip,i1
        real(8) :: lbfr(0:imx,0:jmx)
        real(8) :: lbfs(0:imx,0:jmx)
        real(8) :: rbfr(0:imx,0:jmx)
        real(8) :: rbfs(0:imx,0:jmx)
        real(8) :: dum
        real(8) :: myrmsphi,rmp(20),myavap(0:imx-1)

	call grad(ip)

        call eqmo(ip)

end subroutine field

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pintef
      use gemx_com
      use equil
      implicit none

      integer :: i,j,k,ip
      real*8 :: dum,dumphi,dum1,dum2,denel,dener,aparl,aparr,ddedt
      real*8 :: ux(0:imx,0:jmx,0:kmx),uz(0:imx,0:jmx,0:kmx)

      phis = phi
      denes = dene
      apars = apar

      call gradz(upar,uz)
      do k = 0,kmx-1
         do i = 1,imx-1
            do j = 1,jmx-1
               ddedt = -uz(i,j,k)*gn0e(i,j)*gbtor(i,j)/((xctr-xdim/2+xg(i))*bmag(i,j))  &
                +gn0e(i,j)*(gcpnex(i,j)*ez(i,j,k)-gcpnez(i,j)*ez(i,j,k))/bmag(i,j)
               dene(i,j,k) = denes(i,j,k)+0.5*dt*ddedt
            end do
         end do
      end do

      call gradz(phi,uz)
      do k = 0,kmx-1
         do i = 1,imx-1
            do j = 1,jmx-1
               ddedt = -uz(i,j,k)*gbtor(i,j)/((xctr-xdim/2+xg(i))*bmag(i,j))
               apar(i,j,k) = apars(i,j,k)+0.5*dt*(ddedt+ezeta(i,j,k))
            end do
         end do
      end do

      return
    end subroutine pintef
    

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cintef(n)
      use gemx_com
      use equil
      implicit none

      integer :: i,j,k,n,ip
      real*8 :: tmpa(0:imx,0:jmx,0:1),tmpd(0:imx,0:jmx,0:1),tmp(0:imx,0:jmx,0:1)
      real*8 :: dum,dumphi,dum1,dum2,denel,dener,aparl,aparr
      REAL*8 :: myrmsapa
      real*8 :: dmnl1(0:imx,0:jmx,0:1),dmnl2(0:imx,0:jmx,0:1),dmnl3(0:imx,0:jmx,0:1),dmnl4(0:imx,0:jmx,0:1)
      real*8 :: ux(0:imx,0:jmx,0:1),uy(0:imx,0:jmx,0:1),exnl(0:imx,0:jmx,0:1)
      real(8) :: mydbr(0:imx-1),v(0:imx-1)


      
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine weight
  
      use gemx_com
      use equil
      implicit none
      INTEGER :: i,j,k,m=10,i1,j1
      real(8) :: ux(0:imx,0:jmx,0:1),uy(0:imx,0:jmx,0:1)
      real(8) :: x,z,zeta,dzt,rdum,dum1,wz0,wz1,wx0,wx1
      real(8) :: bp,bxp,bzp,radiusp,btorp


     
      dzt = dzeta/float(m)
      do i = 1,imx-1
!         x = xg(i)
!         if (i==257)write(*,*)xg(i)
         do j = 1,jmx-1
            x= xg(i)
            z = zg(j)
            do k = 0,m-1
               i1 = int(x/dxeq)
               j1 = int(z/dzeq)

               i1 = int(x/dxeq)
               wx0 = ((i1+1)*dxeq-x)/dxeq
               wx1 = 1.-wx0

               j1 = int(z/dzeq)
               wz0 = ((j1+1)*dzeq-z)/dzeq
               wz1 = 1-wz0

               bp = wx0*wz0*b0(i1,j1)+wx0*wz1*b0(i1,j1+1) &
                 +wx1*wz0*b0(i1+1,j1)+wx1*wz1*b0(i1+1,j1+1) 
               bxp = wx0*wz0*b0x(i1,j1)+wx0*wz1*b0x(i1,j1+1) &
                 +wx1*wz0*b0x(i1+1,j1)+wx1*wz1*b0x(i1+1,j1+1) 
               bzp = wx0*wz0*b0z(i1,j1)+wx0*wz1*b0z(i1,j1+1) &
                 +wx1*wz0*b0z(i1+1,j1)+wx1*wz1*b0z(i1+1,j1+1) 
               btorp = wx0*wz0*b0zeta(i1,j1)+wx0*wz1*b0zeta(i1,j1+1) &
                 +wx1*wz0*b0zeta(i1+1,j1)+wx1*wz1*b0zeta(i1+1,j1+1) 
               radiusp = xctr-xdim/2+x
 !              if(i==257.and.j==257) then
 !              write(*,*)x
 !              write(*,*)bxp
 !              write(*,*)bzp
 !              write(*,*)btorp
 !              write(*,*)radiusp
 !              write(*,*)dx
 !              write(*,*)dz
 !              write(*,*)dzt*radiusp*bxp/btorp
 !              end if
            
               x = x+dzt*radiusp*bxp/btorp
 !             if(i==257.and.j==257) then
 !              write(*,*)x-xg(i)
 !            end if
               x = min(x,lx)
               x = max(x,0.)
 !            if(i==257.and.j==257) then
 !              write(*,*)x-xg(i)
 !            end if
            
               z = z+dzt*radiusp*bzp/btorp
               z = min(z,lz)
               z = max(z,0.)
            end do
            iright(i,j) = int(x/dx)
            jright(i,j) = int(z/dz)
            xforw (i,j) =x
            zforw(i,j)=z
  !          write(*,*)i,iright(i,j),j,jright(i,j)
         end do
      end do

      dzt = -dzeta/float(m)
      do i = 1,imx-1
!         x = xg(i)
         do j = 1,jmx-1
            x = xg(i)
            z = zg(j)
            do k = 0,m-1
               i1 = int(x/dxeq)
               j1 = int(z/dzeq)

               i1 = int(x/dxeq)
               wx0 = ((i1+1)*dxeq-x)/dxeq
               wx1 = 1.-wx0

               j1 = int(z/dzeq)
               wz0 = ((j1+1)*dzeq-z)/dzeq
               wz1 = 1-wz0

               bp = wx0*wz0*b0(i1,j1)+wx0*wz1*b0(i1,j1+1) &
                 +wx1*wz0*b0(i1+1,j1)+wx1*wz1*b0(i1+1,j1+1) 
               bxp = wx0*wz0*b0x(i1,j1)+wx0*wz1*b0x(i1,j1+1) &
                 +wx1*wz0*b0x(i1+1,j1)+wx1*wz1*b0x(i1+1,j1+1) 
               bzp = wx0*wz0*b0z(i1,j1)+wx0*wz1*b0z(i1,j1+1) &
                 +wx1*wz0*b0z(i1+1,j1)+wx1*wz1*b0z(i1+1,j1+1) 
               btorp = wx0*wz0*b0zeta(i1,j1)+wx0*wz1*b0zeta(i1,j1+1) &
                 +wx1*wz0*b0zeta(i1+1,j1)+wx1*wz1*b0zeta(i1+1,j1+1) 
               radiusp = xctr-xdim/2+x               
               x = x+dzt*radiusp*bxp/btorp
               x = min(x,lx)
               x = max(x,0.)               
               z = z+dzt*radiusp*bzp/btorp
               z = min(z,lz)
               z = max(z,0.)
            end do
            ileft(i,j) = int(x/dx)
            jleft(i,j) = int(z/dz)
            xbackw(i,j)=x
            zbackw(i,j)=z
 !           write(*,*)i,ileft(i,j),j,jleft(i,j)
         end do
      end do

      return
    end subroutine weight
    
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



 
       subroutine ComputeInitialGuess(ksp,init_guess,ctx,ierr)
       use petscksp
       implicit none
       PetscErrorCode  ierr
       KSP ksp
       PetscInt ctx(*)
       Vec init_guess
       PetscScalar  h

       h=0.0
       PetscCall(VecSet(init_guess,h,ierr))
       end subroutine
       !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ComputeMatrix(ksp,AA,BB,dummy,ierr)
      use petscksp
      use gemx_com
      use equil
       implicit none
       PetscErrorCode  ierr
       KSP ksp
       Mat AA,BB
       integer dummy(*)
       DM dm
       integer :: ii,jj

      PetscInt    i,j,mx,my,xm
      PetscInt    ym,xs,ys,i1,i5
      PetscScalar  v(5),Hx,Hy
      PetscScalar  Hx2,Hy2,tmp_r,a_value
      MatStencil   row(4),col(4,5)

      i1 = 1
      i5 = 5
      a_value = 0.5
      PetscCall(KSPGetDM(ksp,dm,ierr))
      PetscCall(DMDAGetInfo(dm,PETSC_NULL_INTEGER,mx,my,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr))

      Hx =dx! (Rgrid(imx)-Rgrid(0)) / real(imx)
      Hy =dz!(Zgrid(jmx)-Zgrid(0)) / real(jmx)
      ! rho_i = sqrt(2*t0i(mid_i,mid_j)/mims(1))*mims(1)/(q(1)*b0(mid_i,mid_j))
      ! write(*,*) rho_i
      Hx2 = Hx**2
      Hy2 = Hy**2
      PetscCall(DMDAGetCorners(dm,xs,ys,PETSC_NULL_INTEGER,xm,ym,PETSC_NULL_INTEGER,ierr))
      do 10,j=ys,ys+ym-1
         do 20, i=xs,xs+xm-1
          row(MatStencil_i) = i
          row(MatStencil_j) = j
          if (mask(i,j) <0.99) then
             v(1)=(c2_over_vA2(i,j) + PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_i(i,j)**2))*(-2.0/Hx2-2.0/Hy2)
             PetscCall(MatSetValuesStencil(BB,i1,row,i1,row,v(1),INSERT_VALUES,ierr))
          else
             if (j > 0) then
                 if (j == jmx) then
                    v(1) = (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_i(i,j)**2))/Hy2-1.0/(2.0*Hy2)*( c2_over_vA2(i,j)- c2_over_vA2(i,j-1))
                    v(1) = v(1) + PADE*(mu0*e*e*rho_i(i,j)**2)*(xn0e(i,j)/t0e(i,j) - xn0e(i,j-1)/t0e(i,j-1))/Hy2
                 else
                    v(1) = (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_i(i,j)**2))/Hy2-1.0/(4.0*Hy2)*( c2_over_vA2(i,j+1)- c2_over_vA2(i,j-1))
                    v(1) = v(1) + PADE*(mu0*e*e*rho_i(i,j)**2)*(xn0e(i,j+1)/t0e(i,j+1) - xn0e(i,j-1)/t0e(i,j-1))/(2*Hy2)
                 end if
             end if
             col(MatStencil_i, 1) = i
             col(MatStencil_j, 1) = j - 1    

             if (i > 0) then
                 if (i == imx) then
                    v(2) =  (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_i(i,j)**2))/Hx2-1.0/(2.0*Hx2)*( c2_over_vA2(i,j)- c2_over_vA2(i-1,j))
                    v(2) = v(2) + PADE*(mu0*e*e*rho_i(i,j)**2)*(xn0e(i,j)/t0e(i,j) - xn0e(i-1,j)/t0e(i-1,j))/Hx2
                 else
                    v(2) =  (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_i(i,j)**2))/Hx2-1.0/(4.0*Hx2)*( c2_over_vA2(i+1,j)- c2_over_vA2(i-1,j))
                    v(2) = v(2) + PADE*(mu0*e*e*rho_i(i,j)**2)*(xn0e(i+1,j)/t0e(i+1,j) - xn0e(i-1,j)/t0e(i-1,j))/(2*Hx2)
                 end if
             end if
             col(MatStencil_i, 2) = i - 1
             col(MatStencil_j, 2) = j

             v(3) = -2.0*(c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_i(i,j)**2))/Hx2 - 2.0*(c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_i(i,j)**2))/Hy2
             col(MatStencil_i, 3) = i
             col(MatStencil_j, 3) = j

 !            write(*,*)v(3), xn0e(i,j)*mu0*e/t0e(i,j)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Boltzmann e!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
             if(iBoltzmann/=0) then
                v(3)=v(3)-xn0e(i,j)*mu0*e*e/t0e(i,j) + PADE*(mu0*e*e*rho_i(i,j)**2)*((xn0e(i-1,j)/t0e(i-1,j) + xn0e(i+1,j)/t0e(i+1,j) - 2*xn0e(i,j)/t0e(i,j))/Hx2 + &
                (xn0e(i,j-1)/t0e(i,j-1) + xn0e(i,j+1)/t0e(i,j+1) - 2*xn0e(i,j)/t0e(i,j))/Hy2)
             end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             

             if (i < imx) then
                 if (i == 0) then
                    v(4) =  (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_i(i,j)**2))/Hx2+1.0/(2.0*Hx2)*( c2_over_vA2(i+1,j)- c2_over_vA2(i,j))
                    v(4) = v(4) + PADE*(mu0*e*e*rho_i(i,j)**2)*(xn0e(i+1,j)/t0e(i+1,j) - xn0e(i,j)/t0e(i,j))/Hx2
                 else
                    v(4) =  (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_i(i,j)**2))/Hx2+1.0/(4.0*Hx2)*( c2_over_vA2(i+1,j)- c2_over_vA2(i-1,j))
                    v(4) = v(4) + PADE*(mu0*e*e*rho_i(i,j)**2)*(xn0e(i+1,j)/t0e(i+1,j) - xn0e(i-1,j)/t0e(i-1,j))/(2*Hx2)
                 end if
             end if
             col(MatStencil_i, 4) = i + 1
             col(MatStencil_j, 4) = j

             if (j < jmx) then
                 if (j == 0) then
                    v(5) =  (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_i(i,j)**2))/Hy2+1.0/(2.0*Hy2)*( c2_over_vA2(i,j+1)- c2_over_vA2(i,j))
                    v(5) = v(5) + PADE*(mu0*e*e*rho_i(i,j)**2)*(xn0e(i,j+1)/t0e(i,j+1) - xn0e(i,j)/t0e(i,j))/Hy2
                 else
                    v(5) =  (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_i(i,j)**2))/Hy2+1.0/(4.0*Hy2)*( c2_over_vA2(i,j+1)- c2_over_vA2(i,j-1))
                    v(5) = v(5) + PADE*(mu0*e*e*rho_i(i,j)**2)*(xn0e(i,j+1)/t0e(i,j+1) - xn0e(i,j-1)/t0e(i,j-1))/(2*Hy2)
                 end if
             end if
             col(MatStencil_i, 5) = i
             col(MatStencil_j, 5) = j + 1
             PetscCall(MatSetValuesStencil(BB, i1, row, i5, col, v, INSERT_VALUES, ierr))
          endif

    
20       continue
10    continue
      PetscCall(MatAssemblyBegin(BB,MAT_FINAL_ASSEMBLY,ierr))
      PetscCall(MatAssemblyEnd(BB,MAT_FINAL_ASSEMBLY,ierr))
      if (AA .ne. BB) then
         PetscCall(MatAssemblyBegin(AA,MAT_FINAL_ASSEMBLY,ierr))
         PetscCall(MatAssemblyEnd(AA,MAT_FINAL_ASSEMBLY,ierr))
      endif
!      PetscCall(MatView(A,PETSC_VIEWER_STDOUT_WORLD,petsc_ierr))
!      PetscCall(MatView(B,PETSC_VIEWER_STDOUT_WORLD,petsc_ierr))
    end subroutine

    
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc


       subroutine ComputeRHS(ksp,bbb,k,ierr)
       use petscksp
       use gemx_com
       use equil
       implicit none
       integer::k,ii,jj,iflag

       PetscErrorCode  ierr
       PetscScalar, POINTER ::b_array(:)

       KSP ksp
       Vec bbb
       PetscScalar  h,Hx,Hy
       PetscInt  mx,my,i,j,xs,xm,ys,ym,vec_start,vec_end
       DM dm
       PetscInt idx
       PetscScalar tmp_value,a_value,tmp_r

       tmp_value=0

  
       PetscCall(KSPGetDM(ksp,dm,ierr))
       PetscCall(DMDAGetInfo(dm,PETSC_NULL_INTEGER,mx,my,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr))
       PetscCallA(VecGetOwnershipRange(bbb,vec_start,vec_end,ierr))
       PetscCall(DMDAGetCorners(dm,xs,ys,PETSC_NULL_INTEGER,xm,ym,PETSC_NULL_INTEGER,ierr))



    idx=vec_start-1
       do 10,j=ys,ys+ym-1
         do 20, i=xs,xs+xm-1
            idx=idx+1           
             if (mask(i,j) <0.99) then
                tmp_value = 0
             else
               if (weightscheme == 0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!2D ni noly now!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if(i3D==0)then
                   if (iBoltzmann==0) then
                      tmp_value = denes(i,j,k)-q(1)*mu0*(den2d2(i,j)-xn0i(i,j))
                   else if (eAdiabatic/=0) then
                     tmp_value = -q(1)*mu0*(den2d2(i,j)-xn0i(i,j)) - (xn0e(i,j)*mu0*e*e/t0e(i,j))*phiavg(i,j) +  &
                              PADE*(mu0*q(1)*rho_i(i,j)**2)*(((den2d2(i-1,j)-xn0i(i-1,j))+(den2d2(i+1,j)-xn0i(i+1,j))-2*(den2d2(i,j)-xn0i(i,j)))/dx**2 + &
                               ((den2d2(i,j-1)-xn0i(i,j-1))+(den2d2(i,j+1)-xn0i(i,j+1))-2*(den2d2(i,j)-xn0i(i,j)))/dz**2) + PADE*(e*e*mu0*rho_i(i,j)**2)*(phiavg(i,j)*(((xn0e(i-1,j)/t0e(i-1,j)) + &
                               (xn0e(i+1,j)/t0e(i+1,j))-2*(xn0e(i,j)/t0e(i,j)))/dx**2 + ((xn0e(i,j-1)/t0e(i,j-1))+(xn0e(i,j+1)/t0e(i,j+1))-2*(xn0e(i,j)/t0e(i,j)))/dz**2) + & 
                               (xn0e(i,j)/t0e(i,j))*((phiavg(i-1,j) + phiavg(i+1,j) -2*phiavg(i,j))/dx**2 + (phiavg(i,j-1)+phiavg(i,j+1) -2*phiavg(i,j))/dz**2) + &
                               (((xn0e(i+1,j)/t0e(i+1,j))-(xn0e(i-1,j)/t0e(i-1,j)))*(phiavg(i+1,j)-phiavg(i-1,j))/(2*dx**2) + &
                               ((xn0e(i,j+1)/t0e(i,j+1))-(xn0e(i,j-1)/t0e(i,j-1)))*(phiavg(i,j+1)-phiavg(i,j-1))/(2*dz**2)))                  
                   else 
                     tmp_value = -q(1)*mu0*(den2d2(i,j)-xn0i(i,j))
                   end if
                
                   !                 tmp_value = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3D case!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                else

                if (iBoltzmann==0) then
                tmp_value = denes(i,j,k)-q(1)*mu0*(den(2,i,j,k)-xn0i(i,j))
                        
                else if (eAdiabatic/=0) then
                  tmp_value = -q(1)*mu0*(den(2,i,j,k)-xn0i(i,j)) - (xn0e(i,j)*mu0*e*e/t0e(i,j))*phiavg(i,j) + &
                  PADE*(mu0*q(1)*rho_i(i,j)**2)*(((den(2,i-1,j,k)-xn0i(i-1,j))+(den(2,i+1,j,k)-xn0i(i+1,j))-2*(den(2,i,j,k)-xn0i(i,j)))/dx**2 + &
                  ((den(2,i,j-1,k)-xn0i(i,j-1))+(den(2,i,j+1,k)-xn0i(i,j+1))-2*(den(2,i,j,k)-xn0i(i,j)))/dz**2) + PADE*(e*e*mu0*rho_i(i,j)**2)*(phiavg(i,j)*(((xn0e(i-1,j)/t0e(i-1,j)) + &
                  (xn0e(i+1,j)/t0e(i+1,j))-2*(xn0e(i,j)/t0e(i,j)))/dx**2 + ((xn0e(i,j-1)/t0e(i,j-1))+(xn0e(i,j+1)/t0e(i,j+1))-2*(xn0e(i,j)/t0e(i,j)))/dz**2) + & 
                  (xn0e(i,j)/t0e(i,j))*((phiavg(i-1,j) + phiavg(i+1,j) -2*phiavg(i,j))/dx**2 + (phiavg(i,j-1)+phiavg(i,j+1) -2*phiavg(i,j))/dz**2) + &
                  (((xn0e(i+1,j)/t0e(i+1,j))-(xn0e(i-1,j)/t0e(i-1,j)))*(phiavg(i+1,j)-phiavg(i-1,j))/(2*dx**2) + &
                  ((xn0e(i,j+1)/t0e(i,j+1))-(xn0e(i,j-1)/t0e(i,j-1)))*(phiavg(i,j+1)-phiavg(i,j-1))/(2*dz**2)))
                else 
                  tmp_value = -q(1)*mu0*(den(2,i,j,k)-xn0i(i,j))
                endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
                endif
               else
                  if(i3D==0)then
                     if (iBoltzmann==0) then
                        tmp_value = denes(i,j,k)-q(1)*mu0*(den2d2(i,j))
                     else if (eAdiabatic/=0) then
                        tmp_value = -q(1)*mu0*(den2d2(i,j)) - (xn0e(i,j)*mu0*e*e/t0e(i,j))*phiavg(i,j) +  &
                        PADE*(mu0*q(1)*rho_i(i,j)**2)*((den2d2(i-1,j)+den2d2(i+1,j)-2*den2d2(i,j))/dx**2 + &
                        (den2d2(i,j-1)+den2d2(i,j+1)-2*den2d2(i,j))/dz**2) + PADE*(e*e*mu0*rho_i(i,j)**2)*(phiavg(i,j)*(((xn0e(i-1,j)/t0e(i-1,j)) + &
                        (xn0e(i+1,j)/t0e(i+1,j))-2*(xn0e(i,j)/t0e(i,j)))/dx**2 + ((xn0e(i,j-1)/t0e(i,j-1))+(xn0e(i,j+1)/t0e(i,j+1))-2*(xn0e(i,j)/t0e(i,j)))/dz**2) + & 
                        (xn0e(i,j)/t0e(i,j))*((phiavg(i-1,j) + phiavg(i+1,j) -2*phiavg(i,j))/dx**2 + (phiavg(i,j-1)+phiavg(i,j+1) -2*phiavg(i,j))/dz**2) + &
                        (((xn0e(i+1,j)/t0e(i+1,j))-(xn0e(i-1,j)/t0e(i-1,j)))*(phiavg(i+1,j)-phiavg(i-1,j))/(2*dx**2) + &
                        ((xn0e(i,j+1)/t0e(i,j+1))-(xn0e(i,j-1)/t0e(i,j-1)))*(phiavg(i,j+1)-phiavg(i,j-1))/(2*dz**2)))                     
                     else 
                        tmp_value = -q(1)*mu0*(den2d2(i,j))
                     end if
                  
                  else

                  if (iBoltzmann==0) then
                     tmp_value = denes(i,j,k)-q(1)*mu0*(den(2,i,j,k))
                          
                  else if (eAdiabatic/=0) then
                     tmp_value = -q(1)*mu0*(den(2,i,j,k)) - (xn0e(i,j)*mu0*e*e/t0e(i,j))*phiavg(i,j) + &
                     PADE*(mu0*q(1)*rho_i(i,j)**2)*((den(2,i-1,j,k)+den(2,i+1,j,k)-2*den(2,i,j,k))/dx**2 + &
                     (den(2,i,j-1,k)+den(2,i,j+1,k)-2*den(2,i,j,k))/dz**2) + PADE*(e*e*mu0*rho_i(i,j)**2)*(phiavg(i,j)*(((xn0e(i-1,j)/t0e(i-1,j)) + &
                     (xn0e(i+1,j)/t0e(i+1,j))-2*(xn0e(i,j)/t0e(i,j)))/dx**2 + ((xn0e(i,j-1)/t0e(i,j-1))+(xn0e(i,j+1)/t0e(i,j+1))-2*(xn0e(i,j)/t0e(i,j)))/dz**2) + & 
                     (xn0e(i,j)/t0e(i,j))*((phiavg(i-1,j) + phiavg(i+1,j) -2*phiavg(i,j))/dx**2 + (phiavg(i,j-1)+phiavg(i,j+1) -2*phiavg(i,j))/dz**2) + &
                     (((xn0e(i+1,j)/t0e(i+1,j))-(xn0e(i-1,j)/t0e(i-1,j)))*(phiavg(i+1,j)-phiavg(i-1,j))/(2*dx**2) + &
                     ((xn0e(i,j+1)/t0e(i,j+1))-(xn0e(i,j-1)/t0e(i,j-1)))*(phiavg(i,j+1)-phiavg(i,j-1))/(2*dz**2)))

                     ! write(*,*) -q(1)*mu0*(den(2,i,j,k)), -(xn0e(i,j)*mu0*e*e/t0e(i,j))*phiavg(i,j), PADE*(mu0*q(1)*rho_i(i,j)**2)*((den(2,i-1,j,k)+den(2,i+1,j,k)-2*den(2,i,j,k))/dx**2 + (den(2,i,j-1,k)+den(2,i,j+1,k)-2*den(2,i,j,k))/dz**2), PADE*(e*e*mu0*rho_i(i,j)**2)*(phiavg(i,j)*(((xn0e(i-1,j)/t0e(i-1,j)) + &
                     ! (xn0e(i+1,j)/t0e(i+1,j))-2*(xn0e(i,j)/t0e(i,j)))/dx**2 + ((xn0e(i,j-1)/t0e(i,j-1))+(xn0e(i,j+1)/t0e(i,j+1))-2*(xn0e(i,j)/t0e(i,j)))/dz**2) + & 
                     ! (xn0e(i,j)/t0e(i,j))*((phiavg(i-1,j) + phiavg(i+1,j) -2*phiavg(i,j))/dx**2 + (phiavg(i,j-1)+phiavg(i,j+1) -2*phiavg(i,j))/dz**2) + &
                     ! (((xn0e(i+1,j)/t0e(i+1,j))-(xn0e(i-1,j)/t0e(i-1,j)))*(phiavg(i+1,j)-phiavg(i-1,j))/(2*dx**2) + &
                     ! ((xn0e(i,j+1)/t0e(i,j+1))-(xn0e(i,j-1)/t0e(i,j-1)))*(phiavg(i,j+1)-phiavg(i,j-1))/(2*dz**2)))
                  else 
                     tmp_value = -q(1)*mu0*(den(2,i,j,k))
                  endif      
                  endif
               end if
                
             end if
             PetscCall(VecSetValues(bbb,1,idx, tmp_value, INSERT_VALUES, ierr))
       
20       continue
10    continue
          
          
        
             PetscCall(VecAssemblyBegin(bbb,ierr))
             PetscCall(VecAssemblyEnd(bbb,ierr))

       
       end subroutine


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine smooth(matrix,mk)
       use gemx_com
       use equil

       IMPLICIT NONE
       real,dimension(0:imx,0:jmx,0:kmx)::matrix,temp
       integer::i,j,k,mk

       if(mk==3)then
       do k=0,kmx
          do i=2,imx-2
             do j=2,jmx-2
                if(mask3(i,j)<2.99)then
                else
                   temp(i,j,k)=(matrix(i,j,k)+matrix(i+1,j,k)+matrix(i,j+1,k)+matrix(i-1,j,k)+matrix(i,j-1,k))*0.2
                endif
                end do
          end do
       end do
     elseif(mk==2)then
       do k=0,kmx
          do i=2,imx-2
             do j=2,jmx-2
                if(mask2(i,j)<1.99)then
                else
                   temp(i,j,k)=(matrix(i,j,k)+matrix(i+1,j,k)+matrix(i,j+1,k)+matrix(i-1,j,k)+matrix(i,j-1,k))*0.2
                endif
                end do
          end do
       end do
    end if
    
       matrix=temp
     end subroutine smooth
     
     !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     subroutine integ(iflag)
       
       use gemx_com
       use equil

       IMPLICIT NONE
       integer::i,j,k,ip,m,iflag,itemp
       real::wx0,wx1,wzeta0,wzeta1,wy0,wy1,x,z,zeta,R_major_over_R,R_major_over_R1, ave_den,avex

       start_integ_tm = MPI_WTIME()
       itemp=2
       den(iflag,:,:,:)=0
       upar=0
       !$acc parallel loop gang vector
       do m=1,mm(1)
         x=x3(m)
         i = int(x/dxeq)
         wx0 = (i+1)-x/dxeq
         wx1 = 1.-wx0

         R_major_over_R=xctr/(xctr-xdim/2+i*dx)
         R_major_over_R1=xctr/(xctr-xdim/2+(i+1)*dx)

         z = z3(m)
         j = int(z/dzeq)
         wy0 = (j+1)-z/dzeq
         wy1 = 1.-wy0

         zeta=modulo(zeta3(m),2*pi)
         k=int(zeta/dzeta)
         wzeta0=(k+1)-zeta/dzeta
         wzeta1=1.-wzeta0

         ! if (weightscheme == 0) then

         !$acc atomic update 
         den(iflag,i,j,k)=den(iflag,i,j,k)+gw(m)*w3(m)*wx0*wy0*wzeta0*R_major_over_R
         !$acc atomic update
         den(iflag,i+1,j,k)=den(iflag,i+1,j,k)+gw(m)*w3(m)*wx1*wy0*wzeta0*R_major_over_R1
         !$acc atomic update
         den(iflag,i,j+1,k)=den(iflag,i,j+1,k)+gw(m)*w3(m)*wx0*wy1*wzeta0*R_major_over_R
         !$acc atomic update
         den(iflag,i+1,j+1,k)=den(iflag,i+1,j+1,k)+gw(m)*w3(m)*wx1*wy1*wzeta0*R_major_over_R1
         !$acc atomic update 
         upar(i,j,k)=upar(i,j,k)+u3(m)*gw(m)*w3(m)*wx0*wy0*wzeta0*R_major_over_R
         !$acc atomic update
         upar(i+1,j,k)=upar(i+1,j,k)+u3(m)*gw(m)*w3(m)*wx1*wy0*wzeta0*R_major_over_R1
         !$acc atomic update
         upar(i,j+1,k)=upar(i,j+1,k)+u3(m)*gw(m)*w3(m)*wx0*wy1*wzeta0*R_major_over_R
         !$acc atomic update
         upar(i+1,j+1,k)=upar(i+1,j+1,k)+u3(m)*gw(m)*w3(m)*wx1*wy1*wzeta0*R_major_over_R1
         
         if(k/=kmx)then
            !$acc atomic update
            den(iflag,i,j,k+1)=den(iflag,i,j,k+1)+gw(m)*w3(m)*wx0*wy0*wzeta1*R_major_over_R
            !$acc atomic update
            den(iflag,i+1,j,k+1)=den(iflag,i+1,j,k+1)+gw(m)*w3(m)*wx1*wy0*wzeta1*R_major_over_R1
            !$acc atomic update
            den(iflag,i,j+1,k+1)=den(iflag,i,j+1,k+1)+gw(m)*w3(m)*wx0*wy1*wzeta1*R_major_over_R
            !$acc atomic update
            den(iflag,i+1,j+1,k+1)=den(iflag,i+1,j+1,k+1)+gw(m)*w3(m)*wx1*wy1*wzeta1*R_major_over_R1
            !$acc atomic update
            upar(i,j,k+1)=upar(i,j,k+1)+u3(m)*gw(m)*w3(m)*wx0*wy0*wzeta1*R_major_over_R
            !$acc atomic update
            upar(i+1,j,k+1)=upar(i+1,j,k+1)+u3(m)*gw(m)*w3(m)*wx1*wy0*wzeta1*R_major_over_R1
            !$acc atomic update
            upar(i,j+1,k+1)=upar(i,j+1,k+1)+u3(m)*gw(m)*w3(m)*wx0*wy1*wzeta1*R_major_over_R
            !$acc atomic update
            upar(i+1,j+1,k+1)=upar(i+1,j+1,k+1)+u3(m)*gw(m)*w3(m)*wx1*wy1*wzeta1*R_major_over_R1
            
         else
            !$acc atomic update
            den(iflag,i,j,0)=den(iflag,i,j,0)+gw(m)*w3(m)*wx0*wy0*wzeta1*R_major_over_R
            !$acc atomic update
            den(iflag,i+1,j,0)=den(iflag,i+1,j,0)+gw(m)*w3(m)*wx1*wy0*wzeta1*R_major_over_R1
            !$acc atomic update
            den(iflag,i,j+1,0)=den(iflag,i,j+1,0)+gw(m)*w3(m)*wx0*wy1*wzeta1*R_major_over_R
            !$acc atomic update
            den(iflag,i+1,j+1,0)=den(iflag,i+1,j+1,0)+gw(m)*w3(m)*wx1*wy1*wzeta1*R_major_over_R1
            !$acc atomic update
            upar(i,j,0)=upar(i,j,0)+u3(m)*gw(m)*w3(m)*wx0*wy0*wzeta1*R_major_over_R
            !$acc atomic update
            upar(i+1,j,0)=upar(i+1,j,0)+u3(m)*gw(m)*w3(m)*wx1*wy0*wzeta1*R_major_over_R1
            !$acc atomic update
            upar(i,j+1,0)=upar(i,j+1,0)+u3(m)*gw(m)*w3(m)*wx0*wy1*wzeta1*R_major_over_R
            !$acc atomic update
            upar(i+1,j+1,0)=upar(i+1,j+1,0)+u3(m)*gw(m)*w3(m)*wx1*wy1*wzeta1*R_major_over_R1

         end if
      end do
!!         !$acc wait

      call MPI_Allreduce(MPI_IN_PLACE, den(iflag,:,:,:), (imx+1)*(jmx+1)*(kmx+1),MPI_Real8, MPI_SUM, MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(MPI_IN_PLACE, upar(:,:,:), (imx+1)*(jmx+1)*(kmx+1),MPI_Real8, MPI_SUM, MPI_COMM_WORLD,ierr)


      den(2,:,:,:) = den(iflag,:,:,:) 
      
         den2d2=0
         do k=0,kmx
            den2d2=den2d2+den(iflag,:,:,k)
            if (i3D==0 .and. k /= 0) upar(:,:,0)=upar(:,:,0)+upar(:,:,k)
         end do
         
         den2d2=den2d2/(kmx+1)
         if(i3D==0)then
            upar(:,:,0)=upar(:,:,0)/(kmx+1)
            do k=1,kmx
               upar(:,:,k)=upar(:,:,0)
            end do
         end if
         
         if(iflag==2)then
             dden2d=den2d2-den2d1
             den2d1=den2d2
 !        den_pre=den(2,:,:,:)
          end if
          
     end_integ_tm = MPI_WTIME()
     integ_tm = integ_tm + end_integ_tm - start_integ_tm 

    end subroutine integ
    
    
    
!     !!!!!!!!!!!!!!!!!!!!!!!!! CALDER Flux Average SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fluxavg(input,output)
   use gemx_com
   use equil

   implicit none 

   ! Input
   real(8), dimension(0:imx,0:jmx,0:kmx) :: input !3D array to be flux averaged
   ! Outputs
   real(8), dimension(0:imx, 0:jmx) :: output ! 2D interpolated array

   ! Local variables
   real(8), dimension(:), allocatable :: phiavg1d, psi1d, psi1d_private,phiavg1d_private
   integer :: gi, xix, yjy, miw, psi_zero, store, k, line_small, line_marker, priv_mark
   integer :: i, j
   real(8) :: weightinput,weightinput3D, phiavggi, psival, wmx0, wmx1, psi_private_min
   allocate(phiavg1d(0:99),psi1d(0:99),psi1d_private(0:99),phiavg1d_private(0:99))
   
   !!!!REMEMBER TO TURN OFF
!    do i = 0, imx
!         do j=0,jmx     
!                 do k=0,kmx
!                    input(i,j,k) = 100*psi_p(i,j)
!               end do
!        end do
!   end do
!!!!!!!!!

   phiavg1d         = 0
   phiavg1d_private = 0
   psi1d            = 0
   psi1d_private    = 0
   
   line_small = 0
   psi_zero   = 1 !Start at one, sets zeroeth index later

   !Revised Computation!!!!!!!!!!!!!!!!
   line_marker = 0

   do gi = 0, 99
      weightinput = 0.0d0
      phiavggi    = 0.0d0 
      store       = 0
      priv_mark   = 0
      do line = line_marker, num_lines
        if (gindex(line) == gi) then
              weightinput = (weight00(line)*input(iarray(line), jarray(line),0) + &
                          weight10(line)*input(iarray(line)+1, jarray(line),0) + &
                          weight01(line)*input(iarray(line), jarray(line)+1,0) + &
                          weight11(line)*input(iarray(line)+1, jarray(line)+1,0))
              if (i3D == 0) then                  
                 phiavggi = phiavggi + (weightinput*jacobian(line))/deno(line)
              else
                 do k=1, kmx
                    weightinput3D = (weight00(line)*input(iarray(line), jarray(line),k) + &
                                   weight10(line)*input(iarray(line)+1, jarray(line),k) + &
                                   weight01(line)*input(iarray(line), jarray(line)+1,k) + &
                                   weight11(line)*input(iarray(line)+1, jarray(line)+1,k))
                    weightinput = weightinput3D + weightinput
                 enddo
                 phiavggi = phiavggi +(weightinput*jacobian(line))/(deno(line)*(kmx+1))
              end if
              store = line
           if (priv(line) == 0) then
              priv_mark = 1
           end if
        end if
        !Remove redundancy from closed loop integration process
        if (phiavggi /= 0) then
           if (gindex(line) /= gi) then
              if (i3D == 0) then
                 phiavggi = phiavggi - (weightinput*jacobian(line-1))/deno(line-1)
              else
                 phiavggi = phiavggi - (weightinput*jacobian(line-1))/(deno(line-1)*(kmx+1))
              end if
              line_marker = line
              exit !Leave contour loop after removing redundancy
           end if
        end if
     end do
     
     if (priv_mark /= 0) then
        phiavg1d(psi_zero) = phiavggi
        psi1d(psi_zero)    = psitab(store)
        psi_zero = psi_zero + 1
     else 
        phiavg1d_private(gi) = phiavggi
        if (line_small == 0) then
           psi_private_min = psitab(store)
        end if
      !   phiavg1d_private(line_small) = phiavggi
      !   psi1d_private(line_small) = psitab(store)
        line_small = line_small + 1
     end if
  end do
  
  !Initialize output to zero
  phiavg1d(0) = phiavg1d(1)
  output = 0.0
  !QUICK FIX
!   phiavg1d_private = 0
  !!
  !INTERPOLATION
  do xix = 0, nx
      do yjy = 0, nz
          psival = psi_p(xix,yjy)
          if (mask(xix,yjy) < 0.99) then 
              output(xix,yjy) = 0
          else
            !   if (yjy < 75 .and. xix < 150 .and. psival > 0.29 .and. psival<0.31) then !Private region under X-point
            !      miw  = int((psival-psi_private_min)/(psi1d(2)-psi1d(1)))
            !    !   wmx0 = ((miw+1)*(psi1d_private(2)-psi1d_private(1))-psival)/(psi1d(2)-psi1d(1))
            !      wmx0 = ((miw+1)*(psi1d(2)-psi1d(1))-psival)/(psi1d(2)-psi1d(1))
            !      wmx1 = 1.-wmx0
            !      output(xix,yjy) = wmx0*phiavg1d_private(miw) + wmx1*phiavg1d_private(miw+1)
            !   else 
                 miw  = int(psival/(psi1d(2)-psi1d(1)))
                 wmx0 = ((miw+1)*(psi1d(2)-psi1d(1))-psival)/(psi1d(2)-psi1d(1))
                 wmx1 = 1. - wmx0  
                 output(xix,yjy) = wmx0*phiavg1d(miw) + wmx1*phiavg1d(miw+1)
            !   end if
          end if
      enddo
   enddo
   ! output = 0.0

end subroutine fluxavg
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! E FIELD SUBROUTINE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine efieldcalc(phi_input)
   use gemx_com
   use equil

   implicit none 
   ! Input
   real(8), dimension(0:imx,0:jmx,0:kmx) :: phi_input !3D phi for calculation of E field

   ! Local variables
   integer :: i, j, k, kminus, kplus

   do k=0, kmx
      do i=2, imx-1
         do j=2, jmx-1
            ex(i,j,k)    = -(phi_input(i+1,j,k) - phi_input(i-1,j,k))/(2*(Rgrid(1)-Rgrid(0)))
            ez(i,j,k)    = -(phi_input(i,j+1,k) - phi_input(i,j-1,k))/(2*(Zgrid(1)-Zgrid(0)))
            if (k==0) then
               kminus = kmx
               ezeta(i,j,k) = -(phi_input(i,j,k+1) - phi_input(i,j,kminus))/(2*Rgrid(i)*(2*pi/(kmx+1)))
            else if (k==kmx) then
               kplus = 0
               ezeta(i,j,k) = -(phi_input(i,j,kplus) - phi_input(i,j,k-1))/(2*Rgrid(i)*(2*pi/(kmx+1)))
            else
               ezeta(i,j,k) = -(phi_input(i,j,k+1) - phi_input(i,j,k-1))/(2*Rgrid(i)*(2*pi/(kmx+1)))
            end if
         end do
      end do
   end do

end subroutine efieldcalc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Boltzmann-Poisson Electron Solver!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine boltzsolve(input_phi)
   use gemx_com
   use equil

   implicit none
   ! Input
   real(8), dimension(0:imx,0:jmx,0:kmx) :: input_phi !3D input phi array

   ! Local Variables
   integer :: i, j, k
   !Need to solve \phi_{k+1} = \phi{k} - \frac{R^k \Delta \phi^k}{R^{k+1}-R^k}


   !Define phi_{k+1}
   do i=0, imx
      do j=0, jmx
         do k=0, kmx
            if (input_phi(i,j,k) == 0) then
               phi_k(i,j,k) = 0.005+(ran2(iseed)*0.005)
            else
               phi_k(i,j,k) = 1.001*input_phi(i,j,k)
            end if
         end do
      end do
   end do
   
   do i=2,imx-1
      do j=2,jmx-1
         do k=0,kmx
            dphidr(i,j,k)   = Rgrid(i)*(input_phi(i+1,j,k)-input_phi(i-1,j,k))/(2*dx)
            dphi_kdr(i,j,k) = Rgrid(i)*(phi_k(i+1,j,k)-phi_k(i-1,j,k))/(2*dx)

            dphidz(i,j,k)   = (input_phi(i,j+1,k)-input_phi(i,j-1,k))/(2*dz)
            dphi_kdz(i,j,k) = (phi_k(i,j+1,k)-phi_k(i,j-1,k))/(2*dz)
         end do
      end do
   end do

   do i=2,imx-1
      do j=2,jmx-1
         do k=0,kmx
            d2phidr2(i,j,k)      = c2_over_vA2(i,j)*(dphidr(i+1,j,k)-dphidr(i-1,j,k))/(2*dx*Rgrid(i))
            d2phi_kdr2(i,j,k)    = c2_over_vA2(i,j)*(dphi_kdr(i+1,j,k)-dphi_kdr(i-1,j,k))/(2*dx*Rgrid(i))

            d2phidz2(i,j,k)      = c2_over_vA2(i,j)*(dphidz(i,j+1,k)-dphidz(i,j-1,k))/(2*dz)
            d2phi_kdz2(i,j,k)    = c2_over_vA2(i,j)*(dphi_kdz(i,j+1,k)-dphi_kdz(i,j-1,k))/(2*dz)
            
            OPPphi(i,j,k)  = (d2phidr2(i,j,k)+d2phidz2(i,j,k))
            OPPphik(i,j,k) = (d2phi_kdr2(i,j,k)+d2phi_kdz2(i,j,k))

            if (i3D == 0) then
               r_hand(i,j,k)  = OPPphi(i,j,k) + q(1)*mu0*den2d2(i,j) - e*mu0*xn0e(i,j)*exp((e/t0e(i,j))*input_phi(i,j,k))
               rk_hand(i,j,k) = OPPphik(i,j,k) + q(1)*mu0*den2d2(i,j) - e*mu0*xn0e(i,j)*exp((e/t0e(i,j))*phi_k(i,j,k))
            else 
               r_hand(i,j,k)  = OPPphi(i,j,k) + q(1)*mu0*den(2,i,j,k) - e*mu0*xn0e(i,j)*exp((e/t0e(i,j))*input_phi(i,j,k))
               rk_hand(i,j,k) = OPPphik(i,j,k) + q(1)*mu0*den(2,i,j,k) - e*mu0*xn0e(i,j)*exp((e/t0e(i,j))*phi_k(i,j,k))
            end if

            input_phi(i,j,k) = phi_k(i,j,k) - (input_phi(i,j,k)-phi_k(i,j,k))*r_hand(i,j,k)/(rk_hand(i,j,k)-r_hand(i,j,k))
            ! write(*,*) (r_hand(i,j,k)/l_hand(i,j,k))
            ! input_phi(i,j,k) = OPPphik(i,j,k)

            if (input_phi(i,j,k) == 0) then
               write(*,*) 'working'
            end if

            if (mask(i,j)<0.99) then
               input_phi(i,j,:) = 0
            end if
         end do
      end do
   end do

   ! if(MyId==0 .and. timestep == 5)then

   ! open(unit=11, file = 'testphi_boltzmann',status='unknown',action='write')
   ! do j=0,jmx
      
   !    write(11,*) input_phi(:,j,0)
   ! enddo
   
   ! close(11)

   ! end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Fourier Solve!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Toroidal Direction
subroutine fourier_modes(input_phi, modes)
   use gemx_com
   use equil

   implicit none
   #include "fftw3.f"

   ! Input
   real(8), dimension(0:imx, 0:jmx, 0:kmx) :: input_phi  ! 3D input phi array
   integer :: modes
   ! integer, dimension(:), intent(in) :: modes  ! Array of modes to keep

   ! Local Variables

   double complex, dimension(:), allocatable :: phi_hat!f_hat!,f_filtered
   real, dimension(:), allocatable :: f_filtered
   real, allocatable :: omega_r, gamma
   !!!
   integer :: i, j, k
   integer*8 :: plan_forward, plan_backward
   allocate(phi_hat(0:((kmx+1)/2 + 1)),f_filtered(0:kmx))
 
   
   do i = 0, imx
      do j = 0, jmx
         !Perform forward FFT
         call dfftw_plan_dft_r2c_1d(plan_forward,kmx+1,input_phi(i,j,:),phi_hat,FFTW_ESTIMATE)
         call dfftw_execute_dft_r2c(plan_forward, input_phi(i,j,:), phi_hat)
         call dfftw_destroy_plan(plan_forward)

         ! if (mod(timestep,10)==0) then
         !    if (j = jmidplane .and. i = grad_peak) then            
         !       ! Calculate real frequency and growth rate
         !       do k = 0, ((kmx+1)/2 + 1)
         !          ! Real frequency: omega_r = 2 * pi * k / L
         !          omega_r = real(k) / R(grad_peak)
   
         !          ! Growth rate: look at imaginary part of phi_hat(k)
         !          gamma = aimag(phi_hat(k))  ! Extract the imaginary part (growth rate)
   
         !          ! Print or store the frequency and growth rate
         !          if (k == 8) then  ! For the selected mode, print values
         !             print*, 'Mode ', k, ' - Frequency: ', omega_r, ' - Growth rate: ', gamma
         !          end if
         !       end do
         !    end if
         ! end if


         ! Zero out unwanted modes
         do k = 0, ((kmx+1)/2 + 1)
            phi_hat(k) = phi_hat(k)/(kmx+1)
            if (k /= modes) then !n toroidal modes
               phi_hat(k) = (0.0,0.0)
            end if
         end do

         ! Perform inverse FFT to reconstruct the input_phi
         call dfftw_plan_dft_c2r_1d(plan_backward, kmx+1, phi_hat, f_filtered, FFTW_ESTIMATE)
         call dfftw_execute_dft_c2r(plan_backward, phi_hat, f_filtered)
         call dfftw_destroy_plan(plan_backward)

         do k = 0, kmx
            input_phi(i,j,k) = f_filtered(k)
            ! if (timestep==10 .and. f_filtered(k)/=0) then
            !    write(*,*) input_phi(i,j,k)
            ! end if
         end do

      end do
   end do

   if(MyId==0 .and. mod(timestep,10)==0)then
      ! write(*,*) phi(:,:,outk)

   open(unit=11, file = 'testphi_fourier',status='unknown',action='write')
   do j=0,jmx
      
      write(11,*) input_phi(:,j,0)
   enddo
   
   close(11)

   end if

   

   ! Deallocate phi_hat if necessary
   deallocate(phi_hat,f_filtered)

end subroutine fourier_modes



! !Poloidal Direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine poloidal_fourier_filter(input_phi)
!    use gemx_com
!    use equil

!    implicit none
!    #include "fftw3.f"

!    ! Input
!    real(8), dimension(0:imx, 0:jmx, 0:kmx) :: input  ! 3D input phi array
! ! 
!    ! Local Variables
!    real(8), dimension(:), allocatable :: contour_phi

!    double complex, dimension(:), allocatable :: phi_hat!f_hat!,f_filtered
!    real, dimension(:), allocatable :: f_filtered
!    !!!
!    integer :: i, j, k
!    integer*8 :: plan_forward, plan_backward
!    allocate(phi_hat(0:((kmx+1)/2 + 1)),f_filtered(0:kmx))

!    line_marker = 0

!    do gi = 0, 101
!       !first reallocate contour phi for each contour group
!       weightinput = 0.0d0
!       phiavggi    = 0.0d0 
!       store       = 0
!       priv_mark   = 0

!       !Calculate length of contour for filtering
!       do line = linemarker, num_lines
!          if (gindex(line) /= gi) then
!             size = gindex(line) - linemarker
!             write(*,*) size
!             exit !Leave loop once contour length has been found
!          end if
!       end do

!       !Allocate contour, fourier transformed contour, and filtered contour arrays using contour length
!       allocate(contour_phi(0:size),contour_hat(0:(size/2)+1),contour_filtered(0:size))
      
!       !Uses jaco.dat to interpolate contour onto psi grid correctly
!       do line = line_marker, num_lines
!         if (gindex(line) == gi) then
!               weightinput = (weight00(line)*input(iarray(line), jarray(line),0) + &
!                           weight10(line)*input(iarray(line)+1, jarray(line),0) + &
!                           weight01(line)*input(iarray(line), jarray(line)+1,0) + &
!                           weight11(line)*input(iarray(line)+1, jarray(line)+1,0))
!              contour_phi(line_marker-line) = weightinput
!         elseif (gindex(line) /= gi) then
!          line_marker = line 
!          exit
!         end if 
!       end do

!       !Now Fourier Filter
!       call dfftw_plan_dft_r2c_1d(plan_forward,size,contour_phi(:),contour_hat,FFTW_ESTIMATE)
!       call dfftw_execute_dft_r2c(plan_forward,contour_phi(:),contour_hat)
!       call dfftw_destroy_plan(plan_forward)

!       ! Zero out unwanted modes
!       do point = 0, (size/2 + 1)
!          contour_hat(point) = contour_hat(point)/size !normalize
!          if (point/=0 .and. point/=1) then
!             contour_hat(point) = 0.0
!          end if
!       end do

!       !Perform inverse FFT to reconstruct the filtered contour line
!       call dfftw_plan_dft_c2r_1d(plan_backward,size,contour_hat,contour_filtered,FFTW_ESTIMATE)
!       call dfftw_execute_dft_c2r(plan_backward,contour_hat,contour_filtered)
!       call dfftw_destroy_plan(plan_backward)

!       !Now replace contour with filtered contour on phi grid
      

!    end do



! end subroutine poloidal_fourier_filter
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
! !!!!GRID FILTERING
subroutine binomial_filter(input_phi)
   use gemx_com
   use equil

   implicit none

   ! Input
   real(8), dimension(0:imx, 0:jmx, 0:kmx) :: input_phi  ! 3D input phi array

   ! Local variables
   integer :: i, j, k, kmin, kmax
   real(8), dimension(:,:,:), allocatable :: Filter
   allocate(Filter(1:imx-1,1:jmx-1,0:kmx))

   !Binomial Filter
   do i = 1, imx-1
      do j = 1, jmx-1
         do k = 0, kmx
            if (k == 0) then
               kmin = kmx
            else
               kmin = k-1
            end if

            if (k == kmx) then
               kmax = 0
            else
               kmax = k+1
            end if

            Filter(i,j,k) = (1.0d0 / 64.0d0) * (input_phi(i-1,j-1,kmin) + 2*input_phi(i,j-1,kmin) + input_phi(i+1,j-1,kmin) + &
                           2*input_phi(i-1,j,kmin) + 4*input_phi(i,j,kmin) + 2*input_phi(i+1,j,kmin) + &
                           input_phi(i-1,j+1,kmin) + 2*input_phi(i,j+1,kmin) + input_phi(i+1,j+1,kmin) + &
                           2*input_phi(i-1,j-1,k) + 4*input_phi(i,j-1,k) + 2*input_phi(i+1,j-1,k) + &
                           4*input_phi(i-1,j,k) + 8*input_phi(i,j,k) + 4*input_phi(i+1,j,k) + &
                           2*input_phi(i-1,j+1,k) + 4*input_phi(i,j+1,k) + 2*input_phi(i+1,j+1,k) + &
                           input_phi(i-1,j-1,kmax) + 2*input_phi(i,j-1,kmax) + input_phi(i+1,j-1,kmax) + &
                           2*input_phi(i-1,j,kmax) + 4*input_phi(i,j,kmax) + 2*input_phi(i+1,j,kmax) + &
                           input_phi(i-1,j+1,kmax) + 2*input_phi(i,j+1,kmax) + input_phi(i+1,j+1,kmax))
         end do
      end do
   end do

   do i = 0, imx
      do j = 0, jmx
         do k = 0, kmx
            if (mask(i,j)<0.99) then
               input_phi(i,j,k) = 0
            else
               input_phi(i,j,k) = Filter(i,j,k)
            end if
         end do
      end do
   end do

end subroutine binomial_filter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine poloidal_filter_methods(input_phi)
   use gemx_com
   use equil

   implicit none

   ! Input
   real(8), dimension(0:imx, 0:jmx, 0:kmx) :: input_phi  ! 3D input phi array

   ! Local variables
   integer :: i, j, k, kmin, kmax
   real(8), dimension(:,:,:), allocatable :: Filter
   allocate(Filter(1:imx-1,1:jmx-1,0:kmx))

   !!!!!!!!!!!!!!!!Binomial Filter in the poloidal plane !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do k = 0, kmx
      do i =1, imx-1
         do j = 1, jmx-1
            ! write(*,*) input_phi(i,j,k)
            Filter(i,j,k) = (1.0d0/16.0d0) * (input_phi(i-1,j-1,k) + 2*input_phi(i,j-1,k) + input_phi(i+1,j-1,k) + &
                           2*input_phi(i-1,j,k) + 4*input_phi(i,j,k) + 2*input_phi(i+1,j,k) + &
                           input_phi(i-1,j+1,k) + 2*input_phi(i,j+1,k) + input_phi(i+1,j+1,k))
         end do
      end do
   end do

   do i = 0, imx
      do j = 0, jmx
         do k = 0, kmx
            if (mask(i,j)<0.99) then
               input_phi(i,j,k) = 0
            else
               ! write(*,*) Filter(i,j,k)
               input_phi(i,j,k) = Filter(i,j,k)
            end if
         end do
      end do
   end do   

end subroutine poloidal_filter_methods

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine ftcamp(input_phi,timestep)

!    use gemx_com
!    use equil

!    implicit none
!    #include "fftw3.f"

!    ! Input
!    real(8), dimension(0:imx, 0:jmx, 0:kmx) :: input_phi  ! 3D input phi array
!    integer :: timestep
!    ! integer, dimension(:), intent(in) :: modes  ! Array of modes to keep

!    ! Local Variables
!    double complex, dimension(:,:,:), allocatable :: phi_hat
!    real, allocatable :: omega_r, gamma
!    !!!
!    integer :: i, j, k
!    integer*8 :: plan_forward, plan_backward
!    ! allocate(phi_hat(0:((imx+1)/2 + 1),0:((jmx+1)/2 + 1),0:((kmx+1)/2 + 1)))
!    allocate(phi_hat(0:((imx+1)/2 + 1),0:jmx,0:kmx))
!    ! write(*,*) 'test placement'
!    !Perform forward FFT
!    call dfftw_plan_dft_r2c_3d(plan_forward, imx+1,jmx+1,kmx+1, input_phi, phi_hat, FFTW_FORWARD,FFTW_ESTIMATE)
!    write(*,*) 'test placement 2'
!    call dfftw_execute_dft_r2c(plan_forward, input_phi, phi_hat)
!    call dfftw_destroy_plan(plan_forward)

!    write(*,*) 'made it'

!    do i=0, imx
!       do j=0, jmx
!          do k=0, ((kmx+1)/2 + 1)
!             write(*,*) phi_hat(i,j,k)
!          end do
!       end do
!    end do
!    ! write(*,*) phi_hat

!    ! if (mod(timestep,10)==0) then
!          !    if (j = jmidplane .and. i = grad_peak) then            
!          !       ! Calculate real frequency and growth rate
!          !       do k = 0, ((kmx+1)/2 + 1)
!          !          ! Real frequency: omega_r = 2 * pi * k / L
!          !          omega_r = real(k) / R(grad_peak)
   
!          !          ! Growth rate: look at imaginary part of phi_hat(k)
!          !          gamma = aimag(phi_hat(k))  ! Extract the imaginary part (growth rate)
   
!          !          ! Print or store the frequency and growth rate
!          !          if (k == 8) then  ! For the selected mode, print values
!          !             print*, 'Mode ', k, ' - Frequency: ', omega_r, ' - Growth rate: ', gamma
!          !          end if
!          !       end do
!          !    end if
!          ! end if


! end subroutine ftcamp
