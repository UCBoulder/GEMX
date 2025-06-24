#include "gemx_c.hpp"
#include "pputil_c.hpp"
#include "fcnt.hpp"
#include "gemx_com_c.hpp"
#include "equil_c.hpp"
#include "mpi.h"
#include "ionPush_c.hpp"
#include "outd_c.hpp"
#include "MultiArraysC.hpp"

#include <cmath>
#include <chrono>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <petsc.h>

using namespace std;

int main() {
   void* kval; //used for compute rhs
   int status,mid_i,mid_j;
   int n,i,j,k,ip,m,outk,ix=135,jx=68;
   int iter; //Calder Edit

   //double random;
   double tmp;
   //double retVal;
   PetscInt is,js,iw,jw,idx;//,n_in_porcs;
   PetscInt one,three,vec_start,vec_end;
   //PetscErrorCode petsc_ierr;
   const PetscScalar* phi_array;
   KSP ksp;
   DM dm;
   //PetscObject  vec;
   Vec petsc_phi;//mpi_phi;
   //PetscViewer viewer;
   //VecScatter   ctx;

   //call init
   initialize_c_();


   outk=0; //(kmx+1)/2

   one = 1;
   three = 3;


   if(eBoltzmann == 0) {

      PETSC_COMM_WORLD = PETSC_COMM;

      PetscCall(PetscInitialize(nullptr, nullptr, nullptr, nullptr));  
      
      PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp));
      PetscCall(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, 
                             DMDA_STENCIL_STAR,imx+1,jmx+1,PETSC_DECIDE,PETSC_DECIDE,
                             one,one, nullptr, nullptr, &dm));
      PetscCall(DMSetFromOptions(dm));
      PetscCall(DMSetUp(dm));
      PetscCall(KSPSetDM(ksp,dm));
      PetscCall(KSPSetComputeInitialGuess(ksp,ComputeInitialGuess,nullptr)); 
      PetscCall(KSPSetComputeOperators(ksp,ComputeMatrix,nullptr));    	
      PetscCall(DMDAGetCorners(dm,&is,&js,nullptr,&iw,&jw,nullptr));
      PetscCall(KSPSetFromOptions(ksp));
      PetscCall(KSPSetUp(ksp)); 
   } 

   if(iget == 0) loadi_c_();
   integ_c_(1); //1st index of den 1-based, since index = 1,2 in ftn, index in c is 0,1 meaning 1 is same index in c as ftn. Funkalicious

   if(myid == 0) {
      //write a bunch of stuff to files
   }

   starttm=MPI_Wtime();
   upar.Clear();

   mid_i=imx/2;
   mid_j=jmx/2;
   mid_i=257;
   mid_j=257;

   tor_n = 1;

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!initialize perturbation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   for(int k = 0; k <= kmx; ++k){
      for(int i = 0; i <= imx; ++i){
         for(int j = 0; j <= jmx; ++j){
//                    call random_number(random)
//                    dene(i,j,k)=mask(i,j)*cos(2*pi*k/(kmx+1))*2*exp(-((i-mid_i)**2+(j-mid_j)**2)/(0.09*min(mid_i,mid_j))**2)




//                     !                    apar(i,j,k)=mask(i,j)*cos(tor_n*2*pi*k/(kmx+1))*2*exp(-((i-mid_i)**2+(j-mid_j)**2)/(0.09*min(mid_i,mid_j))**2)
//                     ! apar(i,j,k)=mask2(i,j)*2*exp(-((i-ix)**2+(j-jx)**2)/(0.04*257)**2)!*cos(tor_n*2*pi*k/(kmx+1))
//                     !                    apar(i,j,k)=mask2(i,j)*(exp(-(psi_p(i,j)-0.13)**2/0.01**2)-exp(-(psi_p(i,j)-0.16)**2/0.01**2))*cos(tor_n*2*pi*k/(kmx+1))
//                    ! if (i==mid_i .and. j==mid_j) then
//                    !    apar(i,j,k)=0
//                    ! else
//                       ! apar(i,j,k)=mask2(i,j)*(exp(-(psi_p(i,j)-0.1905)**2/0.01**2))*cos(tor_n*2*pi*k/(kmx+1))*(2*(j-mid_j)**2*dz**2/((j-mid_j)**2*dz**2+(i-mid_i)**2*dx**2)-1)!cos(2*pi*2*ATAN((j-mid_j)/(i-mid_i)))
//                    ! endif
                    
                    
// !                    apar(i,j,k)=mask(i,j)*(ran2(iseed)-0.5)
// !                    apar(i,j,k)=0
            apars(i,j,k) = 0; //could use Clear in future
            //apar(i,j,k) = 0;
            jpar(i,j,k)=0;
            dene(i,j,k)=0;//ran(-0.5); 
            phi(i,j,k)=0;//-j*0.01;
            ez(i,j,k)=0;//0.01/dz;
         }
      }
   }

   phiavg.Clear();

   get_jpar_(apar); 
   get_ne_c_(1);
   if(i3D==0){
      apar.Clear();
      dene.Clear();
      integ_c_(1);
   }


   if(myid==0){
      //write loads of stuff
   }
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end of init perturbation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            

   if(ifield_solver == 1) ncurr=1;
      
   start_total_tm = MPI_Wtime();
   for(timestep=ncurr; timestep<=nm; ++timestep) {
      for(int q = 0; q <=10006; ++q){
         if(ran2_c_(iseed)-0.5 > 0){
            rand_table[q]=1;
         } else {
            rand_table[q]=-1;
         }
      }
      tcurr = tcurr+dt;

   //	   accumulate(timestep-1,0)
   //	   ezamp()
   //	   gkps()
   //    field(timestep-1,0)


      if(ifield_solver == 1) {
         phi.Clear();
         phiavg.Clear();
         denes=dene;



         if(i3D != 0) {
            for(int k = myid*(kmx+1)/(numprocs); k < (myid+1)*(kmx+1)/(numprocs); ++k){
               for(iter = 0; iter<=iterations; ++iter){
                  fluxavg_c_(phi, phiavg);

                  if(eBoltzmann == 1){
                     BoltzSolve_c_(phi);
                     cout << "Boltz" << endl;
                  } else {
                     kval = (void*)k;
                     PetscCall(KSPSetComputeRHS(ksp,ComputeRHS,kval));
                     PetscCall(KSPSetComputeRHS(ksp,ComputeRHS,kval));
                     PetscCall(KSPSolve(ksp,nullptr,nullptr));
                     PetscCall(KSPGetSolution(ksp,&petsc_phi));
                     PetscCall(VecGetArrayRead(petsc_phi,&phi_array));
                     PetscCall(VecGetOwnershipRange(petsc_phi,&vec_start,&vec_end));
                     for(idx = 1; idx <= (vec_end-vec_start); ++idx) {
                        //cout << idx << endl;
                        i=((idx-1)%iw)+is;
                        j=(idx-1)/(iw)+js;
                        phi(i,j,k)=phi_array[idx];//*mask(i,j);  //right here officer
                     }
                     //cout << "after loop" << endl;
                     PetscCall(VecRestoreArrayRead(petsc_phi,&phi_array));
                  }
               }
            }
         ierr = MPI_Allreduce(MPI_IN_PLACE, phi.start(), (imx+1)*(jmx+1)*(kmx+1),MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      } else {
         k = 0;

         for(iter = 0; iter <= iterations; ++iter) {
            fluxavg_c_(phi, phiavg);
            if(eBoltzmann == 1) {
               BoltzSolve_c_(phi);
            } else {
               kval = (void*)k;
               PetscCall(KSPSetComputeRHS(ksp,ComputeRHS,kval));
               PetscCall(KSPSolve(ksp,NULL,NULL));
               PetscCall(KSPGetSolution(ksp,&petsc_phi));
               PetscCall(VecGetOwnershipRange(petsc_phi,&vec_start,&vec_end));
               PetscCall(VecGetArrayRead(petsc_phi,&phi_array));
               for(idx = 1; idx <=(vec_end-vec_start); ++idx) {
                  i=idx-1%(iw)+is;
                  j=(idx-1)/(iw)+js;
                  phi(i,j,k)=phi_array[idx];//*mask(i,j);
               }

               PetscCall(VecRestoreArrayRead(petsc_phi,&phi_array));

               ierr = MPI_Allreduce(MPI_IN_PLACE, phi.start(), (imx+1)*(jmx+1)*(kmx+1),MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

               for(int i = 0; i <= imx; ++i){
                  for(int j = 0; j <= jmx; ++j){
                     for(int k = 0; k <= kmx; ++k){
                        phi(i,j,k) = phi(i,j,0);
                     }
                  }
               }
            }
         }
      }
   }


      efieldcalc_c_(phi);



   //             if (Myid==0) then
   //                open(unit=11, file = 'testphis',status='unknown',action='write')
   //                do j=0,jmx
                     
   //                   write(11,*) phi(:,j,outk)-phi(:,j,outk+1)
   //                   enddo
   //                   close(11)
   //              endif

   //               open(unit=11, file = 'testphis1',status='unknown',action='write')
   //               do j=0,jmx
                     
   //                  write(11,*) phi(:,j,outk+1)
   //                  enddo
   //               close(11)
      

      get_apar_(-1);
      //smooth(apars,2);
      get_jpar_(apars);
      //smooth(jpar,3)
      get_ne_c_(-1);

   // !           if(myid==0)then
   // !               open(unit=11, file = 'testapars',status='unknown',action='write')
   // !               do j=0,jmx
                     
   // !                  write(11,*) apars(:,j,outk)-apars(:,j,outk+1)
   // !                  enddo
   // !                  close(11)

   // !               open(unit=11, file = 'testjpars',status='unknown',action='write')
   // !               do j=0,jmx
                     
   // !                  write(11,*) jpar(:,j,outk)-jpar(:,j,outk+1)
   // !                  enddo
   // !               close(11)   


   //  !              open(unit=11, file = 'testnes',status='unknown',action='write')
   //  !              do j=0,jmx
                     
   //  !                 write(11,*) denes(:,j,outk)-denes(:,j,outk+1)
   //  !                 enddo
   //  !                 close(11)

   //  !              open(unit=11, file = 'testBR',status='unknown',action='write')
   //  !              do j=0,jmx
                     
   //  !                 write(11,*) b0x(:,j)
   //  !                 enddo
   //  !                 close(11)
   //  !              end if

      if(ision==1) ppush_c_(timestep);
      if(ifluid==1){ 
         integ_c_(0); //again 1-index shananiganery
      } else {
         if(ision==1) ppush_c_(timestep);
            //if(ifluid==1)call pintef
         if(ifluid==1) integ_c_(0);
   // !             if(myid==0)then
   // !                open(unit=11, file = 'testden',status='unknown',action='write')
   // !                do j=0,jmx                 
   // !                  write(11,*) den2d2(:,j)
   // !                enddo
   // !                  close(11)
   // !              end if
      }
   // ! write(*,*)'dx=', dx, 'dz=',dz


   // !	   call accumulate(timestep,1)
   // !	   call ezamp
   // !	   call gkps
   // !	   call field(timestep,1)
      if(ifield_solver == 1) {
         phi.Clear();
         if(i3D != 0){
            for(k=myid*(kmx+1)/(numprocs); k <= (myid+1)*(kmx+1)/(numprocs)-1; ++k){
               for(iter = 0; iter<=iterations; ++iter){
                  fluxavg_c_(phi, phiavg);
                  if(eBoltzmann == 1){
                     BoltzSolve_c_(phi);
                  } else {
                     kval = (void*)k;
                     PetscCall(KSPSetComputeRHS(ksp,ComputeRHS,kval));
                     PetscCall(KSPSolve(ksp,NULL,NULL));
                     PetscCall(KSPGetSolution(ksp,&petsc_phi));
                     PetscCall(VecGetOwnershipRange(petsc_phi,&vec_start,&vec_end));
                     PetscCall(VecGetArrayRead(petsc_phi, &phi_array));
                     for(idx=1; idx <= (vec_end-vec_start); ++idx) {
                        i=((idx-1)%(iw))+is;
                        j=(idx-1)/(iw)+js;
                        phi(i,j,k)=phi_array[idx];//*mask(i,j);
                     }
                     PetscCall(VecRestoreArrayRead(petsc_phi,&phi_array));
                  }
               }
            }
            ierr = MPI_Allreduce(MPI_IN_PLACE, phi.start(), (imx+1)*(jmx+1)*(kmx+1),MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         } else {

            k = 0;

            for(iter = 0; iter <= iterations; ++iter) {
               fluxavg_c_(phi, phiavg);
               if(eBoltzmann == 1){
                  BoltzSolve_c_(phi);
               } else {
                  kval = (void*)k;
                  PetscCall(KSPSetComputeRHS(ksp,ComputeRHS,nullptr));
                  PetscCall(KSPSolve(ksp,nullptr,nullptr));
                  PetscCall(KSPGetSolution(ksp,&petsc_phi));
                  PetscCall(VecGetArrayRead(petsc_phi, &phi_array));

                  for(idx=1; idx <= (vec_end-vec_start); ++idx) {
                     i=(idx-1%(iw))+is;
                     j=(idx-1)/(iw)+js;
                     phi(i,j,k)=phi_array[idx];//*mask(i,j);
                  }



                  PetscCall(VecRestoreArrayRead(petsc_phi,&phi_array));

                  ierr = MPI_Allreduce(MPI_IN_PLACE, phi.start(), (imx+1)*(jmx+1)*(kmx+1),MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); CHKERRQ(ierr);

                  for(int i = 0; i <= imx; ++i){
                     for(int j = 0; j <= jmx; ++j){
                        for(int k = 0; k <= kmx; ++k){
                           phi(i,j,k) = phi(i,j,0);
                        }
                     }
                  }
               }
            }
         }
      }

      efieldcalc_c_(phi);
      if(i3D == 1){
         growthdiag_c_(phi);
      }
      if(myid == 0 && (timestep%10)==0){
         //write more stuff
      }
      cout << "stuff before" << endl;
      if(ision==1) cpush_c_(timestep);
      cout << "stuff after" << endl;
      //cintef(timestep);
      if(ifluid==1){
         integ_c_(2);
      } else {
         if(ision==1) cpush_c_(timestep);
         //if(ifluid==1) cintef(timestep)
         if(ifluid==1) integ_c_(2);
         //MPI_BARRIER(MPI_COMM_WORLD)
      }

      if(myid == 0 && (timestep%10)==0){
         //write more stuff
      }

      outd_c_(timestep);

      if(myid==0 && ifield_solver==1 && (timestep%10) == 0){
         //write a shit load of stuff
      }


      if(myid==master && (timestep%xnplt)==0){
         //plot some stuff
      }

   }
   end_total_tm = MPI_Wtime();
   total_tm = total_tm + end_total_tm - start_total_tm;
   
   ierr = MPI_Reduce(&ppush_tm, &tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); CHKERRQ(ierr);
   if(myid==0)ppush_tm = tmp/std::real(numprocs);
   ierr = MPI_Reduce(&cpush_tm, &tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); CHKERRQ(ierr);
   if(myid==0)cpush_tm = tmp/std::real(numprocs);
   ierr = MPI_Reduce(&integ_tm, &tmp, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD); CHKERRQ(ierr);
   if(myid==0)integ_tm = tmp/std::real(numprocs);
   MPI_Reduce(&total_tm, &tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); CHKERRQ(ierr);
   if(myid==0)total_tm = tmp/std::real(numprocs);
   //More writing stuff to files, work on after

   lasttm=MPI_Wtime();
   tottm=lasttm-starttm;

   if(eBoltzmann==0){
      PetscCall(PetscFinalize());
   }

   ierr = MPI_Finalize();
   cleanUpEquil();
   cleanupCom();
   return 0;
}

void initialize_c_(){
   double  dum, jacp; 
   //double dum1, dum2,  xndum, r, wx0, wx1
   // double x[2];
   // double y[2]; //0-1
   
   ppinit_c(myid, numprocs, ntube, kmx, i3D, GRID_COMM, TUBE_COMM, PETSC_COMM, petsc_color, petsc_rank);
  
   //reset timestep counter
   last = numprocs-1;
   timestep = 0;
   tcurr = 0.;

   for(int i = 0; i <= last; ++i){
      if(myid == i) init();
      ierr = MPI_Barrier(MPI_COMM_WORLD);
   }
   dum = 0.;
   for(int i = 0; i < imx; ++i){
      dum = dum+(jac[i]+jac[i+1])/2;
   }
   ierr = MPI_Allreduce(&dum, &jacp, 1, MPI_DOUBLE, MPI_SUM, TUBE_COMM);
   totvol = lx*lz*pi2*xctr;
   n0 = static_cast<float>(tmm[0]/totvol);


//         do k=0,kmx
//            den_pre(:,:,k)=xn0i
//         end do


   ierr = MPI_Barrier(MPI_COMM_WORLD);

   ncurr = 1;
}  

void init(){
   int ns, i, k, idum;
   double x, z, dum; //zdum
   double wx0, wx1, wz0, wz1, b;
   double bfldp, btorp, bxp, bzp, gt0ip, gt0ep, gn0ip, gn0ep, capnxp, capnzp;
   double upae0p; //??? seems weird to me, not calculated anywhere, maybe old/ got deleted?
   //std::complex<double> IU[2] = {0., 1.};
   //read values from gemx.in
   FILE *in_file = fopen("gemx.in", "r");
   if(in_file == NULL)
   {
      printf("Error! Could not open file: gemx.in\n");
      exit(-1); 
   } else 
   {
      fscanf(in_file, "%*[^\n]\n");
      fscanf(in_file, "%d %d %d %d %d %d %d", &imx, &jmx, &kmx, &mmx, &nmx, &nsmx, &ntube);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%lf %d %d %d", &dt, &nm, &nsm, &iez);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%d %d %d %d", &iput, &iget, &ision, &peritr);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%d %d", &nplot, &xnplt);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%lf %lf %lf", &cut, &amp, &tor);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%lf", &etaohm);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%d %lf %lf", &ifluid, &amie, &rneu);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%lf %d %d %lf", &betaVal, &nonlin, &nonline, &vcut);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%d %d %d %d %d %d %d %d", &ntracer, &ifield_solver, &i3D, &iBoltzmann, &eAdiabatic, &iterations, &icollision);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%d", &eBoltzmann);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%lf %d %lf %lf %lf %lf %lf", &psi_max, &psi_min, &R_min, &R_min, &Z_internal, &psi_div, &psi_a);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%d", &dbg);
   }
   fclose(in_file);

   nsm = 1;

   new_gemx_com(); //initializes arrays in gemx_com_c
   
   ns = 0;
   tmm[ns] = mmx;
   mm[ns] = mmx;
   mims[ns] = 2.0*1.67e-27;
   q[ns] = 1.0*1.6e-19; 
   lr[ns] = 4;
   
   emass = 1./amie;
   qel = -1;

   new_equil_c(); 
   
   lx = xdim;
   lz = zdim;

   if(myid == master){
      //plot - do later
   }

   iadi = 0;

   if(iget == 1) amp = 0.;

   dx =lx/std::real(imx);
   dz=lz/std::real(jmx);
   dzeta=pi2/(kmx+1);

   for(int i = 0; i <= imx; ++i){
      xg[i] = i*dx;
   }
   for(int k = 0; k <= kmx; ++k){
      zg[k] = k*dz; 
   }

   for(int i1 = 0; i1 <= imx; ++i1){
      x = i1*dx+xctr-xdim/2;
      i = int(x/dxeq);
      i = std::min(i,nx-1);
      wx0 = ((i+1)*dxeq-x)/dxeq;
      wx1 = 1.-wx0;

      for(int k1 = 0; k1 <= jmx; ++k1){
         z = k1*dz;
         k = int(z/dzeq);
         k = std::min(k,nz-1);            
         wz0 = ((k+1)*dzeq-z)/dzeq;
         wz1 = 1-wz0;

         bfldp = wx0*wz0*b0(i,k)+wx0*wz1*b0(i,k+1) 
                 +wx1*wz0*b0(i+1,k)+wx1*wz1*b0(i+1,k+1); 
         btorp = wx0*wz0*b0zeta(i,k)+wx0*wz1*b0zeta(i,k+1) 
               +wx1*wz0*b0zeta(i+1,k)+wx1*wz1*b0zeta(i+1,k+1); 
         bxp = wx0*wz0*b0x(i,k)+wx0*wz1*b0x(i,k+1) 
               +wx1*wz0*b0x(i+1,k)+wx1*wz1*b0x(i+1,k+1); 
         bzp = wx0*wz0*b0z(i,k)+wx0*wz1*b0z(i,k+1) 
               +wx1*wz0*b0z(i+1,k)+wx1*wz1*b0z(i+1,k+1); 
         gt0ip = wx0*wz0*t0i(i,k)+wx0*wz1*t0i(i,k+1) 
               +wx1*wz0*t0i(i+1,k)+wx1*wz1*t0i(i+1,k+1); 
         gt0ep = wx0*wz0*t0e(i,k)+wx0*wz1*t0e(i,k+1) 
               +wx1*wz0*t0e(i+1,k)+wx1*wz1*t0e(i+1,k+1); 
         gn0ip = wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) 
               +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1); 
         gn0ep = wx0*wz0*xn0e(i,k)+wx0*wz1*xn0e(i,k+1) 
               +wx1*wz0*xn0e(i+1,k)+wx1*wz1*xn0e(i+1,k+1); 
         capnxp = wx0*wz0*capnex(i,k)+wx0*wz1*capnex(i,k+1) 
               +wx1*wz0*capnex(i+1,k)+wx1*wz1*capnex(i+1,k+1); 
         capnzp = wx0*wz0*capnez(i,k)+wx0*wz1*capnez(i,k+1) 
               +wx1*wz0*capnez(i+1,k)+wx1*wz1*capnez(i+1,k+1); 

         b = 1.-tor+tor*bfldp;
         bmag(i1, k1) = b;
         gbtor(i1,k1) = btorp;
         gbx(i1,k1) = bxp;
         gbz(i1,k1) = bzp;              

         gt0i(i1,k1) = gt0ip;
         gt0e(i1,k1) = gt0ep;
         gn0e(i1,k1) = gn0ep;
         gn0i(i1,k1) = gn0ip;
         gcpnex(i1,k1) = capnxp;
         gcpnez(i1,k1) =  capnzp;          

         gupae0(i1,k1) = upae0p;
//         gnuoby(i1,k1) = (-dydrp*dnuobdtp+r0/q0*qhatp*dnuobdrp)*fp/radiusp*grcgtp
//         gnuobx(i1,k1) = dnuobdtp*fp/radiusp*grcgtp
      }
   }
   iseed = -(1777+myid*13);
   idum = ran2_c_(iseed);
   phi.Clear();
   apar.Clear();
   dene.Clear();
   upar.Clear();


   for(int i = 0; i <= imx; ++i){
      for(int j = 0; j <= jmx; ++j){
         for(int k = 0; k <= kmx; ++k){
            phi(i,j,k) = amp*(ran2_c_(idum)-0.5)*ifluid*1e-8;
            dene(i,j,k) = amp*(ran2_c_(idum)-0.5)*ifluid*1e-8;
            apar(i,j,k) = amp*(ran2_c_(idum)-0.5)*ifluid*1e-10; 
         }
      }
   }

   if(myid == master){
      //more plots
   }
}

void parperp_c_(double& vpar,double& vperp2, const int& m, const int& cnt){ 
   double r1 = 0.0;
   double r2 = 0.0;
   double t = 0.0;

   const double c0 = 2.515517;
   const double c1 = 0.802853;
   const double c2 = 0.010328;

   const double d1 = 1.432788;
   const double d2 = 0.189269;
   const double d3 = 0.001308;
   double temp = 0.0;
   int iflag = 1;


   r1 = revers_c_(m+myid*cnt, 7); //sets values for r1 and r2 (random numbers)
   r2 = revers_c_(m+myid*cnt, 11); 

   //.....quiet start---see denavit pf '71(?) & abramowitz hand book
   //.....fibonacci start---see denavit comm. pla. phy. & con. fus. '81
   // warning: we have g1=1 in the x-direction. This surpresses all odd
   //          modes in the x-direction!!!

   if (r1 <= 0.5) goto jump;

   r1 = 1.0 - r1;
   iflag = -1;

   jump:
   if(r1 >= 1.0e-6){
      t = sqrt(log(1.0/(r1*r1)));
   }
   else{
      t = 5.0;
      std::cout << "parperp2 warning m= " << m << std::endl;
   }

   temp = t-(c0+c1*t+c2*(t*t))/(1.+d1*t+d2*(t*t)+d3*(t*t*t));
   vpar = temp*iflag;

   vperp2 = -2.0*log(r2); 
   return;
}

void get_jpar_(CArray3D<double> &matrix){
   int i, j, k;

   for(k = 0; k <= kmx; ++k){     
      for(i = 2; i <= imx-2; ++i){
         for(j = 2; j <= jmx-2; ++j){
            if(mask3(i,j)>=2.99){ 
               jpar(i,j,k)=(-(matrix(i+1,j,k)+matrix(i-1,j,k)-2*matrix(i,j,k))/(dx*dx)     
                                 -(matrix(i,j+1,k)+matrix(i,j-1,k)-2*matrix(i,j,k))/(dz*dz)  
                                 -(matrix(i+1,j,k)-matrix(i-1,j,k))*0.5/(dx*Rgrid[i]/xu))  
                                 -q[0]*mu0*upar(i,j,k);
            }
         }
      }
   }
}

double ran2_c_(int& idum){ 
   const int IM1 = 2147483563;
   const int IM2 = 2147483399;
   const int IMM1 = IM1-1;
   const int IA1 = 40014;
   const int IA2 = 40692;
   const int IQ1 = 53668;
   const int IQ2 = 52774;
   const int IR1 = 12211;
   const int IR2 = 3791;
   const int NTAB = 32; 
   
   const double AM = 1.0/IM1; 
   const double NDIV = 1+IMM1/NTAB;
   const double EPS = 1.2e-7;
   const double RNMX = 1.0-EPS;

   int j = 0;
   int k = 0;

   static long idum2 = 123456789;
   static long iy = 0;
   static long iv[NTAB];

   double temp;
   double retVal;

   if(idum <= 0){
      if((-1 * idum) < 1){ 
         idum = 1;
      }
      else{ 
         idum = (-1 * idum);
      }
      idum2 = idum;
      for(j = NTAB+7; j >= 0; j--){
         k = idum/IQ1;
         idum = IA1*(idum-k*IQ1)-k*IR1;

         if(idum < 0) idum = idum+IM1;
         if(j < NTAB) iv[j] = idum;
      }
      iy = iv[0];
   }

   k = idum/IQ1;
   idum = IA1*(idum-k*IQ1)-k*IR1;
   if(idum < 0) idum = idum+IM1;
   k = idum2/IQ2;
   idum2 = IA2*(idum2-k*IQ2)-k*IR2;
   if(idum2 < 0) idum2 = idum2+IM2;  
   j = iy/NDIV;
   iy = iv[j]-idum2; 
   iv[j] = idum;
   if(iy<1) iy = iy+IMM1;
   temp = AM*iy;

   if(temp>RNMX){
      retVal = RNMX;
   }else{
      retVal = temp;
   }
   return retVal;
}

void loadi_c_(){ 
   int i = 0; 
   int j = 0;
   int k = 0;
   int m = 0;

   double vpar = 0.0;
   double vperp2 = 0.0;
   double r = 0;
   double x = 0;
   double z = 0;
   // double b = 0;
   double ter = 0;
   double bfldp = 0;

   double avgv = 0;
   double myavgv = 0;
   // double avgw = 0;
   double myavgw = 0;
   
   double dumx, dumy, dumz, jacp; //jacp used in initialize, not sure if value is supposed to be updated here since jacp not passed to function - currently does nothing
   double wx0, wx1, wz0, wz1;

   const long double pi2 = M_PI*2;

   cnt = static_cast<int>(tmm[0]/numprocs);
   cnt = mmx; 
   while(m < mm[0]){
   //load a slab of ions...

      //dumx=xdim*(ran2(iseed)+0.01)*0.9
      //dumy=zdim*(ran2(iseed)+0.01)*0.9
      //revers(MyId*cnt+j,2) !ran2(iseed)
      dumx=2*dxeq+(xdim-4*dxeq)*ran2_c_(iseed); //!revers(MyId*cnt+j,2) !ran2(iseed)
      dumy=2*dzeq+(zdim-4*dzeq)*ran2_c_(iseed); //!revers(MyId*cnt+j,3) !ran2(iseed)
      dumz=pi2*ran2_c_(iseed); //revers(MyId*cnt+j,5) !ran2(iseed)

      //dumx=dxeq+(xdim-2*dxeq)*m/((mm(1)))

      r = xctr-xdim/2+dumx;   
      jacp = r/(xctr+xdim/2);
//    if(ran2(iseed)<jacp){
//       x2(m)=min(dumx,xdim-dxeq)
//       z2(m)=min(dumy,zdim-dzeq)
//       x2(m)=max(dumx,dxeq)
//       z2(m)=max(dumz,dzeq)
//       }
      zeta2[m] = dumz;
      x2[m] = dumx;
      z2[m] = dumy;

      parperp_c_(vpar, vperp2, m+1, cnt);

      x = x2[m];
      i = static_cast<int>(x/dxeq);
      wx0 = ((i+1)*dxeq-x)/dxeq;
      wx1 = 1 - wx0;

      z = z2[m];
      k = static_cast<int>(z/dzeq);
      wz0 = ((k+1)*dzeq-z)/dzeq;
      wz1 = 1-wz0;

      bfldp = wx0*wz0*b0(i,k)+wx0*wz1*b0(i,k+1)+wx1*wz0*b0(i+1,k)+wx1*wz1*b0(i+1,k+1);
      ter = wx0*wz0*t0i(i,k)+wx0*wz1*t0i(i,k+1)+wx1*wz0*t0i(i+1,k)+wx1*wz1*t0i(i+1,k+1);
      u2[m] = vpar/sqrt(mims[0]/ter);
      mu[m] = 0.5*vperp2/bfldp*ter;

      myavgv = myavgv+u2[m];
//    LINEAR: perturb w(m) to get linear growth...
//       w2(m)=2.*amp*ran2(iseed)
      w2[m] = (wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) 
               +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1))*r/xctr*((imx-3)*(jmx-3)*(kmx+1))/(numprocs*mmx); //*xctr/(x+xctr-xdim/2.)
//    w2(m) = r/xctr*((imx-1)*(jmx-1)*(kmx+1))/(numprocs*mmx)

      myavgw += w2[m];
      m++;
   }
   //do i=1,mmx
   //avex=avex+x2(i)
   //end do
   //write(*,*)avex/mmx
   if(myid==0){ 
      std::ofstream myFile;
      std::string fileName = "testdepo_posi"; 
      myFile.open(fileName, std::ios::app);

      j = mmx-10001;
      while(j < mmx){
         myFile << std::setprecision(16) << x2[j] << "      ";
         myFile << std::setprecision(16) << z2[j] << "      ";
         myFile << std::setprecision(16) << zeta2[j] << "      " << std::endl;
         j++;
      }
      myFile.close();
   }

   myavgw = myavgw/mm[0];

   ierr = MPI_Allreduce(&myavgv, &avgv, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD);

   if(idg == 1) std::cout << "all reduce" << std::endl;
   avgv = avgv/static_cast<float>(tmm[0]);

   m = 0;
   do{
      u2[m] = u2[m]-avgv;
      x3[m] = x2[m];
      z3[m] = z2[m];
      zeta3[m] = zeta2[m];
      u3[m] = u2[m];
//    w2(m) = w2(m)-myavgw
      w3[m] = w2[m];
      m++;
   }while(m < mm[0]);

   return;
}

void gradu_c_(CArray3D<double> &u_, CArray3D<double> &ux_, CArray3D<double> &uz_){

   int ju = 0;
   int jl = 0;
   double ul = 0;
   
   for(int j = 0; j < jmx; ++j){
      ju = j+1;
      jl = j-1;
      if(j == 0) jl = jmx-1;
      for (int i = 0; i <= imx-1; ++i){ 
         for (int k = 0; k <= kmx; ++k){
            uz_(i,j,k) = (u_(i,ju,k)-u_(i,jl,k))/(2*dz);
         }
      }  
   }

   for(int i = 1; i <= imx-1; ++i){
      for(int j = 0; j <= jmx-1; ++j){
         for(int k = 0; k <= kmx; ++k){
            ux_(i,j,k) = (u_(i+1,ju,k)-u_(i-1,j,k))/(2*dx);
         }
      }
   }

   for(int j = 0; j <= jmx-1; ++j){
      for(int k = 0; k <= kmx; ++k){
         ul = u_(imx-1, j, k);
         ux_(0,j,k) = (u_(1,j,k)-ul)/(2*dx);
      }
   }
   return;
}

void gradz_c_(CArray3D<double> &u, CArray3D<double> &uz){
   int kleft, kright;
   double wx0, wx1, wz0, wz1, uleft, uright;

   for(int k = 0; k <= kmx; ++k)
   {
      kleft = k-1;
      if(k==0) kleft = kmx;
      kright = k+1;
      if(k==kmx) kright = 0;
      for(int i = 1; i <= imx; ++i){
         for(int j = 1; j <= jmx; ++j){
            wx0 = ((ileft(i,j)+1)*dx-xbackw(i,j))/dx;
            wx1 = 1-wx0;
            wz0 = ((jleft(i,j)+1)*dz-zbackw(i,j))/dz;
            wz1 = 1-wz0;
            uleft = wx0*wz0*u(ileft(i,j),jleft(i,j),kleft) 
                  +wx1*wz0*u(ileft(i,j)+1,jleft(i,j),kleft) 
                  +wx0*wz1*u(ileft(i,j),jleft(i,j)+1,kleft) 
                  +wx1*wz1*u(ileft(i,j)+1,jleft(i,j)+1,kleft);
            wx0 = ((iright(i,j)+1)*dx-xforw(i,j))/dx;
            wx1 = 1-wx0;
            wz0 = ((jright(i,j)+1)*dz-zforw(i,j))/dz;
            wz1 = 1-wz0;
            uright = wx0*wz0*u(iright(i,j),jright(i,j),kright) 
                  +wx1*wz0*u(iright(i,j)+1,jright(i,j),kright) 
                  +wx0*wz1*u(iright(i,j),jright(i,j)+1,kright) 
                  +wx1*wz1*u(iright(i,j)+1,jright(i,j)+1,kright);
            uz(i,j,k)=(uright-uleft)/(2*b0(i,j)/b0zeta(i,j)*dzeta*(Rgrid[i])/xu);
         }
      }
   }
}

void smooth_c_(CArray3D<double> &matrix_c, int &mk){
   CArray3D<double> temp_c;
   temp_c.resize(imx+1, jmx+1, kmx+1);
   
   if(mk == 3){
      for(int k = 0; k <= kmx; ++k){
         for(int i = 2; i <= imx-2; ++i){
            for(int j = 2; j <= jmx-2; ++j){
               if(!(mask3(i,j) < 2.99)){
                  temp_c(i,j,k) = (matrix_c(i,j,k)+matrix_c(i+1,j,k)+matrix_c(i,j+1,k)+matrix_c(i-1,j,k)+matrix_c(i,j-1,k))*0.2;
               } 
            }
         }
      }
   } else if (mk == 2) {
      for(int k = 0; k <= kmx; ++k){
         for(int i = 2; i <= imx-2; ++i){
            for(int j = 2; j <= jmx; ++j){
               if(!(mask2(i,j)<1.99)){
                  temp_c(i,j,k) =(matrix_c(i,j,k)+matrix_c(i+1,j,k)+matrix_c(i,j+1,k)+matrix_c(i-1,j,k)+matrix_c(i,j-1,k))*0.2;
               }  
            }
         }
      }
   } 
   matrix_c = temp_c;
}

void integ_c_(int iflag) { //fix in future
   int i, j, k;
   double wx0,wx1,wzeta0,wzeta1,wy0,wy1,x,z,zeta,R_major_over_R,R_major_over_R1;
   int start_integ_tm = MPI_Wtime();
   for(int i = 0; i <= imx; ++i) {
      for (int j = 0; j <= jmx; ++j) {
         for(int k = 0; k <= kmx; ++k) {
            den(iflag,i,j,k) = 0;
         }
      }
   }
   upar.Clear(); 

   #pragma acc parallel loop gang vector
   for(int m = 0; m < mm[0]; ++m) {

      x = x3[m];
      i = static_cast<int>(x/dxeq);
      wx0 = (i+1)-x/dxeq;
      wx1 = 1-wx0;

      R_major_over_R=xctr/(xctr-xdim/2+i*dx);
      R_major_over_R1=xctr/(xctr-xdim/2+(i+1)*dx);

      z = z3[m];
      j = static_cast<int>(z/dzeq);
      wy0 = (j+1)-z/dzeq;
      wy1 = 1-wy0;
      
      zeta= fmod(zeta3[m], pi2);
      k=static_cast<int>(zeta/dzeta);
      wzeta0=(k+1)-zeta/dzeta;
      wzeta1=1-wzeta0;

      #pragma acc atomic update 
      den(iflag,i,j,k)=den(iflag,i,j,k)+w3[m]*wx0*wy0*wzeta0*R_major_over_R;
      #pragma acc atomic update
      den(iflag,i+1,j,k)=den(iflag,i+1,j,k)+w3[m]*wx1*wy0*wzeta0*R_major_over_R1;
      #pragma acc atomic update
      den(iflag,i,j+1,k)=den(iflag,i,j+1,k)+w3[m]*wx0*wy1*wzeta0*R_major_over_R;
      #pragma acc atomic update
      den(iflag,i+1,j+1,k)=den(iflag,i+1,j+1,k)+w3[m]*wx1*wy1*wzeta0*R_major_over_R1;
      #pragma acc atomic update 
      upar(i,j,k)=upar(i,j,k)+u3[m]*w3[m]*wx0*wy0*wzeta0*R_major_over_R;
      #pragma acc atomic update
      upar(i+1,j,k)=upar(i+1,j,k)+u3[m]*w3[m]*wx1*wy0*wzeta0*R_major_over_R1;
      #pragma acc atomic update
      upar(i,j+1,k)=upar(i,j+1,k)+u3[m]*w3[m]*wx0*wy1*wzeta0*R_major_over_R;
      #pragma acc atomic update
      upar(i+1,j+1,k)=upar(i+1,j+1,k)+u3[m]*w3[m]*wx1*wy1*wzeta0*R_major_over_R1;

      if(k != kmx) {
         #pragma acc atomic update
         den(iflag,i,j,k+1)=den(iflag,i,j,k+1)+w3[m]*wx0*wy0*wzeta1*R_major_over_R;
         #pragma acc atomic update
         den(iflag,i+1,j,k+1)=den(iflag,i+1,j,k+1)+w3[m]*wx1*wy0*wzeta1*R_major_over_R1;
         #pragma acc atomic update
         den(iflag,i,j+1,k+1)=den(iflag,i,j+1,k+1)+w3[m]*wx0*wy1*wzeta1*R_major_over_R;
         #pragma acc atomic update
         den(iflag,i+1,j+1,k+1)=den(iflag,i+1,j+1,k+1)+w3[m]*wx1*wy1*wzeta1*R_major_over_R1;
         #pragma acc atomic update
         upar(i,j,k+1)=upar(i,j,k+1)+u3[m]*w3[m]*wx0*wy0*wzeta1*R_major_over_R;
         #pragma acc atomic update
         upar(i+1,j,k+1)=upar(i+1,j,k+1)+u3[m]*w3[m]*wx1*wy0*wzeta1*R_major_over_R1;
         #pragma acc atomic update
         upar(i,j+1,k+1)=upar(i,j+1,k+1)+u3[m]*w3[m]*wx0*wy1*wzeta1*R_major_over_R;
         #pragma acc atomic update
         upar(i+1,j+1,k+1)=upar(i+1,j+1,k+1)+u3[m]*w3[m]*wx1*wy1*wzeta1*R_major_over_R1;
   } else {
         #pragma acc atomic update
         den(iflag,i,j,0)=den(iflag,i,j,0)+w3[m]*wx0*wy0*wzeta1*R_major_over_R;
         #pragma acc atomic update
         den(iflag,i+1,j,0)=den(iflag,i+1,j,0)+w3[m]*wx1*wy0*wzeta1*R_major_over_R1;
         #pragma acc atomic update
         den(iflag,i,j+1,0)=den(iflag,i,j+1,0)+w3[m]*wx0*wy1*wzeta1*R_major_over_R;
         #pragma acc atomic update
         den(iflag,i+1,j+1,0)=den(iflag,i+1,j+1,0)+w3[m]*wx1*wy1*wzeta1*R_major_over_R1;
         #pragma acc atomic update
         upar(i,j,0)=upar(i,j,0)+u3[m]*w3[m]*wx0*wy0*wzeta1*R_major_over_R;
         #pragma acc atomic update
         upar(i+1,j,0)=upar(i+1,j,0)+u3[m]*w3[m]*wx1*wy0*wzeta1*R_major_over_R1;
         #pragma acc atomic update
         upar(i,j+1,0)=upar(i,j+1,0)+u3[m]*w3[m]*wx0*wy1*wzeta1*R_major_over_R;
         #pragma acc atomic update
         upar(i+1,j+1,0)=upar(i+1,j+1,0)+u3[m]*w3[m]*wx1*wy1*wzeta1*R_major_over_R1;
      }
   }
   #pragma acc wait
   
  
   ierr = MPI_Allreduce(MPI_IN_PLACE, &den(iflag,0,0,0), (imx+1)*(jmx+1)*(kmx+1), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   ierr = MPI_Allreduce(MPI_IN_PLACE, &upar(0,0,0), (imx+1)*(jmx+1)*(kmx+1),MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        

   for(int i = 0 ; i <= imx; ++i) {
      for(int j = 0; j <= jmx; ++j) {
         for(int k = 0; k <= kmx; ++k) {
            den(1,i,j,k) = den(iflag, i, j, k);
         }
      }
   }

   for(int i = 0; i <= imx; ++i) {
      for(int j = 0; j <= jmx; ++j) {
         den2d2(i,j) = 0;
      }
   }
   
   for(int i = 0 ; i <= imx; ++i) {
      for(int j = 0; j <= jmx; ++j) {
         for(int k = 0; k <= kmx; ++k) {
            den2d2(i,j)=den2d2(i,j)+den(iflag,i,j,k);
            if (i3D==0 && k != 0) upar(i,j,0) += upar(i,j,k);
         }
      }
   }

   for(i = 0; i <= imx; ++i) {
      for(j = 0; j <= jmx; ++j) {
         den2d2(i,j) = den2d2(i,j)/(kmx+1);
         if(i3D == 0) {
            upar(i,j,0) = upar(i,j,0)/(kmx+1);
            for(k = 1; k <= kmx; ++k) {
               upar(i,j,k) = upar(i,j,0);
            }
         }
      }
   }

   if(iflag==2) {
      for(i = 0; i <= imx; ++i) {
         for(j = 0; j <= jmx; ++j) {
            dden2d(i,j)-=den2d1(i,j);
            den2d1(i,j)=den2d2(i,j);
            // for(k = 0 ; k <= kmx; ++k){
            //    //den_pre=den(2,i,j,k); 
            // }
         }
      }
   }

   int end_integ_tm = MPI_Wtime();
   integ_tm = integ_tm + end_integ_tm - start_integ_tm; 
}

void get_apar_(const int &flagnumber) {
    CArray3D<double> gradPar;
    gradPar.resize(imx+1, jmx+1, kmx+1);
    gradpar_c_(phi, gradPar); 

   if(flagnumber == -1) {
      for(int i = 0; i <= imx; ++i){
         for(int j = 0; j <= jmx; ++j){
            for(int k = 0; k <= kmx; ++k){
               apars(i,j,k) = apar(i,j,k)-0.5*dt*(gradPar(i,j,k));  
            }
         }
      }
   } else {
      for(int i = 0; i <= imx; ++i){
         for(int j = 0; j <= jmx; ++j){
            for(int k = 0; k <= kmx; ++k){
               apar(i,j,k) = apar(i,j,k)-0.5*dt*(gradPar(i,j,k));
            }
         }
      }
   }
}

void get_ne_c_(int flagnumber){
   CArray3D<double> gradPar;
   gradPar.resize(imx+1, jmx+1, kmx+1);
   gradpar_c_(jpar, gradPar);

   if(flagnumber == -1) {
      for(int i = 0; i <= imx; ++i) {
         for(int j = 0; j <= jmx; ++j) {
            for(int k = 0; k <= kmx; ++k) {
               denes(i,j,k) = dene(i,j,k)+0.5*dt*(gradPar(i,j,k));  
            }
         }
      }
   } else if(flagnumber == 1) {
      for(int i = 0; i <= imx; ++i) {
         for(int j = 0; j <= jmx; ++j) {
            for(int k = 0; k <= kmx; ++k) {
               dene(i,j,k) = dene(i,j,k)+dt*(gradPar(i,j,k));
            }
         }
      }
   } else if(flagnumber == 0) {
      for(int i = 2; i <= imx-2; ++i) {
         for(int j = 2; j <= jmx-2; ++j) {
            for(int k = 0; k <= kmx; ++k) {
               if(mask4(i,j) == 4){
                  dene(i,j,k)=jpar(i,j,k)*sqrt(c2_over_vA2(i,j));
               }
            }
         }
      }
   }
   if(i3D == 0){
      for(int k = 1; k <= kmx; ++k) {
         if (flagnumber ==-1 ){
            for(int i = 0; i <= imx; ++i) {
               for(int j = 0; j <= jmx; ++j) {
                  denes(i,j,0) = denes(i,j,0)+denes(i,j,k);
               }
            }
         } else {
            for(int i = 0; i <= imx; ++i) {
               for(int j = 0; j <= jmx; ++j) {
                  dene(i,j,0) = dene(i,j,0)+dene(i,j,k);
               }
            }
         }
      }
      if(flagnumber == -1){
         for(int i = 0; i <= imx; ++i) {
            for(int j = 0; j <= jmx; ++j) {
               for(int k = 0; k <= kmx; ++k) {
                  denes(i,j,k) = denes(i,j,k)/(kmx+1);
               }
            }
         }
      } else {
         for(int i = 0; i <= imx; ++i) {
            for(int j = 0; j <= jmx; ++j) {
               for(int k = 0; k <= kmx; ++k) {
                  dene(i,j,k) = dene(i,j,k)/(kmx+1);
               }
            }
         }
      }
   }
}

void gradparz_c_(double *matrix){ //need const values for args, using const int later will work, issue with extern atm  
   double *gradparz = new double[(imx+1) * (jmx+1) * (kmx+1)];
   CArray3D<double> gradparz_c;
   CArray3D<double> matrix_c;
   gradparz_c.CreateArray3D(gradparz, imx+1, jmx+1, kmx+1);
   matrix_c.CreateArray3D(matrix, imx+1, jmx+1, kmx+1);
   gradz_c_(matrix_c, gradparz_c);
   for(int k = 0; k <= kmx; ++k){
      for(int i = 0; i <= imx; ++i){
         for(int j = 0; j <= jmx; ++j){
            gradparz_c(i,j,k) = 0.5*gradparz_c(i,j,k)*mask2(i,j);
         }
      }
   }
}

void gradpar_c_(CArray3D<double> &matrix, CArray3D<double> &gradPar){ 
   for(int k = 1; k <= kmx-1; ++k){
      for(int i = 2; i <= imx-2; ++i){
         for(int j = 2; j <= jmx-2; ++j){
            if (!(mask2(i,j)<1.99)){
               gradPar(i,j,k)=(b0x(i,j)/b0(i,j)*(matrix(i+1,j,k)-matrix(i-1,j,k))*0.5/dx  
               +b0z(i,j)/b0(i,j)*(matrix(i,j+1,k)-matrix(i,j-1,k))*0.5/dz                
               +b0zeta(i,j)/b0(i,j)*(matrix(i,j,k+1)-matrix(i,j,k-1))/((Rgrid[i]/xu)*2*dzeta));
            }
         }
      }
   }
   for(int i = 2; i <= imx-2; ++i){
      for(int j = 2; j <= jmx-2; ++j){
            if(!(mask2(i,j)<1.99)){
            gradPar(i,j,0)=(b0x(i,j)/b0(i,j)*(matrix(i+1,j,0)-matrix(i-1,j,0))*0.5/dx  
               +b0z(i,j)/b0(i,j)*(matrix(i,j+1,0)-matrix(i,j-1,0))*0.5/dz                
               +b0zeta(i,j)/b0(i,j)*(matrix(i,j,1)-matrix(i,j,kmx))/((Rgrid[i]/xu)*2*dzeta));
         }
      }
   }
   for(int i = 2; i <= imx-2; ++i){
      for(int j = 2; j <= jmx-2; ++j){
         if(!(mask2(i,j)<1.99)){
            gradPar(i,j,kmx)=(b0x(i,j)/b0(i,j)*(matrix(i+1,j,kmx)-matrix(i-1,j,kmx))*0.5/dx
               +b0z(i,j)/b0(i,j)*(matrix(i,j+1,kmx)-matrix(i,j-1,kmx))*0.5/dz    
               +b0zeta(i,j)/b0(i,j)*(matrix(i,j,0)-matrix(i,j,kmx-1))/((Rgrid[i]/xu)*2*dzeta));
         }
      }
   }
}

void pintef_c_(){
   int i,j,k;
   double ddedt;
   CArray3D<double> uz;
   uz.resize(imx+1, jmx+1, kmx+1);

   for(i = 0; i <= imx; ++i){
      for(j = 0; j <= jmx; ++j){
         for(k = 0; k <= kmx; ++k){
            phis(i,j,k) = phi(i,j,k);
            denes(i,j,k) = dene(i,j,k);
            apars(i,j,k) = apar(i,j,k);
         }
      }
   }

   gradz_c_(upar, uz); 
   for(k = 0; k <= kmx; ++k){
      for(i = 1; i <= imx; ++i){
         for(j = 1; j <= jmx; ++j){
            ddedt = -uz(i,j,k)*gn0e(i,j)*gbtor(i,j)/((xctr-xdim/2+xg[i])*bmag(i,j))
            +gn0e(i,j)*(gcpnex(i,j)*ez(i,j,k)-gcpnez(i,j)*ez(i,j,k))/bmag(i,j);
            dene(i,j,k) = denes(i,j,k)+0.5*dt*ddedt;
         }
      }
   }

   gradz_c_(phi, uz);
   for(k = 0; k <= kmx; ++k){
      for(i = 1; i <= imx; ++i){
         for(j = 1; j <= jmx; ++j){
            ddedt = -uz(i,j,k)*gbtor(i,j)/((xctr-xdim/2+xg[i])*bmag(i,j));
            apar(i,j,k) = apars(i,j,k)+0.5*dt*(ddedt+ezeta(i,j,k));
         }
      }
   }
}

PetscErrorCode ComputeInitialGuess(KSP ksp, Vec init_guess, void* ctx_void) {
    PetscFunctionBegin;

    //PetscInt* ctx = (PetscInt*) ctx_void;  // Cast void* to expected type

    PetscScalar h = 0.0;
    PetscCall(VecSet(init_guess, h));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode ComputeMatrix(KSP ksp, Mat AA, Mat BB, void* dummy) {
   DM dm;
   //int ii,jj;
   dummy = static_cast<int*>(dummy);
   PetscInt i,j,mx,my,xm;
   PetscInt    ym,xs,ys,i1, i5; 
   PetscScalar  v[5],Hx,Hy;
   PetscScalar  Hx2,Hy2; //tmp_r,a_value
   MatStencil   row[1],col[5]; 

   PetscInt ncols = 0;

   i1 = 1;
   i5 = 5;
   //a_value = 0.5;
   ierr = KSPGetDM(ksp,&dm); CHKERRQ(ierr);
   ierr = DMDAGetInfo(dm,nullptr,&mx,&my,nullptr,nullptr,nullptr,nullptr,
                      nullptr,nullptr,nullptr,nullptr,nullptr,nullptr); CHKERRQ(ierr);

   Hx = dx;//! (Rgrid(imx)-Rgrid(0)) / real(imx)
   Hy = dz;//(Zgrid(jmx)-Zgrid(0)) / real(jmx)

   Hx2 = Hx*Hx;
   Hy2 = Hy*Hy;
   PetscCall(DMDAGetCorners(dm,&xs,&ys,nullptr,&xm,&ym,nullptr));

   for(j=ys; j < ys+ym; ++j){
      for(i=xs; i < xs+xm; ++i) {
         ncols = 0;
         row[0].i = i;
         row[0].j = j;
         if(mask(i,j)<0.99){
            v[0] = c2_over_vA2(i,j)*(-2.0/Hx2-2.0/Hy2);
            ierr = MatSetValuesStencil(BB,i1,row,i1,row,v,INSERT_VALUES); CHKERRQ(ierr);
         } else {
            if(j > 0) {
               if(j == jmx) {
                  v[ncols] = c2_over_vA2(i,j)/Hy2-1.0/(2.0*Hy2)*( c2_over_vA2(i,j)- c2_over_vA2(i,j-1));
               } else {
                  v[ncols] = c2_over_vA2(i,j)/Hy2-1.0/(4.0*Hy2)*( c2_over_vA2(i,j+1)- c2_over_vA2(i,j-1));
               }
               col[ncols].i = i;
               col[ncols].j = j - 1;
               ncols++;
            }
            
            if(i > 0) {
               if(i == imx) {
                  v[ncols] =  c2_over_vA2(i,j)/Hx2-1.0/(2.0*Hx2)*( c2_over_vA2(i,j)- c2_over_vA2(i-1,j));
               } else {
                  v[ncols] =  c2_over_vA2(i,j)/Hx2-1.0/(4.0*Hx2)*( c2_over_vA2(i+1,j)- c2_over_vA2(i-1,j));
               }
               col[ncols].i = i - 1;
               col[ncols].j = j;
               ncols++;
            }

            v[ncols] = -2.0* c2_over_vA2(i,j) / Hx2 - 2.0* c2_over_vA2(i,j) / Hy2;
            col[ncols].i = i;
            col[ncols].j = j;
            //  write(*,*)v(3), xn0e(i,j)*mu0*e/t0e(i,j)
//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Boltzmann e!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
            if(eBoltzmann != 0) {
               v[ncols]-= xn0e(i,j)*mu0*e*e/t0e(i,j);
            }
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
            ncols++;

            if(i < imx){
               if(i == 0){
                  v[ncols] = c2_over_vA2(i,j)/Hx2+1.0/(2.0*Hx2)*( c2_over_vA2(i+1,j)- c2_over_vA2(i,j));
               } else {
                  v[ncols] =  c2_over_vA2(i,j)/Hx2+1.0/(4.0*Hx2)*( c2_over_vA2(i+1,j)- c2_over_vA2(i-1,j));
               }
               col[ncols].i = i + 1;
               col[ncols].j = j;
               ncols++;
            }
            
            if(j < jmx) {
               if(j == 0){
                  v[ncols] =  c2_over_vA2(i,j)/Hy2+1.0/(2.0*Hy2)*( c2_over_vA2(i,j+1)- c2_over_vA2(i,j));
               } else {
                  v[ncols] =  c2_over_vA2(i,j)/Hy2+1.0/(4.0*Hy2)*( c2_over_vA2(i,j+1)- c2_over_vA2(i,j-1));
               }
               col[ncols].i = i;
               col[ncols].j = j + 1;
               ncols++;
            }
   
            ierr = MatSetValuesStencil(BB, i1, row, ncols, col, v, INSERT_VALUES); CHKERRQ(ierr);
         }
      }
   }

   ierr = MatAssemblyBegin(BB,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   ierr = MatAssemblyEnd(BB,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   if(AA != BB) {
      ierr = MatAssemblyBegin(AA,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(AA,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   }
   //   PetscCall(MatView(AA,PETSC_VIEWER_STDOUT_WORLD));
   //   PetscCall(MatView(BB,PETSC_VIEWER_STDOUT_WORLD));

   return 0;
}

PetscErrorCode ComputeRHS(KSP ksp, Vec bbb, void* ctx) {
   int ii,jj,iflag;
   //PetscScalar* b_array = new PetscScalar[]; 
   PetscScalar  h,Hx,Hy;
   PetscInt  mx,my,i,j,xs,xm,ys,ym,vec_start,vec_end;
   DM dm;
   PetscInt idx;
   PetscScalar tmp_value = 0.0;
   // PetscScalar a_value,tmp_r;

   PetscInt k = (PetscInt)ctx; 
   tmp_value = 0;
   
   PetscCall(KSPGetDM(ksp,&dm));
   PetscCall(DMDAGetInfo(dm,nullptr,&mx,&my,nullptr,nullptr,nullptr,
                        nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr));
   PetscCall(VecGetOwnershipRange(bbb,&vec_start,&vec_end));
   PetscCall(DMDAGetCorners(dm,&xs,&ys,nullptr,&xm,&ym,nullptr));

   idx = vec_start-1;
   for(j = ys; j < ys+ym; ++j) {
      for(i = xs; i < xs+xm; ++i) {
         idx+=1;
         if(mask(i,j) < 0.99) {
            tmp_value = 0;
         } else {
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!2D ni noly now!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if(i3D == 0) {
               if(iBoltzmann == 0){
                  tmp_value = (denes(i,j,k) - q[0]*mu0*(den2d2(i,j) - xn0i(i,j)));
               } else {
                  tmp_value = -q[0]*mu0*(den2d2(i,j)-xn0i(i,j));
               }
//                tmp_value = 1
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3D case!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            } else {
               if(iBoltzmann == 0) {
                  tmp_value = denes(i,j,k)-q[0]*mu0*(den(1,i,j,k)-xn0i(i,j));
               }

               if(eAdiabatic != 0){
                  tmp_value = -q[0]*mu0*(den(1,i,j,k)-xn0i(i,j)) - (xn0e(i,j)*mu0*e*e/t0e(i,j))*phiavg(i,j);
               } else {
                  tmp_value = -q[0]*mu0*(den(1,i,j,k)-xn0i(i,j));
               }
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            }
         }
         PetscCall(VecSetValues(bbb,1,&idx, &tmp_value, INSERT_VALUES));
      }
   }

   PetscCall(VecAssemblyBegin(bbb));
   PetscCall(VecAssemblyEnd(bbb));
   
   return(PETSC_SUCCESS);
}

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALDER Flux Average SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void fluxavg_c_(CArray3D<double> &phi, CArray2D<double> &phiavg_in){
   //Currently only good for 2D case
   
   //Local Variables
   double phiavg1d[102];
   double psi1d[102];
   double phiavg1d_private[102];
   int gi = 0, xix, yjy, miw, psi_zero, store, k; 
   double weightinput,weightinput3D, phiavggi, psival, wmx0, wmx1;

   // //set all arrays to zero here.
   std::fill(std::begin(phiavg1d), std::end(phiavg1d), 0);
   std::fill(std::begin(psi1d), std::end(psi1d), 0);
   std::fill(std::begin(phiavg1d_private), std::end(phiavg1d_private), 0);

   psi_zero = 1;

   //Computation
   for (int line = 1; line <= 101; ++line) {
      phiavggi = 0;
      store = 0;
      for(gi = 0; gi < num_lines; ++gi) {
         if(gindex[gi] == line-1) {
            weightinput = (weight00[gi]*phi(iarray[gi], jarray[gi], 0) + 
                           weight10[gi]*phi(iarray[gi]+1, jarray[gi], 0) + 
                           weight01[gi]*phi(iarray[gi], jarray[gi]+1, 0) + 
                           weight11[gi]*phi(iarray[gi]+1, jarray[gi]+1, 0));

            if(i3D == 0) {
               phiavggi += (weightinput*jacobian[gi])/deno[gi]; 
            } else {
               for(k = 1; k <= kmx; ++k) {
                  weightinput3D = (weight00[gi]*phi(iarray[gi], jarray[gi],k) + 
                           weight10[gi]*phi(iarray[gi]+1, jarray[gi],k) + 
                           weight01[gi]*phi(iarray[gi], jarray[gi]+1,k) + 
                           weight11[gi]*phi(iarray[gi]+1, jarray[gi]+1,k));
               }
               weightinput += weightinput3D;
               phiavggi += (weightinput*jacobian[gi])/(deno[gi]*(kmx+1)); 
            }

            if(priv[gi] == 0) {
               store = gi;
            }
         }
         //Remove redundancy from closed loop integration process
         if(phiavggi != 0) {
            if(gindex[gi] != line-1) {
               if(i3D == 0){
                  phiavggi -= (weightinput*jacobian[gi-1])/deno[gi-1]; 
               }else{
                  phiavggi -= (weightinput*jacobian[gi-1])/(deno[gi-1]*(kmx+1));               
               }
               break;
            }
         }
      }
      if(store != 0) {
            phiavg1d[line] = phiavggi;
            psi1d[psi_zero] = psitab[store];
            psi_zero += 1;
         } else {
            phiavg1d_private[line] = phiavggi; 
         }
   }
   phiavg1d[0] = phiavg1d[1];
   //  phiavg1d(0) = input(268,254,0) //phiavg1d[0] = phi(268,254,0);
   
   // Save psi1d and timesteps of phiavg1d to understand convergence
   //  if (timestep == 10) then
   //     open(unit=11, file = 'psi1d',status='unknown',action='write')
   //                 write(11,*) psi1d(:)
   //              close(11)
   //  endif

   //  open(unit=11, file = 'phiavg1d',status='unknown',position='append')                
   //  write(11,*) phiavg1d(:)
   //  close(11)

   // Initialize output to zero

   //INTERPOLATION
   for(xix = 0; xix <= nx; ++xix) {
      for(yjy = 0; yjy <= nz; ++yjy) {
         psival = psi_p(xix, yjy);
         if(mask(xix,yjy) < 0.99) { 
            phiavg(xix,yjy) = 0;
         } else {
            miw  = int(psival/(psi1d[2]-psi1d[1]));
            wmx0 = ((miw+1)*(psi1d[2]-psi1d[1])-psival)/(psi1d[2]-psi1d[1]);
            wmx1 = 1-wmx0;
            if (yjy < 75 && xix < 150 && psival > 0.29 && psival < 0.31) { //Private region under X-point
               phiavg(xix,yjy) = wmx0*phiavg1d_private[miw] + wmx1*phiavg1d_private[miw+1];
            } else {
               phiavg(xix,yjy) = wmx0*phiavg1d[miw] + wmx1*phiavg1d[miw+1]; 
            }
         }
      }
   }
}

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALDER E FIELD SUBROUTINE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void efieldcalc_c_(CArray3D<double> &input_phi){ 
   //Input is phi array - labeled phi by extern

   //local vars
   int i, j, k, kminus, kplus;

   for(k = 0; k <= kmx; ++k){
      for(i = 2; i <= imx-1; ++i){
         for(j = 2; j <= jmx-1; ++j){
            ex(i,j,k) = -(input_phi(i+1,j,k) - input_phi(i-1,j,k))/(2*(Rgrid[1]-Rgrid[0]));
            ez(i,j,k) = -(input_phi(i,j+1,k) - input_phi(i,j-1,k))/(2*(Zgrid[1]-Zgrid[0]));
            if(k == 0){
               kminus = kmx;
               ezeta(i,j,k) = -(input_phi(i,j,k+1) - input_phi(i,j,kminus))/(2*Rgrid[i]*(2*M_PI/(kmx+1)));
            }
            if(k == kmx){
               kplus = 0;
               ezeta(i,j,k) = -(input_phi(i,j,kplus) - input_phi(i,j,k-1))/(2*Rgrid[i]*(2*M_PI/(kmx+1)));
            }else{
               ezeta(i,j,k) = -(input_phi(i,j,k+1) - input_phi(i,j,k-1))/(2*Rgrid[i]*(2*M_PI/(kmx+1)));
            }
         }
      }
   }
}

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Growth Rate Diagnostic!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void growthdiag_c_(CArray3D<double> &input_phi){
   //Local variables
   int store, k, gi, peak;
   double phiavgsq, weightinput, phiavggi;

   phiavgsq = 0;
   /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
   peak     = 45; //Manually set contour number of peak temperature gradient
   /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
   phiavggi = 0.0;
   store    = 0;

   for(gi = 21759; gi < 22540; ++gi){ //offset by one for gindex ptr indexing starting at 0 -Dom
      if(gindex[gi] == peak){
         for(k = 0; k <= kmx; ++k){
            weightinput = (weight00[gi]*input_phi(iarray[gi], jarray[gi],k)*input_phi(iarray[gi], jarray[gi],k) + 
                           weight10[gi]*input_phi(iarray[gi]+1, jarray[gi],k)*input_phi(iarray[gi]+1, jarray[gi],k) + 
                           weight01[gi]*input_phi(iarray[gi], jarray[gi]+1,k)*input_phi(iarray[gi], jarray[gi]+1,k) + 
                           weight11[gi]*input_phi(iarray[gi]+1, jarray[gi]+1,k)*input_phi(iarray[gi]+1, jarray[gi]+1,k));
         }
         phiavggi = phiavggi + (weightinput*jacobian[gi])/(deno[gi]*(kmx+1));

         if(priv[gi] == 0){
            store = gi;
         }
      }

      if(phiavggi != 0){
         if(gindex[gi] != peak){
            phiavggi = phiavggi - (weightinput*jacobian[gi-1])/(deno[gi-1]*(kmx+1));
            break;
         }
      }
   }

   if(store != 0){
      phiavgsq = phiavggi;
   }

   std::ofstream myFile("testphiavgsq", std::ios::app);
   if(myFile.is_open()){
      myFile << "            " << timestep << "    "; //looks weird, just making it look identical to the old test files
      myFile << std::setprecision(16) << phiavgsq << std::endl;
   }else{
      std::cerr << "Error opening testphiavgsq" << std::endl;
   }
   myFile.close();
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Boltzmann-Poisson Electron Solver!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void BoltzSolve_c_(CArray3D<double> &input_phi){ 
   //3D input phi array
   //Local Variables
   int i, j, k;

   for(i = 0; i <= imx; ++i){
      for(j = 0; j <= jmx; ++j){
         for(k = 0; k <= kmx; ++k){
            if(input_phi(i,j,k) == 0){
               phi_k(i,j,k) = 0.005+(ran2_c_(iseed)*0.005);
            }else{
               phi_k(i,j,k) = 1.01*input_phi(i,j,k);
            }
         }
      }
   }

   for(i = 0; i <= imx; ++i){
      for(j = 0; j <= jmx; ++j){
         for(k = 0; k <= kmx; ++k){
            dphidr(i,j,k)   = c2_over_vA2(i,j)*(input_phi(i+1,j,k)-input_phi(i-1,j,k))/(2*dx);
            dphi_kdr(i,j,k) = c2_over_vA2(i,j)*(phi_k(i+1,j,k)-phi_k(i-1,j,k))/(2*dx);

            dphidz(i,j,k)   = c2_over_vA2(i,j)*(input_phi(i,j+1,k)-input_phi(i,j-1,k))/(2*dz);
            dphi_kdz(i,j,k) = c2_over_vA2(i,j)*(phi_k(i,j+1,k)-phi_k(i,j-1,k))/(2*dz);
         }
      }
   }

    for(i = 0; i <= imx; ++i){
      for(j = 0; j <= jmx; ++j){
         for(k = 0; k <= kmx; ++k){
            d2phidr2(i,j,k)      = (dphidr(i+1,j,k)-dphidr(i-1,j,k))/(2*dx);
            d2phi_kdr2(i,j,k)    = (dphi_kdr(i+1,j,k)-dphi_kdr(i-1,j,k))/(2*dx);

            d2phidz2(i,j,k)      = (dphidz(i,j+1,k)-dphidz(i,j-1,k))/(2*dz);
            d2phi_kdz2(i,j,k)    = (dphi_kdz(i,j+1,k)-dphi_kdz(i,j-1,k))/(2*dz);
            
            OPPphi(i,j,k)  = (d2phidr2(i,j,k)+d2phidz2(i,j,k));
            OPPphik(i,j,k) = (d2phi_kdr2(i,j,k)+d2phi_kdz2(i,j,k));
   
            l_hand(i,j,k) = (OPPphik(i,j,k)-OPPphi(i,j,k))/(phi_k(i,j,k)-input_phi(i,j,k)) - (xn0e(i,j)*mu0*e*e/t0e(i,j)*exp(input_phi(i,j,k)*e/t0e(i,j)));

            if(i3D == 0){
               r_hand(i,j,k) = (OPPphi(i,j,k)+q[0]*mu0*den2d2(i,j)-e*mu0*xn0e(i,j)*exp((e/t0e(i,j))*(input_phi(i,j,k)))); //q set in gem_com externs check in on this
            }else{
               r_hand(i,j,k) = (OPPphi(i,j,k)+q[0]*mu0*den(1,i,j,k)-e*mu0*xn0e(i,j)*exp((e/t0e(i,j))*(input_phi(i,j,k)))); //den is a 4D array
            }

            input_phi(i,j,k) = phi_k(i,j,k) - (r_hand(i,j,k)/l_hand(i,j,k));
            // input_phi(i,j,k) = OPPphik(i,j,k)

            if(mask(i,j) < 0.99){
                  input_phi.Clear(); //TODO - DOUBLE CHECK THIS IS WHAT WE NEED TO DO
            }
         }
      }
   }
   // if (timestep == 1) then
   //    open(unit=11, file = 'testboltzmann',status='unknown',action='write')
   //    do j=0, jmx
   //    write(11,*) input_phi(:,j,0)
   //   end do
   //    close(11)
   // end if
}