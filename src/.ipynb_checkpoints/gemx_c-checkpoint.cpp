#include "gemx_c.hpp"
#include "pputil_c.hpp"
#include "fcnt.hpp"
#include "gemx_com_c.hpp"
#include "equil_c.hpp"
#include "ionPush_c.hpp"
#include "outd_c.hpp"

#include <fstream>
#include <iomanip>
#include <cmath>  //c++ version of math.h in c
#include <chrono> //used for timing, leaving in case anyone wants to use
#include <unistd.h>
#include <iostream>

#include <vector>

using namespace std;

int main() {
   ofstream file;
   PetscInt kval;
   int status,mid_i,mid_j;
   int n,i,j,k,ip,m,outk,ix=135,jx=68;
   int iter; //Calder Edit

   double tmp;
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
   
   while(dbg == 1){
      sleep(1);         //if debug option set sleep for forever. To release type "dbg = 0" into debug consol once attatched. Happy Hunting!
   }

   outk=0; //(kmx+1)/2

   one = 1;
   three = 3;


   if(eBoltzmann == 0) {

      PETSC_COMM_WORLD = PETSC_COMM;
      PetscCall(PetscInitialize(nullptr, nullptr, nullptr, nullptr));  
      
      PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp));
      // 5-point Stencil
      PetscCall(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, 
                             DMDA_STENCIL_STAR,imx+1,jmx+1,PETSC_DECIDE,PETSC_DECIDE,
                             one,one, nullptr, nullptr, &dm));
      // 9-point Stencil
      // PetscCall(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, 
                           //   DMDA_STENCIL_BOX,imx+1,jmx+1,PETSC_DECIDE,PETSC_DECIDE,
                           //   one,one, nullptr, nullptr, &dm));
      PetscCall(DMSetFromOptions(dm));
      PetscCall(DMSetUp(dm));
      PetscCall(KSPSetDM(ksp,dm));
      PetscCall(KSPSetComputeInitialGuess(ksp,ComputeInitialGuess,0)); 
      PetscCall(KSPSetComputeOperators(ksp,ComputeMatrix,0));    	
      PetscCall(DMDAGetCorners(dm,&is,&js,nullptr,&iw,&jw,nullptr));
      PetscCall(KSPSetFromOptions(ksp));
      PetscCall(KSPSetUp(ksp)); 
   } 

   if(iget == 0) loadi_c_();
   prepareDeviceData();
   integ_c_(1);
   if(myid == 0) {
      
      file.open("testden", ios::app);
      for(int j = 0; j <= jmx; ++j) {
         for(int i = 0; i <= imx; ++i) {
            file << den2d2(i,j) << "   ";
         }
         file << "\n";
      }
      file.close();
   }
   if(i3D == 0) {
      xn0i = den2d2;
   }

   starttm=MPI_Wtime();
   upar.Clear();

   // mid_i=imx/2; 
   // mid_j=jmx/2;
   // mid_i=257;
   // mid_j=257;
   mid_i = (imx+1)/2;
   mid_j = (jmx+1)/2;

   tor_n = 1;

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!initialize perturbation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (checkpoint == 0) {
      
   } else {
      //read files from out/checkpoint - written later (do both at same time is probably easiest) 
   }

   // phiavg.Clear(); phiavg all 0s by default

   get_jpar_(apar); 
   get_ne_c_(0);

   if(i3D==0){
      apar.Clear();
      dene.Clear();
      integ_c_(1);
   }


   if(myid==0){
      file.open("testj0");
      for(int i = 0 ; i <= imx; ++i){
         for(int j = 0; j <= jmx; ++j){
            file << jpar(i,j,outk) << "   ";
         }
         file << "\n";
      }
      file.close();

      file.open("test_particleload");
      for(int i = 0 ; i <= imx; ++i){
         for(int j = 0 ; j <= jmx; ++j){
            for(int k = 0 ; k <= kmx; ++k){
               if(mask(i,j)<0.99)continue;
               if(xn0i(i,j)<0.0)continue;
               double rel = (den(1,i,j,k) - xn0i(i,j)) / xn0i(i,j);
               file << rel << "  ";
               file << "\n";
            }
         }
      }
      file.close();

      file.open("test_Rweight_particleload_num");
      for(int i = 0 ; i <= imx; ++i){
         for(int j = 0 ; j <= jmx; ++j){
            for(int k = 0 ; k <= kmx; ++k){
               if(mask(i,j)<0.99)continue;
               if(xn0i(i,j)<0.0)continue;
               // double rel = (den(1,i,j,k) - xn0i(i,j)) / xn0i(i,j);
               file << (den(1,i,j,k) - xn0i(i,j))*Rgrid[i] << "  ";
               file << "\n";
            }
         }
      }
      file.close();

      file.open("test_Rweight_particleload_den");
      for(int i = 0 ; i <= imx; ++i){
         for(int j = 0 ; j <= jmx; ++j){
            for(int k = 0 ; k <= kmx; ++k){
               if(mask(i,j)<0.99)continue;
               if(xn0i(i,j)<0.0)continue;
               // double rel = (den(1,i,j,k) - xn0i(i,j)) / xn0i(i,j);
               file << xn0i(i,j)*Rgrid[i] << "  ";
               file << "\n";
            }
         }
      }
      file.close();

      file.open("testne0");
      for(int i = 0; i <= imx; ++i){
         for(int j = 0; j <= jmx; ++j){
            file << dene(i,j,outk) << "    ";
         }
         file << "\n";
      }
      file.close();

      file.open("testapar0");
      for(int i = 0; i <= imx; ++i){
         for(int j = 0; j <= jmx; ++j){
            file << apar(i,j,outk) << "    ";
         }
         file << "\n";
      }
      file.close();

      file.open("testne0_zeta");
      for(int k = 0; k <= kmx; ++k){
            file << dene(mid_i,mid_j,k) << "   \n";
      }
      file.close();

      file.open("testjpar0_zeta");
      for(int k = 0; k <= kmx; ++k){
            file << jpar(mid_i,mid_j,k) << "   \n";
      }
      file.close();
   }
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end of init perturbation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            

   if(ifield_solver == 1) ncurr=1;
       start_total_tm = MPI_Wtime();
   for(timestep=ncurr; timestep<=nm; ++timestep) {
      // start_total_tm = MPI_Wtime();
      for(int randTabInd = 0; randTabInd <= 10006; ++randTabInd){
         if(ran2_c_(iseed)-0.5 > 0){
            rand_table[randTabInd]=1;
         } else {
            rand_table[randTabInd]=-1;
         }
      }
      tcurr = tcurr+dt;

   //	   accumulate(timestep-1,0)
   //	   ezamp()
   //	   gkps()
   //    field(timestep-1,0)

   // ofstream weightFile;
   // if(myid==0){
   //    weightFile.open("testweights", ios::app);
   //    double weight_diag = 0.0;
   //    #pragma acc update self(zeta2[0:mmx], w2[0:mmx])
   //    for(int m = 0; m < mm[0]; ++m) {
   //     //  weight_diag =+ w2[m]/mm[0];
   //     weightFile << gw[m] << "\n";
   //    }
   //   // weightFile << timestep << "   " << weight_diag << "   " << zeta2[mmx-1] <<"\n";
     
   //    weightFile.close();
   // }
   // if (weightscheme == 0){
   //     // density_filter(den);
   // }    

   if(ifield_solver == 1) {
      phi.Clear();
      phiavg.Clear();
      denes=dene;

      if(eBoltzmann == 1) {
         k = 0;
         for(iter = 0; iter <= iterations; ++iter) {
            //  do iter=0, iterations
            //     call fluxavg(phi,phiavg)
            //  end do
            BoltzSolve_c_(phi);
         }
      
      // call  MPI_Allreduce(MPI_IN_PLACE, phi, (imx+1)*(jmx+1)*(kmx+1),MPI_Real8, MPI_SUM, MPI_COMM_WORLD,ierr)

      } else if(i3D != 0) {
         // cout << "i3d != 0" << endl;
         // phi.Clear();
         phiavg.Clear();
         for(iter = 0; iter <= iterations; ++iter) {
            phi.Clear(); // 9/12/2025 test
            // phiavg.Clear();
            // fluxavg_c_(phi,phiavg); //Uses previous time step phi
            for(k = myid*(kmx+1)/(numprocs); k < (myid+1)*(kmx+1)/(numprocs); ++k) {
               
               // fluxavg_c_(phi,phiavg); //Uses previous time step phi
               //phi.Clear();
               kval = (PetscInt)k;

               PetscCall(KSPSetComputeRHS(ksp,ComputeRHS,&kval));
               PetscCall(KSPSolve(ksp,nullptr,nullptr));
               PetscCall(KSPGetSolution(ksp,&petsc_phi));
               PetscCall(VecGetArrayRead(petsc_phi, &phi_array));
               PetscCall(VecGetOwnershipRange(petsc_phi,&vec_start,&vec_end));

               for(idx = 0; idx < (vec_end-vec_start); ++idx) {
                  i=((idx)%(iw))+is;
                  j=(idx)/(iw)+js;
                  phi(i,j,k)=phi_array[idx];//*mask(i,j);
               }
               PetscCall(VecRestoreArrayRead(petsc_phi,&phi_array));
            } 
            fluxavg_c_(phi,phiavg); //Uses previous time step phi
            MPI_Allreduce(MPI_IN_PLACE, phi.start(), (imx+1)*(jmx+1)*(kmx+1),MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            // fluxavg(phi,phiavg) //Turned off for CBC
         }
      
      } else {
         k = 0;
         phi.Clear();
         
         for(iter = 0; iter <= iterations; ++iter) {
            fluxavg_c_(phi, phiavg);
         
            kval = (PetscInt)k;
            PetscCall(KSPSetComputeRHS(ksp,ComputeRHS,&kval));
            PetscCall(KSPSolve(ksp,NULL,NULL));
            PetscCall(KSPGetSolution(ksp,&petsc_phi));
            PetscCall(VecGetOwnershipRange(petsc_phi,&vec_start,&vec_end));
            PetscCall(VecGetArrayRead(petsc_phi,&phi_array));
            for(idx = 0; idx < (vec_end-vec_start); ++idx) {
               i=((idx)%(iw))+is;
               j=(idx)/(iw)+js;
               phi(i,j,k)=phi_array[idx];//*mask(i,j);
            }

            PetscCall(VecRestoreArrayRead(petsc_phi,&phi_array));

            MPI_Allreduce(MPI_IN_PLACE, phi.start(), (imx+1)*(jmx+1)*(kmx+1),MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         }
         // MPI_Allreduce(MPI_IN_PLACE, phi.start(), (imx+1)*(jmx+1)*(kmx+1),MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD,ierr)
         for(int i = 0; i <= imx; ++i){ //Put outside Calder Loops
            for(int j = 0; j <= jmx; ++j){
               for(int kk = 1; kk <= kmx; ++kk){
                  phi(i,j,kk) = phi(i,j,0);
               }
            }
         }
      }

      //PING FUNCTION
      // if(timestep <= 10 && nonlin != 1) {
      //    for(int i = 0; i <= imx; ++i) {
      //       for (int j = 0; j <= jmx; ++j) {
      //          for(int k=0; k <= kmx; ++k) {
      //             //  phi(i,j,k) = 100
      //             // phi(i,j,k) = 1e-8*cos(modes*((pi2*k)/(kmx+1)-1.3*atan2(Zgrid[j]-Zgrid[jmx/2],Rgrid[i]-Rgrid[imx/2])))* 
      //             // exp(-pow((sqrt(pow((Zgrid[j]-Zgrid[jmx/2]),2)+pow((Rgrid[i]-Rgrid[imx/2]),2))-0.25),2)/(2*0.05*0.05)); //*cos(atan2(Zgrid(j)-Zgrid(jmx/2),Rgrid(i)-Rgrid(imx/2))/2)
      //             // phi(i,j,k) = 1e-5*cos(modes*(-1.3*atan2(Zgrid(j)+Zgrid(0)-Zgrid(jmx/2),Rgrid(i)+Rgrid(0)-Rgrid(imx/2)))) * &
      //             // exp(-((sqrt((Rgrid(i)+Rgrid(0)-Rgrid(imx/2))**2+(Zgrid(j)+Zgrid(0)-Zgrid(jmx/2))**2)-0.25)**2)/(2*(0.15)**2))
      //             if (mask(i,j) < 0.99) {
      //                // phi(i,j,k) = 0;
      //             }
      //          }
      //       }
      //    }
      // }
   
      if(modes >= 0) {
         fourier_modes(phi,modes);
      }

      //binomial_filter(phi); //TODO

      if(radial_filter == 1){
         // flux_fourier_filter(phi);
         radial_binomial_filter(phi);
      }
      if(hyper_filter == 1){
         hyperdiffusion_filter(phi);
      }

      if(filtering_iterations) {
         for(int filter_int = 1; filter_int <= filtering_iterations; ++filter_int) {
            poloidal_filter_methods(phi);
            // radial_binomial_filter(phi); 
         }
      }

      if(low_filter == 1){
         low_mode_filter(phi);
      }

      // //EFIELD TESTING
      efieldcalc_c_(phi);

      
      // get_apar_(-1);
      //smooth(apars,2);
      // get_jpar_(apars);
      //smooth(jpar,3)
      // get_ne_c_(-1);


      if(ision==1) ppush_c_(timestep);
      if(ifluid==1) integ_c_(0);


   } else {
         if(ision==1) ppush_c_(timestep);
         //if(ifluid==1)call pintef
         if(ifluid==1) integ_c_(0);
   }
   
     // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PHI DIAGNOSTIC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      phi_diag = 0.0;
      phi_diag_freq = 0.0;

      if(myid==0){
         for(int i = 0; i <= imx; ++i) {
            for(int j = 0; j <= jmx; ++j) {
               for(int k = 0; k <= kmx; ++k) {
                  phi_diag = phi_diag + pow(abs(phi(i,j,k)),2)/(imx*jmx*kmx);
               // phi_diag_freq = phi_diag_freq + phi(i,j,k)/(imx*jmx*kmx);
               }
            }
         }
         ofstream file;
         file.open("testPhiDiag", ios::app);
         file << timestep << "		" << phi_diag << endl;
         file.close();
      
         file.open("testPhiFreq", ios::app);
         file << timestep << "		" << phi(192,128,0) << "    " << phi(180,128,0) << endl;
         file.close();
      
         file.open("testPhiFreq2", ios::app);
         file << timestep << "		" << phi(170,128,0) << "    " << phi(200,128,0) << endl;
         file.close();
      }

   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (weightscheme == 0){
       // density_filter(den);
   } 

   if(ifield_solver == 1) {
      phi.Clear();

      if (eBoltzmann == 1) {
         k = 0;
         for(iter = 0; iter <= iterations; ++iter) {
            BoltzSolve_c_(phi);
         }
         //MPI_Allreduce(MPI_IN_PLACE, phi,start(), (imx+1)*(jmx+1)*(kmx+1),MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      } else if(i3D != 0) {
         // phi.Clear();
         phiavg.Clear();
         for(iter = 0; iter <= iterations; ++iter) {
            phi.Clear();
            // phiavg.Clear();
            for(k=myid*(kmx+1)/(numprocs); k < (myid+1)*(kmx+1)/(numprocs); ++k){

               if((timestep%1)==0){ // DELETE LATER
                  file.open("testmyid", ios::app);
                  file << myid << "    ";
                  file << "\n";
                  file.close();
               }

               // fluxavg_c_(phi, phiavg);

               kval = (PetscInt)k;
               PetscCall(KSPSetComputeRHS(ksp,ComputeRHS,&kval));
               PetscCall(KSPSolve(ksp,NULL,NULL));
               PetscCall(KSPGetSolution(ksp,&petsc_phi));
               PetscCall(VecGetOwnershipRange(petsc_phi,&vec_start,&vec_end));
               PetscCall(VecGetArrayRead(petsc_phi, &phi_array));

               for(idx=0; idx < (vec_end-vec_start); ++idx) {
                  i=((idx)%(iw))+is;
                  j=(idx)/(iw)+js;
                  phi(i,j,k)=phi_array[idx];//*mask(i,j);
               }
               PetscCall(VecRestoreArrayRead(petsc_phi,&phi_array));
            }
            fluxavg_c_(phi,phiavg); //Uses previous time step phi
            MPI_Allreduce(MPI_IN_PLACE, phi.start(), (imx+1)*(jmx+1)*(kmx+1),MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); CHKERRQ(ierr);
         }
      } else {

         k = 0;
         phi.Clear();

         for(iter = 0; iter <= iterations; ++iter) {
            fluxavg_c_(phi, phiavg);

            kval = (PetscInt)k;
            PetscCall(KSPSetComputeRHS(ksp,ComputeRHS,&kval));
            PetscCall(KSPSolve(ksp,nullptr,nullptr));
            PetscCall(KSPGetSolution(ksp,&petsc_phi));
            PetscCall(VecGetArrayRead(petsc_phi, &phi_array));

            for(idx=0; idx < (vec_end-vec_start); ++idx) {
               i=(idx%(iw))+is;
               j=(idx)/(iw)+js;
               phi(i,j,k)=phi_array[idx];//*mask(i,j);
            }

            PetscCall(VecRestoreArrayRead(petsc_phi,&phi_array));

            MPI_Allreduce(MPI_IN_PLACE, phi.start(), (imx+1)*(jmx+1)*(kmx+1),MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); CHKERRQ(ierr);
         }
         for(int i = 0; i <= imx; ++i){
            for(int j = 0; j <= jmx; ++j){
               for(int kk = 1; kk <= kmx; ++kk){
                  phi(i,j,kk) = phi(i,j,0); 
               }
            }
         }
      }

      //PING FUNCTION
      // if(timestep <= 10 && nonlin == 0) {
      //    for(int i = 0; i <= imx; ++i) {
      //       for(int j = 0; j <= jmx; ++j) {
      //          for(int k = 0; k <= kmx; ++k) {
      //             // phi(i,j,k) = 1e-8*cos(modes*((pi2*k)/(kmx+1)-1.3*atan2(Zgrid[j]-Zgrid[jmx/2],Rgrid[i]-Rgrid[imx/2])))* 
      //             // exp(-pow((sqrt(pow((Zgrid[j]-Zgrid[jmx/2]),2)+pow((Rgrid[i]-Rgrid[imx/2]),2))-0.25),2)/(2*0.05*0.05));//*cos(atan2(Zgrid(j)-Zgrid(jmx/2),Rgrid(i)-Rgrid(imx/2))/2)
      //          // ! phi(i,j,k) = 1e-5*cos(modes*(pi2*k-1.3*atan2(Zgrid(j)+Zgrid(0)-Zgrid(jmx/2),Rgrid(i)+Rgrid(0)-Rgrid(imx/2)))) * &
      //          // ! exp(-((sqrt((Rgrid(i)+Rgrid(0)-Rgrid(imx/2))**2+(Zgrid(j)+Zgrid(0)-Zgrid(jmx/2))**2)-0.25)**2)/(2*(0.15)**2))
      //             if (mask(i,j) < 0.99) {
      //                // phi(i,j,k) = 0;
      //             }
      //          }
      //       }
      //    }
      // }

      if(modes >= 0) {
         fourier_modes(phi,modes);
      }

   //    !!!!!!!!!!!!!!!!!!!!!!!
   // ! call ftcamp(phi,timestep)
   // !!!!!!!!!!!!!!!!!!!!!!!
   // ! call binomial_filter(phi)


      if(radial_filter == 1){
         // flux_fourier_filter(phi);
         radial_binomial_filter(phi);
      }
      if(hyper_filter == 1){
         hyperdiffusion_filter(phi);
      }

      if(filtering_iterations) {
         for(int filter_int = 1; filter_int <= filtering_iterations; ++filter_int) {
            poloidal_filter_methods(phi);
            // radial_binomial_filter(phi);
         }
      }

      if(low_filter == 1){
         low_mode_filter(phi);
      }

      //EFIELD TESTING
      efieldcalc_c_(phi);
      

      if(i3D == 1){
         // growthdiag_c_(phi); //deprecated
      }


      if(myid == 0 && (timestep%1)==0){
         cout << "outk=" << outk << "\n";

         file.open("testphi");
         for(int j = 0; j <= jmx; ++j)  {
            for(int i = 0; i <= imx; ++i) {
               file << phi(i,j,outk) << "    ";
            }
            file << "\n";
         }
         file.close();
      }

      // get_apar_(1);
      //  !call smooth_c(apar,2)
      // !call get_jpar(apar)
      // get_jpar_(apar);
      //call smooth(jpar,3)
      // get_ne_c_(1);

      if(myid == 0 && (timestep%1) == 0) {
         file.open("testphiavg");
         for(int j = 0; j <= jmx; ++j)  {
            for(int i = 0; i <= imx; ++i) {
               file << phiavg(i,j) << "    ";
            }
            file << "\n";
         }
         file.close();

         file.open("testER");
         for(int j = 0; j <= jmx; ++j)  {
            for(int i = 0; i <= imx; ++i) {
               file << ex(i,j,0) << "    ";
            }
            file << "\n";
         }
         file.close();

         file.open("testEZ");
         for(int j = 0; j <= jmx; ++j)  {
            for(int i = 0; i <= imx; ++i) {
               file << ez(i,j,0) << "    ";
            }
            file << "\n";
         }
         file.close();
      }

      if(ision==1) cpush_c_(timestep); 
      //cintef(timestep);
      if(ifluid==1) integ_c_(1);
   } else {
      if(ision==1) cpush_c_(timestep);
      //if(ifluid==1) cintef(timestep)
      if(ifluid==1) integ_c_(1);
      //MPI_BARRIER(MPI_COMM_WORLD)
   }

      if(myid == 0 && (timestep%1)==0){
         file.open("testden2");
         for(int i = 0; i <= imx; ++i) {
            for(int j = 0; j <= jmx; ++j){
               file << den2d2(i,j) << "    ";
            }
            file << "\n";
         }
         file.close();

         file.open("testdiffden");
         for(int i = 0; i <= imx; ++i) {
            for(int j = 0; j <= jmx; ++j){
               file << dden2d(i,j) << "   ";
            }
            file << "\n";
         }
         file.close();

         file.open("testupar");
         for(int i = 0; i <= imx; ++i) {
            for(int j = 0; j <= jmx; ++j){
               file << upar(i,j,0) << "   ";
            }
            file << "\n";
         }
         file.close();
      }


      outd_c_(timestep);

      if(myid==0 && ifield_solver==1 && (timestep%1) == 0){
         file.open("testapar");
         for(int i = 0; i <= imx; ++i) {
            for(int j = 0; j <= jmx; ++j){
               file << apar(i,j,outk) << "   ";
            }
            file << "\n";
         }
         file.close();

         file.open("testjpar");
         for(int i = 0; i <= imx; ++i) {
            for(int j = 0; j <= jmx; ++j){
               file << jpar(i,j,outk) << "   ";
            }
            file << "\n";
         }
         file.close();

         file.open("testne");
         for(int i = 0; i <= imx; ++i) {
            for(int j = 0; j <= jmx; ++j){
               file << dene(i,j,outk) << "   ";
            }
            file << "\n";
         }
         file.close();

         file.open("testphi_r_phi");
         for(int i = 0; i <= imx; ++i) {
            for(int k = 0; k <= kmx; ++k){
               file << phi(i,mid_j,k) << "   ";
            }
            file << "\n";
         }
         file.close();
      }


      if(myid==master && (timestep%xnplt)==0 && (timestep%50)==0){
         file.open("testPhit", ios::app);
         for(int i = 0; i <= imx; ++i) {
            for(int j = 0; j <= jmx; ++j){
               file << phi(i,j,outk) << "   ";
            }
            file << "\n";
         }
         file.close();
      }
      
      if(myid == master && ifield_solver == 1) {
         cout << "time_step=" << timestep << "\n";
         cout << "dx=" << dx << "dz=" << dz << "   dzeta=" << dzeta << "   omega_A0=" << tor_n/(Rgrid[mid_i]/xu*sqrt(c2_over_vA2(mid_i,mid_j))) << "\n";
         cout << "v_A=" << 1/sqrt(c2_over_vA2(mid_i,mid_j)) << "  Omega_i=" << q[0]*b0(mid_i,mid_j)/mims[0] << "\n";
      }
      end_total_tm = MPI_Wtime();
      double totTime = end_total_tm - start_total_tm;
      // cout << "totTIme = " << totTime << "\n";
     
   }
   total_tm = total_tm + end_total_tm - start_total_tm;
   
   ierr = MPI_Reduce(&ppush_tm, &tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); CHKERRQ(ierr);
   if(myid==0)ppush_tm = tmp/std::real(numprocs);
   ierr = MPI_Reduce(&cpush_tm, &tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); CHKERRQ(ierr);
   if(myid==0)cpush_tm = tmp/std::real(numprocs);
   ierr = MPI_Reduce(&integ_tm, &tmp, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD); CHKERRQ(ierr);
   if(myid==0)integ_tm = tmp/std::real(numprocs);
   MPI_Reduce(&total_tm, &tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); CHKERRQ(ierr);
   if(myid==0)total_tm = tmp/std::real(numprocs);

   if (myid == 0) {
    file.open("gemx_timing.txt");
    file << ppush_tm << " "
         << cpush_tm << " "
         << integ_tm << " "
         << total_tm << "\n";
    file.close();
   }

   lasttm=MPI_Wtime();
   tottm=lasttm-starttm;

   if(eBoltzmann==0){
      PetscCall(PetscFinalize());
   }

   ierr = MPI_Finalize();
   cleanUpEquil();
   cleanupCom();
   freeDeviceData();
   return 0;
}

void initialize_c_(){
   double  dum, jacp; 
  
   //reset timestep counter
   last = numprocs-1;
   timestep = 0;
   tcurr = 0.;
   init();
   ppinit_c(myid,numprocs,ntube,kmx,i3D,TUBE_COMM,GRID_COMM, PETSC_COMM,petsc_color,petsc_rank);
   iseed = -(1777+myid*13);
   // for(int i = 0; i <= last; ++i){
   //    if(myid == i) {
         
   //    }
   //       ierr = MPI_Barrier(MPI_COMM_WORLD);
   // }

   dum = 0.;
   for(int i = 0; i < imx; ++i){
      dum = dum+(jac[i]+jac[i+1])/2;
   }
   MPI_Allreduce(&dum, &jacp, 1, MPI_DOUBLE, MPI_SUM, TUBE_COMM);
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
   } else {
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
      fscanf(in_file, "%lf %d %lf", &betaVal, &nonlin, &vcut);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%d %d %d %d %d %d %d", &ntracer, &ifield_solver, &i3D, &iBoltzmann, &eAdiabatic, &iterations, &icollision);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%d", &eBoltzmann);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%d %d", &iflr, &PADE);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%d", &CST);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%d %d", &weightscheme, &loadingscheme);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%d %d %d", &modes, &filtering_iterations, &low_filter);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%d %d %d", &cold_start, &hyper_filter, &radial_filter);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%lf %d %lf %lf %lf %lf %lf", &psi_max, &psi_min, &R_min, &Z_min, &Z_internal, &psi_div, &psi_a);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%d", &checkpoint);
      fscanf(in_file, " %*[^\n]\n");
      fscanf(in_file, "%d", &dbg);
   }
   fclose(in_file);
   nsm = 1;

   new_gemx_com(); //initializes arrays in gemx_com_c
   nx = imx;
   nz = jmx;
   nzeta = kmx;  
   
   ns = 0; //set to nsm
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
         bmag(i1,k1) = b;
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
      std::cout << "parperp2 warning m= " << m << "\n";
   }

   temp = t-(c0+c1*t+c2*(t*t))/(1.+d1*t+d2*(t*t)+d3*(t*t*t));
   vpar = temp*iflag;

   vperp2 = -2.0*log(r2);

   if (vpar > 3) {
      vpar = 3;
   }
   if (vpar < -3) {
      vpar = -3;
   }

   if (vperp2 > 9) {
      vperp2 = 9;
   }
}

void get_jpar_(CArray3D<double> &matrix){
   int i, j, k;

   for(i = 2; i <= imx-2; ++i) {     
      for(j = 2; j <= jmx-2; ++j) {
         for(k = 0; k <= kmx; ++k) {
            if(!(mask3(i,j) < 2.99)){
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
      if((-idum) < 1){ 
         idum = 1;
      }
      else{ 
         idum = (-idum);
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

   double realparticles = 0;
   
   double q_safety = 0;

   double avgv = 0;
   double myavgv = 0;
   // double avgw = 0;
   double myavgw = 0;
   
   double dumx, dumy, dumz, jacp; //jacp used in initialize, not sure if value is supposed to be updated here since jacp not passed to function - currently does nothing
   double wx0, wx1, wz0, wz1;

   cnt = static_cast<int>(tmm[0]/numprocs);
   cnt = mmx; 

   if(CST != 0) {
      read1D("rdata.dat", x2, 0);
      read1D("zdata.dat", z2, 0);
   }

   for(i = 2; i <= imx-2; ++i) {
      for(j = 2; j <= jmx-2; ++j){
         realparticles = realparticles + xn0i(i,j)*Rgrid[i];
      }
   }

   while(m < mm[0]) {
   //load a slab of ions...

      //dumx=xdim*(ran2(iseed)+0.01)*0.9
      //dumy=zdim*(ran2(iseed)+0.01)*0.9
      //revers(MyId*cnt+j,2) !ran2(iseed)
      dumx=2*dxeq+(xdim-4*dxeq)*ran2_c_(iseed); //!revers(MyId*cnt+j,2) !ran2(iseed)
      dumy=2*dzeq+(zdim-4*dzeq)*ran2_c_(iseed); //!revers(MyId*cnt+j,3) !ran2(iseed)
      dumz=pi2*ran2_c_(iseed); //revers(MyId*cnt+j,5) !ran2(iseed)

      //dumx=dxeq+(xdim-2*dxeq)*m/((mm(1)))

      if(CST != 0) {
         dumx = x2[m]-xctr+xdim/2;
         dumy = z2[m]-zctr+zdim/2;
      }

      r = xctr-xdim/2+dumx;   
      jacp = r/(xctr+xdim/2);
     
      i = static_cast<int>(dumx/dxeq);
      wx0 = ((i+1)*dxeq-dumx)/dxeq;
      wx1 = 1.-wx0;

      k = static_cast<int>(dumy/dzeq);
      wz0 = ((k+1)*dzeq-dumy)/dzeq;
      wz1 = 1.-wz0;

      if(loadingscheme==1){
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

         bfldp = wx0*wz0*b0(i,k)+wx0*wz1*b0(i,k+1) 
                  +wx1*wz0*b0(i+1,k)+wx1*wz1*b0(i+1,k+1); 
         ter = wx0*wz0*t0i(i,k)+wx0*wz1*t0i(i,k+1) 
                  +wx1*wz0*t0i(i+1,k)+wx1*wz1*t0i(i+1,k+1);

         u2[m] = vpar/sqrt(mims[0]/ter);
         mu[m] = 0.5*vperp2/bfldp*ter;

         myavgv = myavgv+u2[m];

         w2[m] = (wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) 
                  +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1))*r/xctr*((imx-3)*(jmx-3)*(kmx+1))/(numprocs*mmx);
         gw[m] = 1;
         if(weightscheme == 1) {
            if(nonlin == 1) {
               // w2[m] = 0;
               // w2[m] = 1e-12*(ran2_c_(iseed)-0.5);
               w2[m] = 1e-3;
            } else {
               q_safety = 2.52*pow(sqrt(pow(dumx + Rgrid[0] - Rgrid[imx / 2], 2) + pow(dumy + Zgrid[0] - Zgrid[jmx / 2], 2)), 2) -
                        0.16*sqrt(pow(dumx + Rgrid[0] - Rgrid[imx / 2], 2) + pow(dumy + Zgrid[0] - Zgrid[jmx / 2], 2)) + 0.86;

               w2[m] = 1e-12 * cos(modes*zeta2[m] - round(modes*q_safety)*atan2(dumy + Zgrid[0] - Zgrid[jmx / 2], dumx + Rgrid[0] - Rgrid[imx / 2])) *
                     exp(-pow(sqrt(pow(dumx + Rgrid[0] - Rgrid[imx / 2], 2) + pow(dumy + Zgrid[0] - Zgrid[jmx / 2], 2)) - 0.3, 2) / (2 * pow(0.05, 2)));

               // w2[m] = 1e-12 * cos(modes * (zeta2[m] - 1.4 * atan2(dumy + Zgrid[0] - Zgrid[jmx / 2], dumx + Rgrid[0] - Rgrid[imx / 2]))) * 
                     // exp(-pow(sqrt(pow(dumx + Rgrid[0] - Rgrid[imx / 2], 2) + pow(dumy + Zgrid[0] - Zgrid[jmx / 2], 2)) - 0.3, 2) / (2 * pow(0.05, 2)));
                     // * cos(atan2(Zgrid[j] - Zgrid[jmx / 2], Rgrid[i] - Rgrid[imx / 2]) / 2);
            }
            gw[m] = (wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) 
                  +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1))*r/xctr*((imx-3)*(jmx-3)*(kmx+1))/(numprocs*mmx);
         }

         myavgw += w2[m];
         m++;
      } else {
         if(ran2_c_(iseed)<jacp && ran2_c_(iseed) < (wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1)+wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1))/xn0i(imx/2,jmx/2)) {

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

            bfldp = wx0*wz0*b0(i,k)+wx0*wz1*b0(i,k+1) 
                     +wx1*wz0*b0(i+1,k)+wx1*wz1*b0(i+1,k+1); 
            ter = wx0*wz0*t0i(i,k)+wx0*wz1*t0i(i,k+1) 
                     +wx1*wz0*t0i(i+1,k)+wx1*wz1*t0i(i+1,k+1);

            u2[m] = vpar/sqrt(mims[0]/ter);
            mu[m] = 0.5*vperp2/bfldp*ter;

            myavgv = myavgv+u2[m];
      //    LINEAR: perturb w(m) to get linear growth...
      //       w2(m)=2.*amp*ran2(iseed)
            // w2[m] = (wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) 
                     // +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1))*r/xctr*((imx-3)*(jmx-3)*(kmx+1))/(numprocs*mmx); //*xctr/(x+xctr-xdim/2.)
            w2[m] = (kmx+1)*realparticles/(numprocs * mmx * xctr);
            // gw[m] = 1;
            gw[m] = 1*1/1.007;

            if(weightscheme == 1) {
               if(nonlin == 1) {
                  // w2[m] = 0;
                  // w2[m] = 1e-12*(ran2_c_(iseed)-0.5);
                  w2[m] = 1e-3 * cos(modes*zeta2[m] - round(modes*q_safety)*atan2(dumy + Zgrid[0] - Zgrid[jmx / 2], dumx + Rgrid[0] - Rgrid[imx / 2])) *
                        exp(-pow(sqrt(pow(dumx + Rgrid[0] - Rgrid[imx / 2], 2) + pow(dumy + Zgrid[0] - Zgrid[jmx / 2], 2)) - 0.3, 2) / (2 * pow(0.05, 2)));
                  if (modes == -1){
                     w2[m] = 1e-3;
                  }
                  // w2[m] = 1e-12 * cos(modes*zeta2[m] - round(modes*q_safety)*atan2(dumy + Zgrid[0] - Zgrid[jmx / 2], dumx + Rgrid[0] - Rgrid[imx / 2])) *
                  //       exp(-pow(sqrt(pow(dumx + Rgrid[0] - Rgrid[imx / 2], 2) + pow(dumy + Zgrid[0] - Zgrid[jmx / 2], 2)) - 0.3, 2) / (2 * pow(0.05, 2)));
               } else {
                  q_safety = 2.52*pow(sqrt(pow(dumx + Rgrid[0] - Rgrid[imx / 2], 2) + pow(dumy + Zgrid[0] - Zgrid[jmx / 2], 2)), 2) -
                        0.16*sqrt(pow(dumx + Rgrid[0] - Rgrid[imx / 2], 2) + pow(dumy + Zgrid[0] - Zgrid[jmx / 2], 2)) + 0.86;

                  w2[m] = 1e-12 * cos(modes*zeta2[m] - round(modes*q_safety)*atan2(dumy + Zgrid[0] - Zgrid[jmx / 2], dumx + Rgrid[0] - Rgrid[imx / 2])) *
                        exp(-pow(sqrt(pow(dumx + Rgrid[0] - Rgrid[imx / 2], 2) + pow(dumy + Zgrid[0] - Zgrid[jmx / 2], 2)) - 0.3, 2) / (2 * pow(0.05, 2)));

                  // w2[m] = 1e-12 * cos(modes * (zeta2[m] - 1.4 * atan2(dumy + Zgrid[0] - Zgrid[jmx / 2], dumx + Rgrid[0] - Rgrid[imx / 2]))) * 
                        // exp(-pow(sqrt(pow(dumx + Rgrid[0] - Rgrid[imx / 2], 2) + pow(dumy + Zgrid[0] - Zgrid[jmx / 2], 2)) - 0.3, 2) / (2 * pow(0.05, 2)));
                        // * cos(atan2(Zgrid[j] - Zgrid[jmx / 2], Rgrid[i] - Rgrid[imx / 2]) / 2);
               }
               // gw[m] = (wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) + wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1))*r/xctr*((imx-3)*(jmx-3)*(kmx+1))/(numprocs*mmx);
               gw[m] = (kmx+1)*realparticles/(numprocs * mmx * xctr);
            }

            myavgw += w2[m];
            m++;
         }
      }      
   }
   
   if(CST != 0) {
      //write to some files
   }

   if(myid==0){
      std::ofstream myFile;
      std::string fileName = "testdepo_posi"; 
      myFile.open(fileName, std::ios::app);

      j = mmx-10001;
      while(j < mmx){
         myFile << std::setprecision(16) << x2[j] << "      ";
         myFile << std::setprecision(16) << z2[j] << "      ";
         myFile << std::setprecision(16) << zeta2[j] << "      " << "\n";
         j++;
      }
      myFile.close();
   }
   myavgw = myavgw/mm[0];

   MPI_Allreduce(&myavgv, &avgv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   if(idg == 1) std::cout << "all reduce" << "\n";
   avgv = avgv/static_cast<float>(tmm[0]);

   for(int m = 0; m < mm[0]; ++m){
      u2[m] = u2[m]-avgv;
      x3[m] = x2[m];
      // if(x3[m] < 0 || x3[m] > lx) cout << "x=" << x3[m] << "\n";
      z3[m] = z2[m];
      // if(z3[m] < 0 || z3[m] > lz) cout << "z=" << z3[m] << "\n";
      zeta3[m] = zeta2[m];
      // if(zeta3[m] < 0 || zeta3[m] > pi2) cout << "zeta3=" << zeta3[m] << "\n";
      u3[m] = u2[m];
//    w2(m) = w2(m)-myavgw
      w3[m] = w2[m]; 
   }
}

void gradu_c_(CArray3D<double> &u_, CArray3D<double> &ux_, CArray3D<double> &uz_){
   int ju = 0;
   int jl = 0;
   double ul = 0;

   for (int i = 0; i <= imx-1; ++i){ //changed order here in case anyone ever uses
      for(int j = 0; j < jmx; ++j){
         ju = j+1;
         jl = j-1;
         if(j == 0) jl = jmx-1;
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

//helper functions for integ parallelization
#pragma acc routine seq
inline double my_fmod(double a, double p) {
    double q = floor(a / p);
    double r = a - p * q;
    return (r < 0.0) ? r + p : r;
}

#pragma acc routine seq
inline size_t get4DIndex(size_t i, size_t j, size_t k, size_t l,
                  size_t ny, size_t nz, size_t nq) {
    return i * (ny * nz * nq) + j * (nz * nq) + k * nq + l;
}
#pragma acc routine seq
inline size_t get3DIndex(size_t i, size_t j, size_t k,
                  size_t ny, size_t nz) {
    return i * (ny * nz) + j * nz + k;
}

void integ_c_(int iflag) {
   int i = 0;
   int j = 0;
   int k = 0;
   int m = 0;
   int l = 0;

   double wx0 = 0;
   double wx1 = 0;
   double wzeta0 = 0;
   double wzeta1 = 0;
   double wy0 = 0;
   double wy1 = 0;
   double x = 0;
   double z = 0;
   double zeta = 0;
   double R_major_over_R = 0;
   double R_major_over_R1 = 0;

   double xt = 0,zt = 0;

   double wz0 = 0;
   double wz1 = 0;
   double bfldp = 0;
   double b = 0;
   double rhog = 0;
   double rhox[4], rhoy[4];

   auto start_integ_tm = MPI_Wtime();

   //using ii,jj,kk to distinguish between i,j, and k used in m loop - maybe add clear for 4d based on iflag - like clear 3d, but start at index based on iflag
   for(int ii = 0; ii <= imx; ++ii) {
      for (int jj = 0; jj <= jmx; ++jj) {
         for(int kk = 0; kk <= kmx; ++kk) {
            den(iflag, ii, jj, kk) = 0;
         }
      }
   }
   upar.Clear(); 
   
   // OPENACC "atomic update" is not compatible with user defined objects. So, we need to maneuver around this: 
   // auto grabs most relavant data type - since den is a double array, auto defines as double. Useful since if den was ever int (for sake of argument), this wouldn't break

   size_t idx = 0;                     // idx used to store index calculated for loop - used in den_ptr and upar_ptr 

   den.updatedev();                        
   upar.updatedev();
   auto* den_ptr = den.start();       
   auto* upar_ptr = upar.start(); 

   #pragma acc parallel loop gang vector present(den_ptr, upar_ptr, x3, z3, zeta3, u3, w3, gw, mu) private(rhox,rhoy) 
   for(m = 0; m < mm[0]; ++m) {
      const auto denY = den.getY();
      const auto denZ = den.getZ();
      const auto denQ = den.getQ();

      const double uparY = upar.getY();
      const double uparZ = upar.getZ();

      x = x3[m];
      if(x < 0 || x > lx) printf("integ x=%lf\n", x);
      i = static_cast<int>(x/dxeq);
      wx0 = (i+1)-x/dxeq;
      wx1 = 1.-wx0;

      z = z3[m];
      k = static_cast<int>(z/dzeq);
      if(z < 0 || z > lz) printf("integ z=%lf\n", z);
      wz0 = (k+1)-z/dzeq;
      wz1 = 1.-wz0;

      bfldp = wx0*wz0*b0(i,k) + wx0*wz1*b0(i,k+1) + wx1*wz0*b0(i+1,k) + wx1*wz1*b0(i+1,k+1);
      b = 1.-tor+tor*bfldp;

      rhog = sqrt(2.*b*mu[m]*mims[0])/(q[0]*b) * iflr;

      rhox[0] = rhog;
      rhoy[0] = 0;
      rhox[1] = -rhox[0];
      rhoy[1] = -rhoy[0];
      rhox[2] = 0;
      rhoy[2] = rhog;
      rhox[3] = 0;
      rhoy[3] = -rhoy[2];

      #pragma acc loop seq
      for(l = 0; l < lr[0]; ++l){
         xt=x3[m]+rhox[l];
         zt=z3[m]+rhoy[l];

         if( (xt<2*dxeq) || (xt>lx-2*dxeq) ) xt=x3[m];
         if( (zt<2*dzeq) || (zt>lz-2*dzeq) ) zt=z3[m];

         zeta= my_fmod(zeta3[m], pi2); 
         if(zeta < 0 || zeta > pi2)printf("integ zeta=%lf\n", zeta);

         i = static_cast<int>(xt/dx);
         j = static_cast<int>(zt/dz);
         k = static_cast<int>(zeta/dzeta);

         wx0 = (i+1)-xt/dx;
         wx1 = 1.-wx0;
         wy0 = (j+1)-zt/dz;
         wy1 = 1.-wy0;
         wzeta0 = (k+1)-zeta/dzeta;
         wzeta1 = 1.-wzeta0;

         R_major_over_R=xctr/(xctr-xdim/2+i*dx);
         R_major_over_R1=xctr/(xctr-xdim/2+(i+1)*dx);
        

         idx = get4DIndex(iflag,i,j,k, denY, denZ, denQ);
         #pragma acc atomic update
         den_ptr[idx] += gw[m]*w3[m]*wx0*wy0*wzeta0*R_major_over_R/4;

         idx = get4DIndex(iflag,i+1,j,k, denY, denZ, denQ);
         #pragma acc atomic update
         den_ptr[idx] += gw[m]*w3[m]*wx1*wy0*wzeta0*R_major_over_R1/4;

         idx = get4DIndex(iflag,i,j+1,k, denY, denZ, denQ);
         #pragma acc atomic update
         den_ptr[idx] += gw[m]*w3[m]*wx0*wy1*wzeta0*R_major_over_R/4;

         idx = get4DIndex(iflag,i+1,j+1,k, denY, denZ, denQ);
         #pragma acc atomic update
         den_ptr[idx] += gw[m]*w3[m]*wx1*wy1*wzeta0*R_major_over_R1/4;

         idx = get3DIndex(i,j,k, uparY, uparZ);
         #pragma acc atomic update 
         upar_ptr[idx] += u3[m]*gw[m]*w3[m]*wx0*wy0*wzeta0*R_major_over_R/4;

         idx = get3DIndex(i+1,j,k, uparY, uparZ);
         #pragma acc atomic update
         upar_ptr[idx] += u3[m]*gw[m]*w3[m]*wx1*wy0*wzeta0*R_major_over_R1/4;

         idx = get3DIndex(i,j+1,k, uparY, uparZ);
         #pragma acc atomic update
         upar_ptr[idx] += u3[m]*gw[m]*w3[m]*wx0*wy1*wzeta0*R_major_over_R/4;

         idx = get3DIndex(i+1,j+1,k, uparY, uparZ);
         #pragma acc atomic update
         upar_ptr[idx] += u3[m]*gw[m]*w3[m]*wx1*wy1*wzeta0*R_major_over_R1/4;

         if(k != kmx) {
            idx = get4DIndex(iflag, i, j, k+1, denY, denZ, denQ);
            #pragma acc atomic update
            den_ptr[idx] += gw[m]*w3[m]*wx0*wy0*wzeta1*R_major_over_R/4;

            idx = get4DIndex(iflag, i+1, j, k+1, denY, denZ, denQ);
            #pragma acc atomic update
            den_ptr[idx] += gw[m]*w3[m]*wx1*wy0*wzeta1*R_major_over_R1/4;

            idx = get4DIndex(iflag, i, j+1, k+1, denY, denZ, denQ);
            #pragma acc atomic update
            den_ptr[idx] += gw[m]*w3[m]*wx0*wy1*wzeta1*R_major_over_R/4;

            idx = get4DIndex(iflag, i+1, j+1, k+1, denY, denZ, denQ);
            #pragma acc atomic update
            den_ptr[idx] += gw[m]*w3[m]*wx1*wy1*wzeta1*R_major_over_R1/4;

            idx = get3DIndex(i, j, k+1, uparY, uparZ);
            #pragma acc atomic update
            upar_ptr[idx] += u3[m]*gw[m]*w3[m]*wx0*wy0*wzeta1*R_major_over_R/4;

            idx = get3DIndex(i+1, j, k+1, uparY, uparZ);
            #pragma acc atomic update
            upar_ptr[idx] += u3[m]*gw[m]*w3[m]*wx1*wy0*wzeta1*R_major_over_R1/4;

            idx = get3DIndex(i, j+1, k+1, uparY, uparZ);
            #pragma acc atomic update
            upar_ptr[idx] += u3[m]*gw[m]*w3[m]*wx0*wy1*wzeta1*R_major_over_R/4;

            idx = get3DIndex(i+1, j+1, k+1, uparY, uparZ);
            #pragma acc atomic update
            upar_ptr[idx] += u3[m]*gw[m]*w3[m]*wx1*wy1*wzeta1*R_major_over_R1/4;

         } else {
            idx = get4DIndex(iflag, i, j, 0, denY, denZ, denQ);
            #pragma acc atomic update
            den_ptr[idx] += gw[m]*w3[m]*wx0*wy0*wzeta1*R_major_over_R/4;

            idx = get4DIndex(iflag, i+1, j, 0, denY, denZ, denQ);
            #pragma acc atomic update
            den_ptr[idx] += gw[m]*w3[m]*wx1*wy0*wzeta1*R_major_over_R1/4;

            idx = get4DIndex(iflag, i, j+1, 0, denY, denZ, denQ);
            #pragma acc atomic update
            den_ptr[idx] += gw[m]*w3[m]*wx0*wy1*wzeta1*R_major_over_R/4;

            idx = get4DIndex(iflag, i+1, j+1, 0, denY, denZ, denQ);
            #pragma acc atomic update
            den_ptr[idx] += gw[m]*w3[m]*wx1*wy1*wzeta1*R_major_over_R1/4;

            idx = get3DIndex(i, j, 0, uparY, uparZ);
            #pragma acc atomic update
            upar_ptr[idx] += u3[m]*gw[m]*w3[m]*wx0*wy0*wzeta1*R_major_over_R/4;

            idx = get3DIndex(i+1, j, 0, uparY, uparZ);
            #pragma acc atomic update
            upar_ptr[idx] += u3[m]*gw[m]*w3[m]*wx1*wy0*wzeta1*R_major_over_R1/4;

            idx = get3DIndex(i, j+1, 0, uparY, uparZ);
            #pragma acc atomic update
            upar_ptr[idx] += u3[m]*gw[m]*w3[m]*wx0*wy1*wzeta1*R_major_over_R/4;

            idx = get3DIndex(i+1, j+1, 0, uparY, uparZ);
            #pragma acc atomic update
            upar_ptr[idx] += u3[m]*gw[m]*w3[m]*wx1*wy1*wzeta1*R_major_over_R1/4;
         }
      }   
   }
   
   den.updatehost();
   upar.updatehost();
   
   MPI_Allreduce(MPI_IN_PLACE, &den(iflag,0,0,0), (imx+1)*(jmx+1)*(kmx+1), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, upar.start(), (imx+1)*(jmx+1)*(kmx+1),MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
   for(int i = 0 ; i <= imx; ++i) {
      for(int j = 0; j <= jmx; ++j) {
         for(int k = 0; k <= kmx; ++k) {
            den(1,i,j,k) = den(iflag, i, j, k);   ///This is the den we want
         }
      }
   }

   den2d2.Clear();
   
   for(int i = 0; i <= imx; ++i) {
      for(int j = 0; j <= jmx; ++j) {
         for(int k = 0; k <= kmx; ++k) {
            den2d2(i,j)=den2d2(i,j)+den(iflag,i,j,k);
            if (i3D==0 && k != 0) upar(i,j,0) += upar(i,j,k);
         }
      }
   }

   for(int i = 0; i <= imx; ++i) {
      for(int j = 0; j <= jmx; ++j) {
         den2d2(i,j)=den2d2(i,j)/(kmx+1);
      }
   }

   for(i = 0; i <= imx; ++i) {
      for(j = 0; j <= jmx; ++j) {
         if(i3D == 0) {
            upar(i,j,0) = upar(i,j,0)/(kmx+1);
            for(k = 1; k <= kmx; ++k) {
               upar(i,j,k) = upar(i,j,0);
            }
         }
      }
   }

   if(iflag==1) {
      for(i = 0; i <= imx; ++i) {
         for(j = 0; j <= jmx; ++j) {
            dden2d(i,j)-=den2d1(i,j);
            den2d1(i,j)=den2d2(i,j);
            // for(k = 0 ; k <= kmx; ++k){
            //    //den_pre=den(1,i,j,k); 
            // }
         }
      }
   }

   
   auto end_integ_tm = MPI_Wtime();
   integ_tm = integ_tm + end_integ_tm - start_integ_tm;  //integ_tm + 
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
               apar(i,j,k) = apar(i,j,k) -dt*(gradPar(i,j,k));  //!+0.04*(jpar+q(1)*mu0*upar)) !+Epar)
            }
         }
      }
   }
}

void get_ne_c_(const int flagnumber){
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

//outdated - Needs update for CARRAY's
void gradparz_c_(double *matrix){
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
   for(int i = 2; i < imx-1; ++i) {
      for(int j = 2; j < jmx-1; ++j) {
         for(int k = 1; k < kmx; ++k) {
            if (!(mask2(i,j)<1.99)) {
               gradPar(i,j,k)=(b0x(i,j)/b0(i,j)*(matrix(i+1,j,k)-matrix(i-1,j,k))*0.5/dx  
               +b0z(i,j)/b0(i,j)*(matrix(i,j+1,k)-matrix(i,j-1,k))*0.5/dz                
               +b0zeta(i,j)/b0(i,j)*(matrix(i,j,k+1)-matrix(i,j,k-1))/((Rgrid[i]/xu)*2*dzeta));
            }
         }
      }
   }

   for(int i = 2; i < imx-1; ++i){
      for(int j = 2; j < jmx-1; ++j){
            if(!(mask2(i,j)<1.99)) {
               gradPar(i,j,0)=(b0x(i,j)/b0(i,j)*(matrix(i+1,j,0)-matrix(i-1,j,0))*0.5/dx  
               +b0z(i,j)/b0(i,j)*(matrix(i,j+1,0)-matrix(i,j-1,0))*0.5/dz                
               +b0zeta(i,j)/b0(i,j)*(matrix(i,j,1)-matrix(i,j,kmx))/((Rgrid[i]/xu)*2*dzeta));
         }
      }
   }

   for(int i = 2; i < imx-1; ++i){
      for(int j = 2; j < jmx-1; ++j){
         if(!(mask4(i,j)<3.99)){
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
	PetscInt i,j,mx,my,xm;
	PetscInt    ym,xs,ys,i1, i5; 
	PetscScalar  v[5],Hx,Hy;
   // PetscScalar v[9],Hx,Hy;
	PetscScalar  Hx2,Hy2; //tmp_r,a_value
	MatStencil   row,col[5];
   // MatStencil   row,col[9]; 
   double rho_squared = 0;
	i1 = 1;
	i5 = 5;
   // i5 = 9;
	//a_value = 0.5;
   // auto start_t = MPI_Wtime();
	PetscCall(KSPGetDM(ksp,&dm)); 
	PetscCall(DMDAGetInfo(dm,nullptr,&mx,&my,nullptr,nullptr,nullptr,nullptr,
						nullptr,nullptr,nullptr,nullptr,nullptr,nullptr)); 

	Hx = dx;//! (Rgrid(imx)-Rgrid(0)) / real(imx)
	Hy = dz;//(Zgrid(jmx)-Zgrid(0)) / real(jmx)

	Hx2 = Hx*Hx;
	Hy2 = Hy*Hy;
	PetscCall(DMDAGetCorners(dm,&xs,&ys,nullptr,&xm,&ym,nullptr));
	for(j=ys; j < ys+ym; ++j){
		for(i=xs; i < xs+xm; ++i) {
         rho_squared = rho_i(i,j) * rho_i(i,j);
			row.i = i;
			row.j = j;
			if(mask(i,j) < 0.99){
				v[0] = (c2_over_vA2(i,j) + PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_squared))*(-2.0/Hx2-2.0/Hy2);
            // v[0] = 1.0;
				PetscCall(MatSetValuesStencil(BB,i1,&row,i1,&row,&v[0],INSERT_VALUES));
			} else {
				if(j > 0) {
					if(j == jmx) {
						v[0] = (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_squared))/Hy2-1.0/(2.0*Hy2)*( c2_over_vA2(i,j)- c2_over_vA2(i,j-1));
						v[0] = v[0] - PADE*(mu0*e*e*rho_squared)*(xn0e(i,j)/t0e(i,j) - xn0e(i,j-1)/t0e(i,j-1))/Hy2;
					} else {
						v[0] = (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_squared))/Hy2-1.0/(4.0*Hy2)*( c2_over_vA2(i,j+1)- c2_over_vA2(i,j-1));
                  v[0] = v[0] - PADE*(mu0*e*e*rho_squared)*(xn0e(i,j+1)/t0e(i,j+1) - xn0e(i,j-1)/t0e(i,j-1))/(2*Hy2);
					}
				}
				col[0].i = i;
				col[0].j = j - 1;

				if(i > 0) {
					if(i == imx) {
						v[1] =  (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_squared))/Hx2-1.0/(2.0*Hx2)*( c2_over_vA2(i,j)- c2_over_vA2(i-1,j));
                  v[1] = v[1] - PADE*(mu0*e*e*rho_squared)*(xn0e(i,j)/t0e(i,j) - xn0e(i-1,j)/t0e(i-1,j))/Hx2;
					} else {
						v[1] =  (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_squared))/Hx2-1.0/(4.0*Hx2)*( c2_over_vA2(i+1,j)- c2_over_vA2(i-1,j));
                  v[1] = v[1] - PADE*(mu0*e*e*rho_squared)*(xn0e(i+1,j)/t0e(i+1,j) - xn0e(i-1,j)/t0e(i-1,j))/(2*Hx2);
					}
				}
				col[1].i = i - 1;
				col[1].j = j;

				v[2] = -2.0*(c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_squared))/Hx2 - 2.0*(c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_squared))/Hy2;
				col[2].i = i;
				col[2].j = j;

				//cout << v[2] << "		" << xn0e(i,j)*mu0*e/t0e(i,j) << endl;
	//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Boltzmann e!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
				if(iBoltzmann != 0) {
					v[2]=v[2]-xn0e(i,j)*mu0*e*e/t0e(i,j) + PADE*(mu0*e*e*rho_squared)*((xn0e(i-1,j)/t0e(i-1,j) + xn0e(i+1,j)/t0e(i+1,j) - 2*xn0e(i,j)/t0e(i,j))/Hx2 + 
					(xn0e(i,j-1)/t0e(i,j-1) + xn0e(i,j+1)/t0e(i,j+1) - 2*xn0e(i,j)/t0e(i,j))/Hy2);
				}

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

				if(i < imx){
					if(i == 0){
						v[3] =  (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_squared))/Hx2+1.0/(2.0*Hx2)*( c2_over_vA2(i+1,j)- c2_over_vA2(i,j));
                  v[3] = v[3] + PADE*(mu0*e*e*rho_squared)*(xn0e(i+1,j)/t0e(i+1,j) - xn0e(i,j)/t0e(i,j))/Hx2;
					} else {
						v[3] =  (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_squared))/Hx2+1.0/(4.0*Hx2)*( c2_over_vA2(i+1,j)- c2_over_vA2(i-1,j));
                  v[3] = v[3] + PADE*(mu0*e*e*rho_squared)*(xn0e(i+1,j)/t0e(i+1,j) - xn0e(i-1,j)/t0e(i-1,j))/(2*Hx2);
					}
				}
				col[3].i = i + 1;
				col[3].j = j;
				
				if(j < jmx) {
					if(j == 0){
						v[4] =  (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_squared))/Hy2+1.0/(2.0*Hy2)*( c2_over_vA2(i,j+1)- c2_over_vA2(i,j));
                  v[4] = v[4] + PADE*(mu0*e*e*rho_squared)*(xn0e(i,j+1)/t0e(i,j+1) - xn0e(i,j)/t0e(i,j))/Hy2;
					} else {
						v[4] =  (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_squared))/Hy2+1.0/(4.0*Hy2)*( c2_over_vA2(i,j+1)- c2_over_vA2(i,j-1));
                  v[4] = v[4] + PADE*(mu0*e*e*rho_squared)*(xn0e(i,j+1)/t0e(i,j+1) - xn0e(i,j-1)/t0e(i,j-1))/(2*Hy2);
					} 
				}
				col[4].i = i;
				col[4].j = j + 1;   

            // Implementation of the 9-point stencil
            // v[5] = (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_squared))/(Hx*Hy) - 1.0/(4.0*(Hx*Hy))*(c2_over_vA2(i+1,j+1)- c2_over_vA2(i-1,j-1));
            // v[5] = v[5] - PADE*(mu0*e*e*rho_squared)*(xn0e(i+1,j+1)/t0e(i+1,j+1) - xn0e(i-1,j-1)/t0e(i-1,j-1))/(2*(Hx*Hy));

            // col[5].i = i-1;
            // col[5].j = j+1;

            // v[6] = (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_squared))/(Hx*Hy) + 1.0/(4.0*(Hx*Hy))*(c2_over_vA2(i+1,j+1)- c2_over_vA2(i-1,j-1));
            // v[6] = v[6] + PADE*(mu0*e*e*rho_squared)*(xn0e(i+1,j+1)/t0e(i+1,j+1) - xn0e(i-1,j-1)/t0e(i-1,j-1))/(2*(Hx*Hy));

            // col[6].i = i+1;
            // col[6].j = j+1;
            
            // v[7] = (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_squared))/(Hx*Hy) + 1.0/(4.0*(Hx*Hy))*(c2_over_vA2(i+1,j+1)- c2_over_vA2(i-1,j-1));
            // v[7] = v[7] + PADE*(mu0*e*e*rho_squared)*(xn0e(i+1,j+1)/t0e(i+1,j+1) - xn0e(i-1,j-1)/t0e(i-1,j-1))/(2*(Hx*Hy));

            // col[7].i = i-1;
            // col[7].j = j-1;

            // v[8] = (c2_over_vA2(i,j)+PADE*((xn0e(i,j)*mu0*e*e/t0e(i,j))*rho_squared))/(Hx*Hy) - 1.0/(4.0*(Hx*Hy))*(c2_over_vA2(i+1,j+1)- c2_over_vA2(i-1,j-1));
            // v[8] = v[8] - PADE*(mu0*e*e*rho_squared)*(xn0e(i+1,j+1)/t0e(i+1,j+1) - xn0e(i-1,j-1)/t0e(i-1,j-1))/(2*(Hx*Hy));

            // col[8].i = i+1;
            // col[8].j = j-1;

            // v[0] = 0.5*v[0];
            // v[1] = 0.5*v[1];
            // v[2] = 0.75*v[2];
            // v[3] = 0.5*v[3];
            // v[4] = 0.5*v[4];
            // v[5] = 0.25*v[5];
            // v[6] = 0.25*v[6];
            // v[7] = 0.25*v[7];
            // v[8] = 0.25*v[8];

				PetscCall(MatSetValuesStencil(BB, i1, &row, i5, col, v, INSERT_VALUES));
			}
		}
   }

	PetscCall(MatAssemblyBegin(BB,MAT_FINAL_ASSEMBLY));
	PetscCall(MatAssemblyEnd(BB,MAT_FINAL_ASSEMBLY));
	if(AA != BB) {
		PetscCall(MatAssemblyBegin(AA,MAT_FINAL_ASSEMBLY));
		PetscCall(MatAssemblyEnd(AA,MAT_FINAL_ASSEMBLY));
	}
	//   PetscCall(MatView(AA,PETSC_VIEWER_STDOUT_WORLD));
	//   PetscCall(MatView(BB,PETSC_VIEWER_STDOUT_WORLD));
   // auto end_t = MPI_Wtime();
   // cout << "ComputeMatrixTm = " << end_t-start_t << endl;
	return 0;
}

PetscErrorCode ComputeRHS(KSP ksp, Vec bbb, void *ctx) {
   int ii,jj,iflag;
   //PetscScalar* b_array = new PetscScalar[]; 
   PetscScalar  h,Hx,Hy;
   PetscInt  mx,my,i,j,xs,xm,ys,ym,vec_start,vec_end;
   DM dm;
   PetscInt idx;
   PetscScalar tmp_value = 0.0;
   double rho_squared = 0;
   // PetscScalar a_value,tmp_r;

   PetscInt k = *(PetscInt*)ctx; 
   // auto start_t = MPI_Wtime();
   tmp_value = 0;
   
   PetscCall(KSPGetDM(ksp,&dm));
   PetscCall(DMDAGetInfo(dm,nullptr,&mx,&my,nullptr,nullptr,nullptr,
                        nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr));
   PetscCall(VecGetOwnershipRange(bbb,&vec_start,&vec_end));
   PetscCall(DMDAGetCorners(dm,&xs,&ys,nullptr,&xm,&ym,nullptr));

   idx = vec_start-1;
   for(j = ys; j < ys+ym; ++j) {
      for(i = xs; i < xs+xm; ++i) {
         rho_squared = rho_i(i,j) * rho_i(i,j);
         idx+=1;
         if(mask(i,j) < 0.99) {
            tmp_value = 0;
         } else {
			if(weightscheme == 0) {
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!2D ni noly now!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				if(i3D == 0) {
					if(iBoltzmann == 0){
						tmp_value = denes(i,j,k) - q[0]*mu0*(den2d2(i,j) - xn0i(i,j));
					} else if(eAdiabatic != 0) {
							tmp_value = -q[0]*mu0*(den2d2(i,j)-xn0i(i,j)) - (xn0e(i,j)*mu0*e*e/t0e(i,j))*phiavg(i,j) +  
										PADE*(mu0*q[0]*rho_squared)*(((den2d2(i-1,j)-xn0i(i-1,j))+(den2d2(i+1,j)-xn0i(i+1,j))-2*(den2d2(i,j)-xn0i(i,j)))/(dx*dx) + 
										((den2d2(i,j-1)-xn0i(i,j-1))+(den2d2(i,j+1)-xn0i(i,j+1))-2*(den2d2(i,j)-xn0i(i,j)))/(dz*dz)) + PADE*(e*e*mu0*rho_squared)*(phiavg(i,j)*(((xn0e(i-1,j)/t0e(i-1,j)) + 
										(xn0e(i+1,j)/t0e(i+1,j))-2*(xn0e(i,j)/t0e(i,j)))/(dx*dx) + ((xn0e(i,j-1)/t0e(i,j-1))+(xn0e(i,j+1)/t0e(i,j+1))-2*(xn0e(i,j)/t0e(i,j)))/(dz*dz)) +  
										(xn0e(i,j)/t0e(i,j))*((phiavg(i-1,j) + phiavg(i+1,j) -2*phiavg(i,j))/(dx*dx) + (phiavg(i,j-1)+phiavg(i,j+1) -2*phiavg(i,j))/(dz*dz)) + 
										(((xn0e(i+1,j)/t0e(i+1,j))-(xn0e(i-1,j)/t0e(i-1,j)))*(phiavg(i+1,j)-phiavg(i-1,j))/(2*(dx*dx)) + 
										((xn0e(i,j+1)/t0e(i,j+1))-(xn0e(i,j-1)/t0e(i,j-1)))*(phiavg(i,j+1)-phiavg(i,j-1))/(2*(dz*dz))));
					} else {
							tmp_value = -q[0]*mu0*(den2d2(i,j)-xn0i(i,j));
					}
//                tmp_value = 1
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3D case!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            	} else {
					if(iBoltzmann == 0) {
						tmp_value = denes(i,j,k)-q[0]*mu0*(den(1,i,j,k)-xn0i(i,j));

						} else if(eAdiabatic != 0) {
						// tmp_value = -q[0]*mu0*(den(1,i,j,k)-xn0i(i,j)) - (xn0e(i,j)*mu0*e*e/t0e(i,j))*phiavg(i,j) + 
						// 			PADE*(mu0*q[0]*rho_squared)*(((den(1,i-1,j,k)-xn0i(i-1,j))+(den(1,i+1,j,k)-xn0i(i+1,j))-2*(den(1,i,j,k)-xn0i(i,j)))/(dx*dx) + 
						// 			((den(1,i,j-1,k)-xn0i(i,j-1))+(den(1,i,j+1,k)-xn0i(i,j+1))-2*(den(1,i,j,k)-xn0i(i,j)))/(dz*dz)) + PADE*(e*e*mu0*rho_squared)*(phiavg(i,j)*(((xn0e(i-1,j)/t0e(i-1,j)) + 
						// 			(xn0e(i+1,j)/t0e(i+1,j))-2*(xn0e(i,j)/t0e(i,j)))/(dx*dx) + ((xn0e(i,j-1)/t0e(i,j-1))+(xn0e(i,j+1)/t0e(i,j+1))-2*(xn0e(i,j)/t0e(i,j)))/(dz*dz)) + 
						// 			(xn0e(i,j)/t0e(i,j))*((phiavg(i-1,j) + phiavg(i+1,j) -2*phiavg(i,j))/(dx*dx) + (phiavg(i,j-1)+phiavg(i,j+1) -2*phiavg(i,j))/(dz*dz)) + 
						// 			(((xn0e(i+1,j)/t0e(i+1,j))-(xn0e(i-1,j)/t0e(i-1,j)))*(phiavg(i+1,j)-phiavg(i-1,j))/(2*(dx*dx)) + 
						// 			((xn0e(i,j+1)/t0e(i,j+1))-(xn0e(i,j-1)/t0e(i,j-1)))*(phiavg(i,j+1)-phiavg(i,j-1))/(2*(dz*dz))));
                  tmp_value = -q[0]*mu0*(den(1,i,j,k)-xn0i(i,j)) - (xn0e(i,j)*mu0*e*e/t0e(i,j))*phiavg(i,j) + 
									PADE*(mu0*q[0]*rho_squared)*(((den(1,i-1,j,k)-xn0i(i-1,j))+(den(1,i+1,j,k)-xn0i(i+1,j))-2*(den(1,i,j,k)-xn0i(i,j)))/(dx*dx) + 
									((den(1,i,j-1,k)-xn0i(i,j-1))+(den(1,i,j+1,k)-xn0i(i,j+1))-2*(den(1,i,j,k)-xn0i(i,j)))/(dz*dz)) + PADE*(e*e*mu0*rho_squared) *
                           (((xn0e(i+1,j)*phiavg(i+1,j)/t0e(i+1,j))+(xn0e(i-1,j)*phiavg(i-1,j)/t0e(i-1,j))-(2*(xn0e(i,j)*phiavg(i,j)/t0e(i,j))))/(dx*dx) + 
                           ((xn0e(i,j+1)*phiavg(i,j+1)/t0e(i,j+1))+(xn0e(i,j-1)*phiavg(i,j-1)/t0e(i,j-1))-(2*(xn0e(i,j)*phiavg(i,j)/t0e(i,j))))/(dz*dz));
                        
					} else {
						tmp_value = -q[0]*mu0*(den(1,i,j,k)-xn0i(i,j));
					}
				}
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            } else {
				//weightscheme on use this!
				if(i3D==0) {
					if(iBoltzmann == 0) {
						tmp_value = denes(i,j,k)-q[0]*mu0*(den2d2(i,j));
					} else if(eAdiabatic != 0) {
						tmp_value = -q[0]*mu0*(den2d2(i,j)) - (xn0e(i,j)*mu0*e*e/t0e(i,j))*phiavg(i,j) +  
									PADE*(mu0*q[0]*rho_squared)*((den2d2(i-1,j)+den2d2(i+1,j)-2*den2d2(i,j))/(dx*dx) + 
									(den2d2(i,j-1)+den2d2(i,j+1)-2*den2d2(i,j))/(dz*dz)) + PADE*(e*e*mu0*rho_squared)*(phiavg(i,j)*(((xn0e(i-1,j)/t0e(i-1,j)) + 
									(xn0e(i+1,j)/t0e(i+1,j))-2*(xn0e(i,j)/t0e(i,j)))/(dx*dx) + ((xn0e(i,j-1)/t0e(i,j-1))+(xn0e(i,j+1)/t0e(i,j+1))-2*(xn0e(i,j)/t0e(i,j)))/(dz*dz)) + 
									(xn0e(i,j)/t0e(i,j))*((phiavg(i-1,j) + phiavg(i+1,j) -2*phiavg(i,j))/(dx*dx) + (phiavg(i,j-1)+phiavg(i,j+1) -2*phiavg(i,j))/(dz*dz)) + 
									(((xn0e(i+1,j)/t0e(i+1,j))-(xn0e(i-1,j)/t0e(i-1,j)))*(phiavg(i+1,j)-phiavg(i-1,j))/(2*(dx*dx)) + 
									((xn0e(i,j+1)/t0e(i,j+1))-(xn0e(i,j-1)/t0e(i,j-1)))*(phiavg(i,j+1)-phiavg(i,j-1))/(2*(dz*dz))));
					} else {
						tmp_value = -q[0]*mu0*(den2d2(i,j));
					}
				} else {

					if(iBoltzmann == 0) {
						tmp_value = denes(i,j,k)-q[0]*mu0*(den(1,i,j,k));
					} else if (eAdiabatic != 0) {
						// tmp_value = -q[0]*mu0*(den(1,i,j,k)) - (xn0e(i,j)*mu0*e*e/t0e(i,j))*phiavg(i,j) + 
						// 			PADE*(mu0*q[0]*rho_squared)*((den(1,i-1,j,k)+den(1,i+1,j,k)-2*den(1,i,j,k))/(dx*dx) + 
						// 			(den(1,i,j-1,k)+den(1,i,j+1,k)-2*den(1,i,j,k))/(dz*dz)) + PADE*(e*e*mu0*rho_squared)*(phiavg(i,j)*(((xn0e(i-1,j)/t0e(i-1,j)) + 
						// 			(xn0e(i+1,j)/t0e(i+1,j))-2*(xn0e(i,j)/t0e(i,j)))/(dx*dx) + ((xn0e(i,j-1)/t0e(i,j-1))+(xn0e(i,j+1)/t0e(i,j+1))-2*(xn0e(i,j)/t0e(i,j)))/(dz*dz)) +  
						// 			(xn0e(i,j)/t0e(i,j))*((phiavg(i-1,j) + phiavg(i+1,j) -2*phiavg(i,j))/(dx*dx) + (phiavg(i,j-1)+phiavg(i,j+1) -2*phiavg(i,j))/(dz*dz)) + 
						// 			(((xn0e(i+1,j)/t0e(i+1,j))-(xn0e(i-1,j)/t0e(i-1,j)))*(phiavg(i+1,j)-phiavg(i-1,j))/(2*(dx*dx)) + 
						// 			((xn0e(i,j+1)/t0e(i,j+1))-(xn0e(i,j-1)/t0e(i,j-1)))*(phiavg(i,j+1)-phiavg(i,j-1))/(2*(dz*dz))));

                  tmp_value = -q[0]*mu0*(den(1,i,j,k)) - (xn0e(i,j)*mu0*e*e/t0e(i,j))*phiavg(i,j) + 
									PADE*(mu0*q[0]*rho_squared)*(((den(1,i-1,j,k))+(den(1,i+1,j,k))-2*(den(1,i,j,k)))/(dx*dx) + 
									((den(1,i,j-1,k))+(den(1,i,j+1,k))-2*(den(1,i,j,k)))/(dz*dz)) + PADE*(e*e*mu0*rho_squared) *
                           (((xn0e(i+1,j)*phiavg(i+1,j)/t0e(i+1,j))+(xn0e(i-1,j)*phiavg(i-1,j)/t0e(i-1,j))-(2*(xn0e(i,j)*phiavg(i,j)/t0e(i,j))))/(dx*dx) + 
                           ((xn0e(i,j+1)*phiavg(i,j+1)/t0e(i,j+1))+(xn0e(i,j-1)*phiavg(i,j-1)/t0e(i,j-1))-(2*(xn0e(i,j)*phiavg(i,j)/t0e(i,j))))/(dz*dz));

                     // if(i == 100 && j == 100) printf("%e\n", tmp_value);

                    //  ! write(*,*) -q(1)*mu0*(den(2,i,j,k)), -(xn0e(i,j)*mu0*e*e/t0e(i,j))*phiavg(i,j), PADE*(mu0*q(1)*rho_i(i,j)**2)*((den(2,i-1,j,k)+den(2,i+1,j,k)-2*den(2,i,j,k))/dx**2 + (den(2,i,j-1,k)+den(2,i,j+1,k)-2*den(2,i,j,k))/dz**2), PADE*(e*e*mu0*rho_i(i,j)**2)*(phiavg(i,j)*(((xn0e(i-1,j)/t0e(i-1,j)) + &
                    //  ! (xn0e(i+1,j)/t0e(i+1,j))-2*(xn0e(i,j)/t0e(i,j)))/dx**2 + ((xn0e(i,j-1)/t0e(i,j-1))+(xn0e(i,j+1)/t0e(i,j+1))-2*(xn0e(i,j)/t0e(i,j)))/dz**2) + & 
                    //  ! (xn0e(i,j)/t0e(i,j))*((phiavg(i-1,j) + phiavg(i+1,j) -2*phiavg(i,j))/dx**2 + (phiavg(i,j-1)+phiavg(i,j+1) -2*phiavg(i,j))/dz**2) + &
                    //  ! (((xn0e(i+1,j)/t0e(i+1,j))-(xn0e(i-1,j)/t0e(i-1,j)))*(phiavg(i+1,j)-phiavg(i-1,j))/(2*dx**2) + &
                    //  ! ((xn0e(i,j+1)/t0e(i,j+1))-(xn0e(i,j-1)/t0e(i,j-1)))*(phiavg(i,j+1)-phiavg(i,j-1))/(2*dz**2)))
					} else {
						tmp_value = -q[0]*mu0*(den(1,i,j,k));
					}
				}
			}
      }
         PetscCall(VecSetValues(bbb,1,&idx, &tmp_value, INSERT_VALUES));
   }
}
   PetscCall(VecAssemblyBegin(bbb));
   PetscCall(VecAssemblyEnd(bbb));
   // auto end_t = MPI_Wtime();
   // cout << "ComputeRHS time = " << end_t-start_t << endl;
   return(PETSC_SUCCESS);
}

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALDER Flux Average SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
// void fluxavg_c_(CArray3D<double> &input, CArray2D<double> &output){
//    //Local Variables
//    double phiavg1d[100];
//    double psi1d[100];
//    double psi1d_private[100];
//    double phiavg1d_private[100];

//    int line_large, line_small, line_marker;
//    int gi = 0, xix=0, yjy=0, miw=0, psi_zero=0, store=0, k=0, priv_mark; 
//    double weightinput=0,weightinput3D=0, phiavggi=0, psival=0, wmx0=0, wmx1=0, psi_private_min = 0.0;

//    // set all arrays to zero here. 
//    std::fill(std::begin(phiavg1d), std::end(phiavg1d), 0.0);
//    std::fill(std::begin(psi1d), std::end(psi1d), 0.0);
//    std::fill(std::begin(psi1d_private), std::end(psi1d_private), 0.0);
//    std::fill(std::begin(phiavg1d_private), std::end(phiavg1d_private), 0.0);

//    line_small = 0;
//    psi_zero = 1;

//    //REVISED COMPUTATION
//    line_marker = 0;

//    for(gi = 0; gi < 100; ++gi) {
//       weightinput = 0.0;
//       phiavggi = 0.0;
//       store = 0;
//       priv_mark = 0;
//       for(line = line_marker; line <=  num_lines; ++line) {
//          if(gindex[line] == gi) {
//             weightinput = (weight00[line]*input(iarray[line], jarray[line],0) + 
//                           weight10[line]*input(iarray[line]+1, jarray[line],0) + 
//                           weight01[line]*input(iarray[line], jarray[line]+1,0) + 
//                           weight11[line]*input(iarray[line]+1, jarray[line]+1,0));
//             if(i3D == 0) {
//                phiavggi += (weightinput*jacobian[line])/deno[line];
//             } else {
//                for(k = 1; k <= kmx; ++k) {
//                   weightinput3D = (weight00[line]*input(iarray[line], jarray[line],k) + 
//                                    weight10[line]*input(iarray[line]+1, jarray[line],k) + 
//                                    weight01[line]*input(iarray[line], jarray[line]+1,k) + 
//                                    weight11[line]*input(iarray[line]+1, jarray[line]+1,k));
//                   weightinput += weightinput3D;
//                }
//                phiavggi = phiavggi +(weightinput*jacobian[line])/(deno[line]*(kmx+1));
//             }
//             store = line;
//             if(priv[line] == 0) {
//                priv_mark = 1;
//             }
//          }
//          //Remove redundancy from closed loop integration process
//          if(phiavggi != 0) {
//             if(gindex[line] != gi) {
//                if(i3D == 0) {
//                   phiavggi = phiavggi - (weightinput*jacobian[line-1])/deno[line-1];
//                } else {
//                   phiavggi = phiavggi - (weightinput*jacobian[line-1])/(deno[line-1]*(kmx+1));
//                }
//                line_marker = line;
//                break;
//             }
//          }
//       }

//       if(priv_mark != 0) {
//         phiavg1d[psi_zero] = phiavggi;
//         psi1d[psi_zero]    = psitab[store];
//         psi_zero = psi_zero + 1;
//       } else {
//          phiavg1d_private[gi] = phiavggi;
//          if(line_small == 0) {
//             psi_private_min = psitab[store];
//          }
//          line_small += 1;
//       }
//    }

//    //Initialize output to zero
//    phiavg1d[0] = phiavg1d[1];
//    output.Clear();
//    //QUICK FIX
//    //phiavg1d_private = 0;

//    //INTERPOLATION
//    for(xix = 0; xix <= nx; ++xix) {
//       for(yjy = 0; yjy <= nz; ++yjy) {
//          psival = psi_p(xix,yjy);
//          if(mask(xix,yjy) < 0.99) {
//             output(xix,yjy) = 0;
//          } else {
//             //    if (yjy < 75 .and. xix < 150 .and. psival > 0.29 .and. psival<0.31) { //!Private region under X-point
//             //       miw  = int((psival-psi_private_min)/(psi1d(2)-psi1d(1)))
//             //     !   wmx0 = ((miw+1)*(psi1d_private(2)-psi1d_private(1))-psival)/(psi1d(2)-psi1d(1))
//             //       wmx0 = ((miw+1)*(psi1d(2)-psi1d(1))-psival)/(psi1d(2)-psi1d(1))
//             //       wmx1 = 1.-wmx0
//             //       output(xix,yjy) = wmx0*phiavg1d_private(miw) + wmx1*phiavg1d_private(miw+1)
//             //  } else {  
//             miw  = static_cast<int>(psival/(psi1d[2]-psi1d[1]));
//             wmx0 = ((miw+1)*(psi1d[2]-psi1d[1])-psival)/(psi1d[2]-psi1d[1]);
//             wmx1 = 1. - wmx0;
//             output(xix,yjy) = wmx0*phiavg1d[miw] + wmx1*phiavg1d[miw+1];
//          }
//       }
//    }
//    // output.Clear();
// }

void fluxavg_c_(CArray3D<double> &input, CArray2D<double> &output){

   int N = 200; //Flux surface resolution; 

   double tmp_num = 0.0;
   double tmp_den = 0.0;
   double max_val = 0.0;
   double val = 0.0;

   double wmx0=0.0, wmx1=0.0, psival=0.0;

   int i,j,k,l,miw;

   // std::vector<double> psi_1d(200);

   double tmp_flux_surface[200];
   double psi_1d[200];

   std::fill(std::begin(tmp_flux_surface), std::end(tmp_flux_surface), 0.0);
   std::fill(std::begin(psi_1d), std::end(psi_1d), 0.0);
   // double max_val = -std::numeric_limits<double>::infinity();

   // Create 1d psi array
   for (i = 0; i<=imx; ++i){
      for (j = 0; j<=jmx; ++j){
         val = psi_p(i,j);
         if (val > max_val){
            max_val = val;
         }
      }
   }

   double dpsi = max_val / (N-1);

   for(l = 1; l <= N; ++l){
      // psi_1d[l] = l*dpsi;
      for(k=0; k<=kmx; ++k){
         for(i=0; i<=imx; ++i){
            for(j=0; j<=jmx; ++j){
               if(psi_p(i,j) <= l*dpsi && psi_p(i,j) >= (l-1)*dpsi){
                  tmp_num += input(i,j,k)*Rgrid[i];
                  tmp_den += Rgrid[i];
               }
            }
         }
      }
      tmp_flux_surface[l] = tmp_num/tmp_den;
      psi_1d[l] = dpsi*l;
      tmp_num = 0.0;
      tmp_den = 0.0;
   }

   for(i=0; i<=imx; ++i){
      for(j=0; j<=jmx; ++j){
         psival = psi_p(i,j);
         if(mask(i,j) < 0.99){
            output(i,j) = 0.0;
         } else{
            miw = static_cast<int>(psival/(psi_1d[2]-psi_1d[1]));
            wmx0 = ((miw+1)*(psi_1d[2]-psi_1d[1])-psival)/(psi_1d[2]-psi_1d[1]);
            wmx1 = 1. - wmx0;
            output(i,j) = wmx0*tmp_flux_surface[miw] + wmx1*tmp_flux_surface[miw+1];
         }
      }
   }
}

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALDER E FIELD SUBROUTINE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void efieldcalc_c_(CArray3D<double> &input_phi){ 
   //local vars
   int i, j, k, kminus, kplus;
   for(i = 2; i < imx; ++i){
      for(j = 2; j < jmx; ++j){
         for(k = 0; k <= kmx; ++k) {
            ex(i,j,k) = -(input_phi(i+1,j,k) - input_phi(i-1,j,k))/(2*(Rgrid[1]-Rgrid[0]));
            ez(i,j,k) = -(input_phi(i,j+1,k) - input_phi(i,j-1,k))/(2*(Zgrid[1]-Zgrid[0]));

            // ex(i,j,k) = -(input_phi(i+1,j,k) - input_phi(i-1,j,k))/((Rgrid[i+1]-Rgrid[i-1]));
            // ez(i,j,k) = -(input_phi(i,j+1,k) - input_phi(i,j-1,k))/((Zgrid[j+1]-Zgrid[j-1]));
            if(k == 0){
               kminus = kmx;
               ezeta(i,j,k) = -(input_phi(i,j,k+1) - input_phi(i,j,kminus))/(2*Rgrid[i]*(2*pi/(kmx+1)));
            }
            else if(k == kmx){
               kplus = 0;
               ezeta(i,j,k) = -(input_phi(i,j,kplus) - input_phi(i,j,k-1))/(2*Rgrid[i]*(2*pi/(kmx+1)));
            } else {
               ezeta(i,j,k) = -(input_phi(i,j,k+1) - input_phi(i,j,k-1))/(2*Rgrid[i]*(2*pi/(kmx+1)));
            }

            if(cold_start == 1){
               if (timestep < 1000){
                  ex(i,j,k)=0.0;
                  ez(i,j,k)=0.0;
                  ezeta(i,j,k)=0.0;
               }
               if (timestep <= 2000){
                  ex(i,j,k)=0.5*(1.0-cos(pi*(timestep-1000)/1000.0))*ex(i,j,k);
                  ez(i,j,k)=0.5*(1.0-cos(pi*(timestep-1000)/1000.0))*ez(i,j,k);
                  ezeta(i,j,k)=0.5*(1.0-cos(pi*(timestep-1000)/1000.0))*ezeta(i,j,k);
               }
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

   for(gi = 21759; gi <= 22539; ++gi){ //all 1d arrays in ftn 1-indexed. 
      if(gindex[gi] == peak){
         for(k = 0; k <= kmx; ++k){
            weightinput = (weight00[gi]*input_phi(iarray[gi], jarray[gi],k)*input_phi(iarray[gi], jarray[gi],k) + 
                           weight10[gi]*input_phi(iarray[gi]+1, jarray[gi],k)*input_phi(iarray[gi]+1, jarray[gi],k) + 
                           weight01[gi]*input_phi(iarray[gi], jarray[gi]+1,k)*input_phi(iarray[gi], jarray[gi]+1,k) + 
                           weight11[gi]*input_phi(iarray[gi]+1, jarray[gi]+1,k)*input_phi(iarray[gi]+1, jarray[gi]+1,k));
         }
         phiavggi += (weightinput*jacobian[gi])/(deno[gi]*(kmx+1));

         if(priv[gi] == 0){
            store = gi;
         }
      }

      if(phiavggi != 0){
         if(gindex[gi] != peak){
            phiavggi -= (weightinput*jacobian[gi-1])/(deno[gi-1]*(kmx+1));
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
      myFile << std::setprecision(16) << phiavgsq << "\n";
   }else{
      std::cerr << "Error opening testphiavgsq" << "\n";
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
            } else {
               phi_k(i,j,k) = 1.001*input_phi(i,j,k);
            }
         }
      }
   }

   for(i = 2; i < imx; ++i){
      for(j = 2; j < jmx; ++j){
         for(k = 0; k <= kmx; ++k){
            dphidr(i,j,k)   = Rgrid[i]*(input_phi(i+1,j,k)-input_phi(i-1,j,k))/(2*dx);
            dphi_kdr(i,j,k) = Rgrid[i]*(phi_k(i+1,j,k)-phi_k(i-1,j,k))/(2*dx);

            dphidz(i,j,k)   = (input_phi(i,j+1,k)-input_phi(i,j-1,k))/(2*dz);
            dphi_kdz(i,j,k) = (phi_k(i,j+1,k)-phi_k(i,j-1,k))/(2*dz);
         }
      }
   }

    for(i = 2; i < imx; ++i){
      for(j = 2; j < jmx; ++j){
         for(k = 0; k <= kmx; ++k){
            d2phidr2(i,j,k)      = c2_over_vA2(i,j)*(dphidr(i+1,j,k)-dphidr(i-1,j,k))/(2*dx*Rgrid[i]);
            d2phi_kdr2(i,j,k)    = c2_over_vA2(i,j)*(dphi_kdr(i+1,j,k)-dphi_kdr(i-1,j,k))/(2*dx*Rgrid[i]);

            d2phidz2(i,j,k)      = c2_over_vA2(i,j)*(dphidz(i,j+1,k)-dphidz(i,j-1,k))/(2*dz);
            d2phi_kdz2(i,j,k)    = c2_over_vA2(i,j)*(dphi_kdz(i,j+1,k)-dphi_kdz(i,j-1,k))/(2*dz);
            
            OPPphi(i,j,k)  = (d2phidr2(i,j,k)+d2phidz2(i,j,k));
            OPPphik(i,j,k) = (d2phi_kdr2(i,j,k)+d2phi_kdz2(i,j,k));

            if(i3D == 0){
               r_hand(i,j,k)  = OPPphi(i,j,k) + q[0]*mu0*den2d2(i,j) - e*mu0*xn0e(i,j)*exp((e/t0e(i,j))*input_phi(i,j,k));
               rk_hand(i,j,k) = OPPphik(i,j,k) + q[0]*mu0*den2d2(i,j) - e*mu0*xn0e(i,j)*exp((e/t0e(i,j))*phi_k(i,j,k));
            }else{
               r_hand(i,j,k)  = OPPphi(i,j,k) + q[0]*mu0*den(1,i,j,k) - e*mu0*xn0e(i,j)*exp((e/t0e(i,j))*input_phi(i,j,k));
               rk_hand(i,j,k) = OPPphik(i,j,k) + q[0]*mu0*den(1,i,j,k) - e*mu0*xn0e(i,j)*exp((e/t0e(i,j))*phi_k(i,j,k));
            }

            input_phi(i,j,k) = phi_k(i,j,k) - (input_phi(i,j,k)-phi_k(i,j,k))*r_hand(i,j,k)/(rk_hand(i,j,k)-r_hand(i,j,k));
            // input_phi(i,j,k) = OPPphik(i,j,k)

            if(input_phi(i,j,k) == 0) {
               cout << "working\n";
            }

            if(mask(i,j) < 0.99){
               for(int kk = 0; kk<= kmx; ++kk) {
                  input_phi(i,j,kk) = 0;
               }
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

void density_filter(CArray4D<double> &input_density){
   CArray3D<double> Filter;

   Filter.resize(imx, jmx, kmx+1);

   for(int i = 1; i < imx; ++i) { //inclusive do loops - imx-1 -> < imx. Also better caching performance for c++ if i,j,k. Not sure for fortran
      for(int j = 1; j < jmx; ++j) {
         for(int k = 0; k <= kmx; ++k) {
            Filter(i,j,k) = (1.0/16.0) * (input_density(1,i-1,j-1,k) + 2*input_density(1,i,j-1,k) + input_density(1,i+1,j-1,k) + 
                           2*input_density(1,i-1,j,k) + 4*input_density(1,i,j,k) + 2*input_density(1,i+1,j,k) + 
                           input_density(1,i-1,j+1,k) + 2*input_density(1,i,j+1,k) + input_density(1,i+1,j+1,k));
         }
      }
   }
   for(int i = 0; i <= imx; ++i) {
      for(int j = 0; j <= jmx; ++j) {
         for(int k = 0; k <= kmx; ++k) {
            input_density(1,i,j,k) = Filter(i,j,k);
         }
      }
   }        
}

void poloidal_filter_methods(CArray3D<double> &input_phi) {
   CArray3D<double> Filter;
   //Quick note, this will get weird due to indexing. From what I understand (what I googled really) there
   // is no way to declare a 1-indexed array in c++ without creating your own class. 
   // It would be much easier to just ignore those indeces, so that's what I'll do
   Filter.resize(imx, jmx, kmx+1);

   for(int i = 1; i < imx; ++i) { //inclusive do loops - imx-1 -> < imx. Also better caching performance for c++ if i,j,k. Not sure for fortran
      for(int j = 1; j < jmx; ++j) {
         for(int k = 0; k <= kmx; ++k) {
            Filter(i,j,k) = (1.0/16.0) * (input_phi(i-1,j-1,k) + 2*input_phi(i,j-1,k) + input_phi(i+1,j-1,k) + 
                           2*input_phi(i-1,j,k) + 4*input_phi(i,j,k) + 2*input_phi(i+1,j,k) + 
                           input_phi(i-1,j+1,k) + 2*input_phi(i,j+1,k) + input_phi(i+1,j+1,k));
         }
      }
   }

   for(int i = 0; i <= imx; ++i) {
      for(int j = 0; j <= jmx; ++j) {
         for(int k = 0; k <= kmx; ++k) {
            if(mask(i,j) < 0.99) {
               input_phi(i,j,k) = 0;
            } else {
               // printf("Filter(%d, %d, %d) = %lf\n", i, j, k, Filter(i,j,k));
               input_phi(i,j,k) = Filter(i,j,k);
            }
         }
      }
   }
}

void radial_binomial_filter(CArray3D<double> &input) {
   CArray3D<double> Filter;
   Filter.resize(imx, jmx, kmx+1);

   int i,j,k,rout,zout,rin,zin;

   double br = 0.0;
   double bz = 0.0;
   double dr = 0.0;

   double R_out,Z_out,R_in,Z_in,inner,outer,area;

   dr = 0.36*1.67 / (45*2); //CBC specific dr for filtering. About 3 times dx
   area = dx*dz;

   for(i=2;i<=imx-2;++i){
      for(j=2;j<=jmx-2;++j){
         for(k=0;k<=kmx;++k){

            if(mask(i,j) < 0.99){
               continue;
            }else{
               br = b0x(i,j)/sqrt(pow(b0x(i,j),2)+pow(b0z(i,j),2));
               bz = b0z(i,j)/sqrt(pow(b0x(i,j),2)+pow(b0z(i,j),2));

               //Outer radial interpolation
               R_out = bz*dr + Rgrid[i];
               Z_out = Zgrid[j] - br*dr;
               rout = static_cast<int>((R_out-Rgrid[0])/dx);
               zout = static_cast<int>((Z_out-Zgrid[0])/dz);

               outer = (input(rout,zout,k)*((Rgrid[rout+1]-R_out)*(Zgrid[zout+1]-Z_out)) + 
                        input(rout+1,zout,k)*((R_out-Rgrid[rout])*(Zgrid[zout+1]-Z_out)) +
                        input(rout,zout+1,k)*((Rgrid[rout+1]-R_out)*(Z_out-Zgrid[zout])) + 
                        input(rout+1,zout+1,k)*((R_out-Rgrid[rout])*(Z_out-Zgrid[zout])))/area;

               //Inner radial interpolation
               R_in = Rgrid[i] - bz*dr;
               Z_in = br*dr + Zgrid[j];
               rin = static_cast<int>((R_in-Rgrid[0])/dx);
               zin = static_cast<int>((Z_in-Zgrid[0])/dz);

               inner = (input(rin,zin,k)*((Rgrid[rin+1]-R_in)*(Zgrid[zin+1]-Z_in)) +
                        input(rin+1,zin,k)*((R_in-Rgrid[rin])*(Zgrid[zin+1]-Z_in)) +
                        input(rin,zin+1,k)*((Rgrid[rin+1]-R_in)*(Z_in-Zgrid[zin])) +
                        input(rin+1,zin+1,k)*((R_in-Rgrid[rin])*(Z_in-Zgrid[zin])))/area;

               Filter(i,j,k) = (1.0/4.0) * (inner + outer + 2*input(i,j,k));

               
               // cout << ((Rgrid[rin+1]-R_in)*(Zgrid[zin+1]-Z_in))/area <<' ' << ((R_in-Rgrid[rin])*(Zgrid[zin+1]-Z_in))/area << ' '<<((Rgrid[rin+1]-R_in)*(Z_in-Zgrid[zin]))/area <<' '<< ((R_in-Rgrid[rin])*(Z_in-Zgrid[zin]))/area << endl;
               
            }
         }
      }
   }
   for(i=2;i<=imx-2;++i){
      for(j=2;j<=jmx-2;++j){
         for(k=0;k<=kmx;++k){
            input(i,j,k) = Filter(i,j,k);
         }
      }
   }
}

void low_mode_filter(CArray3D<double> &input_phi) {
   CArray3D<double> mode_holder;
   mode_holder.resize(imx, jmx, kmx+1);

   CArray3D<double> filtered;
   filtered.resize(imx, jmx, kmx+1);

   CArray3D<double> Filter;
   Filter.resize(imx, jmx, kmx+1);

   int rout,zout,rin,zin;

   double br = 0.0;
   double bz = 0.0;
   double dr = 0.0;
   double hyper_alpha = 0.0;
   double r_1 = 0.0;
   double r_2 = 0.0;
   double r_3 = 0.0;
   double r_4 = 0.0;
   double z_4 = 0.0;
   double rz_2 = 0.0;
   double r_2_z_2 = 0.0;

   double R_out,Z_out,R_in,Z_in,inner,outer,area;

   // hyper_alpha = 0.02*pow(dx,4)/dt;
   hyper_alpha = 0.06*pow(dx,4)/dt;

   dr = 0.36*1.67 / (45*2); //CBC specific dr for filtering. About 3 times dx
   area = dx*dz;

   fftw_complex* phi_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ((kmx + 1) / 2 + 1));
   double* f_filtered = (double*) fftw_malloc(sizeof(double) * (kmx + 1));

   fftw_plan plan_forward = fftw_plan_dft_r2c_1d(kmx + 1, &input_phi(0,0,0), phi_hat, FFTW_ESTIMATE);
   fftw_plan plan_backward = fftw_plan_dft_c2r_1d(kmx + 1, phi_hat, f_filtered, FFTW_ESTIMATE);
   // printf("Testing");

   for(int mode=0; mode <= ((kmx+1)/2); ++mode){
      for(int i = 2; i<=imx-2; ++i){
         for(int j = 2; j<=jmx-2; ++j){
            fftw_execute_dft_r2c(plan_forward, &input_phi(i,j,0), phi_hat);

            for (int m = 0; m <= ((kmx + 1) / 2); ++m) {
               phi_hat[m][0] /= (kmx + 1);
               phi_hat[m][1] /= (kmx + 1);
               
               if (m != mode || mode > ((kmx+1)/2)*0.8) {
                  phi_hat[m][0] = 0.0;
                  phi_hat[m][1] = 0.0;
               }
            }

            fftw_execute_dft_c2r(plan_backward, phi_hat, f_filtered);

            for (int k = 0; k <= kmx; ++k) {
               mode_holder(i,j,k) = f_filtered[k];
               // if (mode == 0){ //Zero mode test, not correct
               //    mode_holder(i,j,k) = 0;
               // }
            }
         }
      }

      if(mode < ((kmx+1)/2)*0.1){
         // radial_binomial_filter(mode_holder);
         // hyperdiffusion_filter(mode_holder);
         for(int i = 3; i<=imx-3; ++i){
            for(int j = 3; j<=jmx-3; ++j){
               for(int k = 0; k<=kmx; ++k){

                  z_4 = (mode_holder(i,j+2,k) - 4*mode_holder(i,j+1,k) + 6*mode_holder(i,j,k) - 4*mode_holder(i,j-1,k) + mode_holder(i,j-2,k))/pow(dz,4);

                  r_4 = (mode_holder(i+2,j,k) - 4*mode_holder(i+1,j,k) + 6*mode_holder(i,j,k) - 4*mode_holder(i-1,j,k) + mode_holder(i-2,j,k))/pow(dx,4);
                  r_3 = (mode_holder(i+2,j,k) - 2*mode_holder(i+1,j,k) + 2*mode_holder(i-1,j,k) - mode_holder(i-2,j,k))/(2*pow(dx,3));
                  r_2 = (mode_holder(i+1,j,k) - 2*mode_holder(i,j,k) + mode_holder(i-1,j,k))/pow(dx,2);
                  r_1 = (mode_holder(i+1,j,k) - mode_holder(i-1,j,k))/(2*dx);

                  rz_2 = ((mode_holder(i+1,j+1,k)-2*mode_holder(i+1,j,k)+mode_holder(i+1,j-1,k))-(mode_holder(i-1,j+1,k)-2*mode_holder(i-1,j,k)+mode_holder(i-1,j-1,k)))/(2*dx*pow(dz,2));

                  r_2_z_2 = ((mode_holder(i+1,j+1,k)-2*mode_holder(i+1,j,k)+mode_holder(i+1,j-1,k)) - 2*(mode_holder(i,j+1,k)-2*mode_holder(i,j,k)+mode_holder(i,j-1,k)) + (mode_holder(i-1,j+1,k)-2*mode_holder(i-1,j,k)+mode_holder(i-1,j-1,k)))/(pow(dx,2)*pow(dz,2));

                  hyper_operator(i,j,k) = z_4 + r_4 + (2/Rgrid[i])*r_3 - (1/pow(Rgrid[i],2))*r_2 + (1/pow(Rgrid[i],3))*r_1 + (2/Rgrid[i])*rz_2 + 2*r_2_z_2;
               }
            }
         }
         for(int i = 3; i<=imx-3; ++i){
            for(int j = 3; j<=jmx-3; ++j){
               for(int k = 0; k<=kmx; ++k){
                  // mode_holder(i,j,k)=Filter(i,j,k);
                  mode_holder(i,j,k)=mode_holder(i,j,k) - dt*hyper_operator(i,j,k)*hyper_alpha;
               }
            }
         }
         for(int i = 3; i<=imx-3; ++i){
            for(int j = 3; j<=jmx-3; ++j){
               for(int k = 0; k<=kmx; ++k){
                  br = b0x(i,j)/sqrt(pow(b0x(i,j),2)+pow(b0z(i,j),2));
                  bz = b0z(i,j)/sqrt(pow(b0x(i,j),2)+pow(b0z(i,j),2));

                  //Outer radial interpolation
                  R_out = bz*dr + Rgrid[i];
                  Z_out = Zgrid[j] - br*dr;
                  rout = static_cast<int>((R_out-Rgrid[0])/dx);
                  zout = static_cast<int>((Z_out-Zgrid[0])/dz);

                  outer = (mode_holder(rout,zout,k)*((Rgrid[rout+1]-R_out)*(Zgrid[zout+1]-Z_out)) + 
                           mode_holder(rout+1,zout,k)*((R_out-Rgrid[rout])*(Zgrid[zout+1]-Z_out)) +
                           mode_holder(rout,zout+1,k)*((Rgrid[rout+1]-R_out)*(Z_out-Zgrid[zout])) + 
                           mode_holder(rout+1,zout+1,k)*((R_out-Rgrid[rout])*(Z_out-Zgrid[zout])))/area;

                  //Inner radial interpolation
                  R_in = Rgrid[i] - bz*dr;
                  Z_in = br*dr + Zgrid[j];
                  rin = static_cast<int>((R_in-Rgrid[0])/dx);
                  zin = static_cast<int>((Z_in-Zgrid[0])/dz);

                  inner = (mode_holder(rin,zin,k)*((Rgrid[rin+1]-R_in)*(Zgrid[zin+1]-Z_in)) +
                           mode_holder(rin+1,zin,k)*((R_in-Rgrid[rin])*(Zgrid[zin+1]-Z_in)) +
                           mode_holder(rin,zin+1,k)*((Rgrid[rin+1]-R_in)*(Z_in-Zgrid[zin])) +
                           mode_holder(rin+1,zin+1,k)*((R_in-Rgrid[rin])*(Z_in-Zgrid[zin])))/area;
                  Filter(i,j,k) = (1.0/4.0) * (inner + outer + 2*mode_holder(i,j,k));
               }
            }
         }
         for(int i = 3; i<=imx-3; ++i){
            for(int j = 3; j<=jmx-3; ++j){
               for(int k = 0; k<=kmx; ++k){
                  mode_holder(i,j,k)=Filter(i,j,k);
               }
            }
         }
      }

      for(int i = 2; i<=imx-2; ++i){
         for(int j = 2; j<=jmx-2; ++j){
            for(int k = 0; k<=kmx; ++k){
               filtered(i,j,k) = filtered(i,j,k) + mode_holder(i,j,k);
            }
         }
      }
      mode_holder.Clear();
   }

   for(int i = 2; i<=imx-2; ++i){
      for(int j = 2; j<=jmx-2; ++j){
         for(int k = 0; k<=kmx; ++k){
            if (mask(i,j) < 0.99){
               input_phi(i,j,k) = 0;
            } else{
               input_phi(i,j,k) = filtered(i,j,k);
            }
         }
      }
   }

   fftw_destroy_plan(plan_forward);
   fftw_destroy_plan(plan_backward);

   fftw_free(phi_hat);
   fftw_free(f_filtered);
}

void binomial_filter(CArray3D<double> &input_phi) {
   int kmin, kmax;
   CArray3D<double> Filter;
   Filter.resize(imx,jmx,kmx+1);

   //Binomial Filter
   for(int i = 1; i < imx; ++i) {
      for(int j = 1; j < jmx; ++j) {
         for(int k = 0; k <= kmx; ++k) {
            if(k == 0) {
               kmin = kmx;
            } else {
               kmin = k-1;
            }

            if(k == kmx) {
               kmax = 0;
            } else {
               kmax = k+1;
            }

            Filter(i,j,k) = (1.0 / 64.0) * (input_phi(i-1,j-1,kmin) + 2*input_phi(i,j-1,kmin) + input_phi(i+1,j-1,kmin) + 
                           2*input_phi(i-1,j,kmin) + 4*input_phi(i,j,kmin) + 2*input_phi(i+1,j,kmin) + 
                           input_phi(i-1,j+1,kmin) + 2*input_phi(i,j+1,kmin) + input_phi(i+1,j+1,kmin) + 
                           2*input_phi(i-1,j-1,k) + 4*input_phi(i,j-1,k) + 2*input_phi(i+1,j-1,k) + 
                           4*input_phi(i-1,j,k) + 8*input_phi(i,j,k) + 4*input_phi(i+1,j,k) + 
                           2*input_phi(i-1,j+1,k) + 4*input_phi(i,j+1,k) + 2*input_phi(i+1,j+1,k) + 
                           input_phi(i-1,j-1,kmax) + 2*input_phi(i,j-1,kmax) + input_phi(i+1,j-1,kmax) + 
                           2*input_phi(i-1,j,kmax) + 4*input_phi(i,j,kmax) + 2*input_phi(i+1,j,kmax) + 
                           input_phi(i-1,j+1,kmax) + 2*input_phi(i,j+1,kmax) + input_phi(i+1,j+1,kmax));
         }
      }
   }

   for(int i = 0; i <= imx; ++i) {
      for(int j = 0; j<= jmx; ++j) {
         for(int k = 0; k <= kmx; ++k) {
            if (mask(i,j) < 0.99) {
               input_phi(i,j,k) = 0.0;
            } else {
               input_phi(i,j,k) = Filter(i,j,k);
            }
         }
      }
   }
}

void fourier_modes(CArray3D<double> &input_phi, const int &modes) {
   fftw_complex* phi_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ((kmx + 1) / 2 + 1));
   double* f_filtered = (double*) fftw_malloc(sizeof(double) * (kmx + 1));

   fftw_plan plan_forward = fftw_plan_dft_r2c_1d(kmx + 1, &input_phi(0,0,0), phi_hat, FFTW_ESTIMATE);
   fftw_plan plan_backward = fftw_plan_dft_c2r_1d(kmx + 1, phi_hat, f_filtered, FFTW_ESTIMATE);
   // printf("Testing");
   for (int i = 0; i <= imx; ++i) {
      for (int j = 0; j <= jmx; ++j) {
         // fftw_plan plan_forward = fftw_plan_dft_r2c_1d(kmx + 1, &input_phi(i,j,0), phi_hat, FFTW_ESTIMATE);
         // fftw_execute(plan_forward);
         // fftw_destroy_plan(plan_forward);

         fftw_execute_dft_r2c(plan_forward, &input_phi(i,j,0), phi_hat);

         for (int k = 0; k <= ((kmx + 1) / 2); ++k) {
               phi_hat[k][0] /= (kmx + 1);
               phi_hat[k][1] /= (kmx + 1);
               // if (k != 0){
               //    if (k != modes) {
               //       phi_hat[k][0] = 0.0;
               //       phi_hat[k][1] = 0.0;
               //    }
               // }
               if (k != modes) {
                  phi_hat[k][0] = 0.0;
                  phi_hat[k][1] = 0.0;
               }
         }

         // fftw_plan plan_backward = fftw_plan_dft_c2r_1d(kmx + 1, phi_hat, f_filtered, FFTW_ESTIMATE);
         // fftw_execute(plan_backward);
         // fftw_destroy_plan(plan_backward);

         fftw_execute_dft_c2r(plan_backward, phi_hat, f_filtered);

         for (int k = 0; k <= kmx; ++k) {
               input_phi(i,j,k) = f_filtered[k];
         }
      }
   }

   fftw_destroy_plan(plan_forward);
   fftw_destroy_plan(plan_backward);

   if(myid == 0 && ((timestep%10) == 0)) {
      // write(*,*) phi(:,:,outk);

      ofstream file;
      file.open("testphi_fourier");
      for(int i = 0; i <= imx; ++i) {
         for(int j = 0; j <= jmx; ++j) {
            file << input_phi(i,j,0) << "    ";
         }
         file << "\n";
      } 
      file.close();
   }

   fftw_free(phi_hat);
   fftw_free(f_filtered);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// void flux_fourier_filter(CArray3D<double> &input) {
//    //Local Variables

//    int i = 0, j = 0, k = 0, l = 0, line = 0;
//    int line_marker, temp_length;
//    int gi = 0;
//    double weightinput=0;
//    double input_mag = 0;
//    double output_mag = 0;

//    std::vector<double> temp_flux_array;
//    std::vector<double> temp_psi_array;
//    std::vector<double> temp_theta_array;

//    //REVISED COMPUTATION
//    // line_marker = 0;

//    for(k = 0; k <= kmx; ++k) {
//       weightinput = 0.0;
//       line_marker = 0;
//       for(i=0; i<=imx; ++i){
//          for(j=0; j<=jmx; ++j){
//             fffoutput(i,j) = input(i,j,k);
//             fffcount(i,j) = 1;
//          }
//       }
//       for(gi = 0; gi < 100; ++gi) {
//          temp_length = -1;
//          temp_flux_array.clear();
//          for(line = line_marker; line<= num_lines; ++line) {
//             if(gindex[line] == gi){
//                weightinput = (weight00[line]*input(iarray[line], jarray[line],k) + 
//                            weight10[line]*input(iarray[line]+1, jarray[line],k) + 
//                            weight01[line]*input(iarray[line], jarray[line]+1,k) + 
//                            weight11[line]*input(iarray[line]+1, jarray[line]+1,k));
//                temp_flux_array.push_back(weightinput);
//                temp_length = temp_length + 1;

//                temp_theta_array.push_back(atan2((Zgrid[jarray[line]])/(Rgrid[iarray[line]]-Rgrid[imx/2])));
//             } else{
//                break;
//             }
//          }

//          if(temp_length < 0){
//             continue;
//          } else{
//             temp_psi_array.push_back(psitab[line_marker]);
//          }

//          ///////// Fourier Filter////////////
//          fftw_complex* flux_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ((temp_length + 1) / 2 +1));
//          double* f_filtered = (double*) fftw_malloc(sizeof(double) * (temp_length + 1));

//          fftw_plan plan_forward = fftw_plan_dft_r2c_1d(temp_length + 1, temp_flux_array.data(), flux_hat, FFTW_ESTIMATE);
//          fftw_plan plan_backward = fftw_plan_dft_c2r_1d(temp_length + 1, flux_hat, f_filtered, FFTW_ESTIMATE);

//          fftw_execute_dft_r2c(plan_forward, &temp_flux_array[0], flux_hat);

//          for (l = 0; l <= ((temp_length + 1)/2); ++l) {
//             flux_hat[l][0] /= (temp_length + 1);
//             flux_hat[l][1] /= (temp_length + 1);
//             if (l > 0.8*((temp_length + 1)/2)){

//                flux_hat[l][0] = 0.0;
//                flux_hat[l][1] = 0.0;
//             }
//          }

//          fftw_execute_dft_c2r(plan_backward, flux_hat, f_filtered);

//          for (l = 0; l <= temp_length; ++l) {
//             temp_flux_array[l] = f_filtered[l];
//          }

//          fftw_destroy_plan(plan_forward);
//          fftw_destroy_plan(plan_backward);

//          fftw_free(flux_hat);
//          fftw_free(f_filtered);
//          ////////////////////////////////////
//          // For interpolation, construct psi,theta grid of contour lines.

//          line_marker = line;
//       }

//       for(i=0; i<=imx; ++i){
//          for(j=0; j=jmx; ++j){
//             psi_lower = temp_psi_array[argmin(psi_p(i,j)-temp_psi_array)];
//             psi_upper = temp_psi_array[argmin(psi_p(i,j)-temp_psi_array)+1];
            
//             theta_lower = argmin(atan2((Zgrid[jarray[line]])/(Rgrid[iarray[line]]-Rgrid[imx/2]))-temp_theta_array);
//             theta_upper = argmin(atan2((Zgrid[jarray[line]])/(Rgrid[iarray[line]]-Rgrid[imx/2]))-temp_theta_array) + 1;


//          }
//       }

//       for(i=0; i<=imx; ++i){
//          for(j=0; j<=jmx; ++j){
//             if(mask(i,j) < 0.99){
//                input(i,j,k) = 0.0;
//             }else{
//                input(i,j,k) = fffoutput(i,j)/fffcount(i,j);
//             }
//          }
//       }
//    }
// }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperdiffusion_filter(CArray3D<double> &input) {

   int i = 0;
   int j = 0;
   int k = 0;

   double r_1 = 0;
   double r_2 = 0;
   double r_3 = 0;
   double r_4 = 0;

   double z_4 = 0;

   double rz_2 = 0;
   double r_2_z_2 = 0;

   double hyper_alpha = 0;

   hyper_alpha = 0.02*pow(dx,4)/dt;

   for(int i = 2; i <= imx-2; ++i) {
      for(int j = 2; j <= jmx-2; ++j) {
         for(int k =0; k <= kmx; ++k){
            z_4 = (input(i,j+2,k) - 4*input(i,j+1,k) + 6*input(i,j,k) - 4*input(i,j-1,k) + input(i,j-2,k))/pow(dz,4);

            r_4 = (input(i+2,j,k) - 4*input(i+1,j,k) + 6*input(i,j,k) - 4*input(i-1,j,k) + input(i-2,j,k))/pow(dx,4);
            r_3 = (input(i+2,j,k) - 2*input(i+1,j,k) + 2*input(i-1,j,k) - input(i-2,j,k))/(2*pow(dx,3));
            r_2 = (input(i+1,j,k) - 2*input(i,j,k) + input(i-1,j,k))/pow(dx,2);
            r_1 = (input(i+1,j,k) - input(i-1,j,k))/(2*dx);

            rz_2 = ((input(i+1,j+1,k)-2*input(i+1,j,k)+input(i+1,j-1,k))-(input(i-1,j+1,k)-2*input(i-1,j,k)+input(i-1,j-1,k)))/(2*dx*pow(dz,2));

            r_2_z_2 = ((input(i+1,j+1,k)-2*input(i+1,j,k)+input(i+1,j-1,k)) - 2*(input(i,j+1,k)-2*input(i,j,k)+input(i,j-1,k)) + (input(i-1,j+1,k)-2*input(i-1,j,k)+input(i-1,j-1,k)))/(pow(dx,2)*pow(dz,2));

            hyper_operator(i,j,k) = z_4 + r_4 + (2/Rgrid[i])*r_3 - (1/pow(Rgrid[i],2))*r_2 + (1/pow(Rgrid[i],3))*r_1 + (2/Rgrid[i])*rz_2 + 2*r_2_z_2;
         }   
      }
   }

   for(int i = 0; i <= imx; ++i) {
      for(int j = 0; j <= jmx; ++j) {
         for(int k = 0; k <= kmx; ++k) {
            if (mask(i,j) < 0.99) {
               input(i,j,k) = 0.0;
            } else {
               input(i,j,k) = input(i,j,k) - dt*hyper_operator(i,j,k)*hyper_alpha;
               // input(i,j,k) = 1;
            }
         }   
      }
   }
}

inline void prepareDeviceData() {
   //ionpush data
   curlb.todev();
   ex.todev();
   ez.todev();
   ezeta.todev();
   dbdx.todev();
   dbdz.todev();
   b0.todev();
   b0x.todev();
   b0z.todev();
   b0zeta.todev();
   captix.todev();
   captiz.todev();
   capnix.todev();
   capniz.todev();
   xn0i.todev();
   delbx.todev();
   delbz.todev();
   t0i.todev();
	dpsi_dr.todev();
	dpsi_dz.todev();

   //integ data
   den.todev();
   upar.todev();

   //particle data arrays
   #pragma acc enter data copyin(mu[0:mmx], u2[0:mmx], u3[0:mmx], x2[0:mmx] ,x3[0:mmx], z2[0:mmx], z3[0:mmx], zeta2[0:mmx], zeta3[0:mmx],  w2[0:mmx], w3[0:mmx], gw[0:mmx])
}

inline void freeDeviceData() {
   //ionpush data
   curlb.fromdev();
   ex.fromdev();
   ez.fromdev();
   ezeta.fromdev();
   dbdx.fromdev();
   dbdz.fromdev();
   b0.fromdev();
   b0x.fromdev();
   b0z.fromdev();
   b0zeta.fromdev();
   captix.fromdev();
   captiz.fromdev();
   capnix.fromdev();
   capniz.fromdev();
   xn0i.fromdev();
   delbx.fromdev();
   delbz.fromdev();
   t0i.fromdev();
	dpsi_dr.fromdev();
	dpsi_dz.fromdev();

   //integ data
   den.fromdev();
   upar.fromdev();

   //particle data arrays
      #pragma acc exit data delete(mu[0:mmx], u2[0:mmx], u3[0:mmx], x2[0:mmx] ,x3[0:mmx], z2[0:mmx], z3[0:mmx], zeta2[0:mmx], zeta3[0:mmx],  w2[0:mmx], w3[0:mmx], gw[0:mmx])

}
