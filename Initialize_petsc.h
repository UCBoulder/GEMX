       one = 1
       three = 3



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
