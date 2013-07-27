program slepc
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/slepcsys.h>
#include <finclude/slepceps.h>

	use eigs
    
    PetscMPIInt    :: rank
    PetscErrorCode :: ierr
    PetscBool      :: flag
    PetscInt       :: i, j, its, n, Istart, Iend, col(3), &
				& nev, maxit, nconv, one, two, three
    PetscScalar    :: value(3), kr, ki
    PetscReal      :: norm, tol, error
    Mat            :: A
    Vec            :: Rhs, Sol
    EPS            :: eps
    EPSType        :: tname
    
    call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
    
    call PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-n",n,flag,ierr)
    CHKERRQ(ierr)
    if (flag .eqv. .false.) n = 10
    
    call MatCreate(PETSC_COMM_WORLD,A,ierr)
    call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
    call MatSetFromOptions(A,ierr)
    call MatSetUp(A,ierr)
    !call MatSetType(A,MATMPIAIJ,ierr)
    
    call MatGetOwnershipRange(A,Istart,Iend,ierr)
    
    if (Istart == 0) then
        col(1:2) = [0,1]
        value(1:2) = [2.0,-1.0]
        call MatSetValues(A,1,0,2,col,value,INSERT_VALUES,ierr)
        Istart = Istart + 1
    end if
    
    if (Iend == n) then
        col(1:2) = [n-2,n-1]
        value(1:2) = [-1.0,2.0]
        call MatSetValues(A,1,n-1,2,col,value,INSERT_VALUES,ierr)
        Iend = Iend - 1
    endif
    
    value(1:3) = [-1.0,2.0,-1.0]
    do i=Istart, Iend-1
        col = [i-1,i,i+1]
        call MatSetValues(A,1,i,3,col,value,INSERT_VALUES,ierr)
    enddo
    
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    
   ! call MatView(A,PETSC_VIEWER_STDOUT_SELF)
    
    call EPSCreate(PETSC_COMM_WORLD,eps,ierr)
    call EPSSetOperators(eps,A,PETSC_NULL_OBJECT,ierr)
    call EPSSetProblemType(eps,EPS_HEP,ierr)
    call EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL)
    call EPSSetFromOptions(eps,ierr)
    
    call EPSSolve(eps,ierr)
    call EPSGetIterationNumber(eps,its,ierr)
    if (rank == 0) write(*,*)'Number of iterations of the method: ',its
    
    call EPSGetType(eps,tname,ierr)
    if (rank == 0) write(*,*)'Soln method: ',tname
    
    call EPSGetDimensions(eps,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)
    if (rank == 0) write(*,*)'Number of requested eigenvalues: ',nev
    
    call EPSGetTolerances(eps,tol,maxit,ierr)
    if (rank == 0) then
    	write(*,*)'Stopping condition tolerance: ',tol
    	write(*,*)'Max number of iterations: ',maxit
    end if
    
    call EPSGetConverged(eps,nconv,ierr)
    if (rank == 0) write(*,*)'Number of converged eigenpairs: ',nconv
    
    if (nconv>0) then
        if (rank==0) then
            write(*,*)
            write(*,*) '         k                Error      '
            write(*,*) ' ----------------- ------------------'
        end if
    endif
    
    do i=0, nconv-1
        call EPSGetEigenpair(eps,i,kr,ki,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
        call EPSComputeRelativeError(eps,i,error,ierr)
        if (rank==0) write(*,'(2E18.10)') PetscRealPart(kr), error
    end do
    
    call EPSDestroy(eps,ierr)

    call MatDestroy(A,ierr)
    call SlepcFinalize(ierr)

end program slepc
