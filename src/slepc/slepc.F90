module eigs
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/slepcsys.h>
#include <finclude/slepceps.h>
#include <finclude/slepcmfn.h>
    
    contains
    
    ! create a sparse hamiltonian matrix of size n
    ! using PETSc's sparse matrix routines
    subroutine hamiltonian1D(A,n)
        PetscInt, intent(in)    :: n
        Mat, intent(out)        :: A
        
        PetscErrorCode :: ierr
        PetscBool      :: flag
        PetscInt       :: i, j, its, Istart, Iend, col(3)
        PetscScalar    :: value(3)
        
        call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
        call MatSetFromOptions(A,ierr)
        call MatSetUp(A,ierr)
        
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
    
    end subroutine hamiltonian1D
    
    ! calculate the matrix exponential using SLEPc    
    subroutine expm(A,t,v,y)
        Vec, intent(in)    :: v
        Vec, intent(out)   :: y
        Mat, intent(in)    :: A
        PetscScalar, intent(in) :: t
        
        PetscErrorCode :: ierr
        PetscScalar    :: alpha
        MFN            :: mfn
        
        call MFNCreate(PETSC_COMM_WORLD, mfn, ierr)
        call MFNSetOperator(mfn,A,ierr)
        call MFNSetFunction(mfn,SLEPC_FUNCTION_EXP,ierr)
        call MFNSetFromOptions(mfn,ierr)
        
        alpha = -t*PETSC_i
        call MatScale(A,alpha,ierr)
        !call MFNSetScaleFactor(mfn, alpha,ierr)
        call MFNSolve(mfn,v,y,ierr)
        
        call MFNDestroy(mfn,ierr)
    
    end subroutine expm
    
    ! calculate the min or max eiganvalue using SLEPc 
    subroutine min_max_eigs(A,rank,Eval,Eval_error,which,verbose)
    
        PetscMPIInt, intent(in)    :: rank
        character(len=3), intent(in) :: which
        Mat, intent(in)            :: A
        PetscScalar, intent(out)   :: Eval
        PetscReal, intent(out)     :: Eval_error
        PetscBool, intent(in)      :: verbose
        
        ! local variables
        PetscErrorCode :: ierr
        PetscInt       :: i, j, its, nev, maxit, nconv
        PetscScalar    :: kr, ki
        PetscReal      :: tol, error
        EPS            :: eps
        EPSType        :: tname
    
        call EPSCreate(PETSC_COMM_WORLD,eps,ierr)
        call EPSSetOperators(eps,A,PETSC_NULL_OBJECT,ierr)
        call EPSSetProblemType(eps,EPS_HEP,ierr)
        
        if (which=='max') then
            write(*,*)'Calculating Emax...'
            call EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL)
        else
            write(*,*)'Calculating Emin...'
            call EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL)
        endif
        
        call EPSSetFromOptions(eps,ierr)
        
        call EPSSolve(eps,ierr)
        
        if (verbose) then
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
        endif
        
        call EPSGetEigenpair(eps,0,kr,ki,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
        call EPSComputeRelativeError(eps,0,Eval_error,ierr)
        
        Eval = PetscRealPart(kr)
        
        call EPSDestroy(eps,ierr)
        
    end subroutine min_max_eigs

end module eigs

program slepc

    use eigs
    
    PetscMPIInt    :: rank
    PetscErrorCode :: ierr
    PetscBool      :: flag
    PetscInt       :: i, j, its, n
    PetscScalar    :: Emin, Emax, t, value(3)
    PetscReal      :: Emin_error, Emax_error
    Mat            :: A
    Vec            :: psi0, psi
    
    ! initialize SLEPc and PETSc
    call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
    
    ! get command line arguments
    call PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-n",n,flag,ierr)
    CHKERRQ(ierr)
    if (flag .eqv. .false.) n = 10
    
    ! matrix creater
    call MatCreate(PETSC_COMM_WORLD,A,ierr)
    call hamiltonian1D(A,n)
    !call MatView(A,PETSC_VIEWER_STDOUT_SELF)
    
    ! Eigenvalue solver
    call min_max_eigs(A,rank,Emin,Emin_error,'min',.false.)
    write(*,*)'min:',PetscRealPart(Emin), Emin_error
    call min_max_eigs(A,rank,Emax,Emax_error,'max',.false.)
    write(*,*)'max:',PetscRealPart(Emax), Emax_error
    
    ! create vectors
    call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,n,psi0,ierr)
    call VecSetFromOptions(psi0,ierr)
    call VecDuplicate(psi0,psi,ierr)
    
    value = [1.0,2.0,-0.5]
    call VecSetValues(psi0,3,[0,1,2],value,INSERT_VALUES,ierr)
    call VecAssemblyBegin(psi0,ierr)
    call VecAssemblyEnd(psi0,ierr)
    
    ! matrix exponential
    t = 5.0
    call expm(A,t,psi0,psi)
    call VecView(psi,PETSC_VIEWER_STDOUT_WORLD,ierr)

    ! destroy matrix/SLEPc
    call MatDestroy(A,ierr)
    call SlepcFinalize(ierr)

end program slepc
