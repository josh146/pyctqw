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
    subroutine hamiltonian_1p(A,d,amp,nd,n)
        PetscInt, intent(in)    :: nd, n
        real(8), intent(in)     :: amp(nd)
        integer, intent(in)     :: d(nd)
        
        Mat, intent(out)        :: A
        
        ! local variables
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
    
    end subroutine hamiltonian_1p
    
    ! calculate y=e^(-iHt).v using SLEPc    
    subroutine expm(A,t,v,y)
        Vec, intent(in)         :: v
        Mat, intent(in)         :: A
        PetscScalar, intent(in) :: t
        
        Vec, intent(out)        :: y
        
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
    subroutine min_max_eigs(A,rank,Eval,Eval_error,which, &
                & eig_solver,worktype,worktypeInt,tolIn,max_it,verbose,error)
    
        PetscMPIInt, intent(in)      :: rank
        PetscInt, intent(in)         :: worktypeInt, max_it
        PetscReal, intent(in)        :: tolIn
        character(len=3), intent(in) :: which, worktype
        character(len=*),intent(in)  :: eig_solver
        Mat, intent(in)              :: A
        PetscBool, intent(in)        :: verbose
        
        PetscScalar, intent(out)     :: Eval
        PetscReal, intent(out)       :: Eval_error
        PetscInt, intent(out)        :: error
        
        ! local variables
        PetscErrorCode    :: ierr
        PetscInt          :: i, j, its, nev, maxit, nconv
        PetscScalar       :: kr, ki
        PetscReal         :: tol
        EPS               :: eps
        EPSType           :: tname
        character(len=12) :: arg
    
        call EPSCreate(PETSC_COMM_WORLD,eps,ierr)
        call EPSSetOperators(eps,A,PETSC_NULL_OBJECT,ierr)
        
        ! set the problem type as EPS_HEP: Hermitian matrix
        call EPSSetProblemType(eps,EPS_HEP,ierr)
        
        ! set tolerance and max number of iterations
        if (tolIn .ne. 0) then
            if (max_it .ne. 0) then
                call EPSSetTolerances(eps,tolIn,max_it,ierr)
            else
                call EPSSetTolerances(eps,tolIn,PETSC_NULL_OBJECT,ierr)
            endif
        else
            if (max_it .ne. 0) call EPSSetTolerances(eps,PETSC_NULL_OBJECT,max_it,ierr)
        endif
        
        ! set the workspace options
        select case(worktype)
            case('mpd');   call EPSSetDimensions(eps,1,PETSC_NULL_OBJECT,worktypeInt,ierr)
            case('ncv');   call EPSSetDimensions(eps,1,worktypeInt,PETSC_NULL_OBJECT,ierr)
            case default;  continue ! SLEPC default workspace
        end select
        
        ! set eigensolver to use
        arg = trim(adjustl(eig_solver))
        select case(arg)
            case('arnoldi');     call EPSSetType(eps,EPSARNOLDI,ierr)
            case('lanczos');     call EPSSetType(eps,EPSLANCZOS,ierr)
            case('krylovschur'); call EPSSetType(eps,EPSKRYLOVSCHUR,ierr)
            case('gd');          call EPSSetType(eps,EPSGD,ierr)
            case('jd');          call EPSSetType(eps,EPSJD,ierr)
            case('lapack');      call EPSSetType(eps,EPSLAPACK,ierr)
            case('arpack');      call EPSSetType(eps,EPSARPACK,ierr)
            case default;        continue ! SLEPC default is Krylov-Schur
        end select        

	! determine whether min or max Re(lambda) is determined
        select case(which)
            case('max');
                if (rank==0) write(*,*)'Calculating Emax...'
                call EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL,ierr)
            case default
                if (rank==0) write(*,*)'Calculating Emin...'
                call EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL,ierr)
        end select
        
        ! get command line arguments
        call EPSSetFromOptions(eps,ierr)
        ! solve the system
        call EPSSolve(eps,ierr)
        
        ! verbose info
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
            if (rank == 0) then
                write(*,*)'Number of converged eigenpairs: ',nconv
                write(*,*)''
            endif
        endif
        
        ! determine if convergence occurs
        call EPSGetConverged(eps,nconv,ierr)
        if (nconv==0) then
            error = 1
            if (rank==0) write(*,*)'No converged eigenpairs. Reduce tolerance, increase &
        				& number of interations, or increase work vector size'
        else
            call EPSGetEigenpair(eps,0,kr,ki,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
            call EPSComputeRelativeError(eps,0,Eval_error,ierr)
            Eval = PetscRealPart(kr)
            error = 0
        endif
        
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
    call hamiltonian_1p(A,[0],[0.d0],1,n)
    !call MatView(A,PETSC_VIEWER_STDOUT_SELF)
    
    ! Eigenvalue solver
    call min_max_eigs(A,rank,Emin,Emin_error,'min','krylovschur','mpd',50,0.d0,0,.true.,ierr)
    if (rank==0 .and. ierr==0) write(*,*)'    ',PetscRealPart(Emin),'+-', Emin_error
    call min_max_eigs(A,rank,Emax,Emax_error,'max','krylovschur','mpd',50,0.d0,0,.true.,ierr)
    if (rank==0 .and. ierr==0) write(*,*)'    ',PetscRealPart(Emax),'+-',Emax_error
    
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
    !call VecView(psi,PETSC_VIEWER_STDOUT_WORLD,ierr)

    ! destroy matrix/SLEPc
    call MatDestroy(A,ierr)
    call SlepcFinalize(ierr)

end program slepc
