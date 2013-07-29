module eigs

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscvec.h90>

#include <finclude/slepcsys.h>
#include <finclude/slepceps.h>
#include <finclude/slepcmfn.h>
    
    contains
    
    ! create a sparse hamiltonian matrix of size n
    ! using PETSc's sparse matrix routines
    subroutine hamiltonian_1p(A,d,amp,nd,n)
        PetscInt, intent(in)    :: nd, n, d(nd)
        PetscScalar, intent(in) :: amp(nd)
        
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
            
            do j = 1, nd
                if (int(d(j)+n/2.0) == 1) then
                    value(1) = 2.0+amp(j)
                    exit
                endif
            end do
            
            call MatSetValues(A,1,0,2,col,value,INSERT_VALUES,ierr)
            Istart = Istart + 1
        end if
        
        if (Iend == n) then
            col(1:2) = [n-2,n-1]
            value(1:2) = [-1.0,2.0]
            
            do j = 1, nd
                if (int(d(j)+n/2.0) == n) then
                    value(2) = 2.0+amp(j)
                    exit
                endif
            end do
            
            call MatSetValues(A,1,n-1,2,col,value,INSERT_VALUES,ierr)
            Iend = Iend - 1
        endif
        
        value(1:3) = [-1.0,2.0,-1.0]
        do i=Istart, Iend-1
            col = [i-1,i,i+1]
            
            value = [-1.0,2.0,-1.0]
            do j = 1, nd
                if (int(d(j)+n/2.0) == i+1) then
                    value(2) = 2.0+amp(j)
                    exit
                endif
            end do
            
            call MatSetValues(A,1,i,3,col,value,INSERT_VALUES,ierr)
        enddo
        
        call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    
    end subroutine hamiltonian_1p
    
    
    ! create a sparse 2P hamiltonian matrix of size n^2 x n^2
    ! using PETSc's sparse matrix routines
    subroutine hamiltonian_2p(H1,H2,n)
        PetscInt, intent(in)    :: n
        Mat, intent(in)         :: H1
        
        Mat, intent(out)        :: H2
        
        ! local variables
        VecScatter     :: ctx
        PetscErrorCode :: ierr
        PetscBool      :: flag
        PetscInt       :: NN, i, j, its, Istart, Iend, col(3)
        PetscScalar    :: value(3)
        Mat            :: temp
        Vec            :: diag, diagAll
        PetscScalar, pointer :: diagArray(:)
        
        NN = n**2
        
        ! create H2 matrix
        call MatSetSizes(H2,PETSC_DECIDE,PETSC_DECIDE,NN,NN,ierr)
        call MatSetFromOptions(H2,ierr)
        call MatSetUp(H2,ierr)
        
        ! create temp matrix
        !call MatDuplicate(H2,MAT_DO_NOT_COPY_VALUES,temp,ierr)
        call MatCreate(PETSC_COMM_WORLD,temp,ierr)
        call MatSetSizes(temp,PETSC_DECIDE,PETSC_DECIDE,NN,NN,ierr)
        call MatSetFromOptions(temp,ierr)
        call MatSetUp(temp,ierr)
        
        ! create vector of diagonal entries
        call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,n,diag,ierr)
        call VecSetFromOptions(diag,ierr)
        call MatGetDiagonal(H1,diag,ierr)
        
        ! scatter the vector to an array on *all* nodes
        call VecScatterCreateToAll(diag,ctx,diagAll,ierr)
        call VecScatterBegin(ctx,diag,diagAll,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(ctx,diag,diagAll,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterDestroy(ctx)

        call VecGetArrayF90(diagAll,diagArray,ierr)
        
        ! create KronProd(H1,I) matrix and store in temp
        call MatGetOwnershipRange(temp,Istart,Iend,ierr)
        
        if (Istart == 0) then
            col(1:2) = [0,1]
            value(1:2) = [diagArray(1),-1.d0+0*PETSC_i]            
            call MatSetValues(temp,1,0,2,col,value,INSERT_VALUES,ierr)
            Istart = Istart + 1
        end if
        
        if (Iend == NN) then
            col(1:2) = [NN-2,NN-1]
            value(1:2) = [-1.d0+0*PETSC_i,diagArray(n)]            
            call MatSetValues(temp,1,NN-1,2,col,value,INSERT_VALUES,ierr)
            Iend = Iend - 1
        endif
        
        do i=Istart, Iend-1
            if (mod(1,n)+1==1) then
            
                col(1:2) = [i,i+1]
                value(1:2) = [diagArray(mod(i,n)+1),-1.d0+0*PETSC_i]            
                call MatSetValues(temp,1,i,2,col,value,INSERT_VALUES,ierr)
            
            elseif (mod(i,n)+1==n) then
            
                col(1:2) = [i-1,i]
                value(1:2) = [-1.d0+0*PETSC_i,diagArray(mod(i,n)+1)]            
                call MatSetValues(temp,1,i,2,col,value,INSERT_VALUES,ierr)
            
            else
                col = [i-1,i,i+1]
                value = [-1.d0+0*PETSC_i,diagArray(mod(i,n)+1),-1.d0+0*PETSC_i]            
                call MatSetValues(temp,1,i,3,col,value,INSERT_VALUES,ierr)
            endif
        enddo
        
        call MatAssemblyBegin(temp,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(temp,MAT_FINAL_ASSEMBLY,ierr)
        
        ! create KronProd(I,H1) matrix and store in H2
        call MatGetOwnershipRange(H2,Istart,Iend,ierr)
        
        do while (Istart<n)
            col(1:2) = [Istart,Istart+n]
            value(1:2) = [diagArray(1),-1.d0+0*PETSC_i]            
            call MatSetValues(H2,1,Istart,2,col,value,INSERT_VALUES,ierr)
            Istart = Istart + 1
        end do
        
        do while (Iend>NN-n)
            col(1:2) = [Iend-2-n,Iend-1]
            value(1:2) = [-1.d0+0*PETSC_i,diagArray(n)]            
            call MatSetValues(H2,1,Iend-1,2,col,value,INSERT_VALUES,ierr)
            Iend = Iend - 1
        end do
        
        do i=Istart, Iend-1
            col = [i-n,i,i+n]
            value = [-1.d0+0*PETSC_i,diagArray(int(real(i)/6.0)+1),-1.d0+0*PETSC_i]
            call MatSetValues(H2,1,i,3,col,value,INSERT_VALUES,ierr)
        enddo
        
        call MatAssemblyBegin(H2,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(H2,MAT_FINAL_ASSEMBLY,ierr)
        
        call VecRestoreArrayF90(diag,diagArray,ierr)
        call VecDestroy(diag,ierr)
        
        ! create 2P hamltonian H2=KronProd(I,H1)+KronProd(H1,I)
        call MatAXPY(H2,1.0+0*PETSC_i,temp,DIFFERENT_NONZERO_PATTERN,ierr)
        
        call MatDestroy(temp,ierr)
    
    end subroutine hamiltonian_2p
    
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
        call MatScale(A,1.0/alpha,ierr)
        
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
                call EPSSetTolerances(eps,tolIn,PETSC_DECIDE,ierr)
            endif
        else
            if (max_it .ne. 0) call EPSSetTolerances(eps,PETSC_DECIDE,max_it,ierr)
        endif
        
        ! set the workspace options
        select case(worktype)
            case('mpd');   call EPSSetDimensions(eps,1,PETSC_DECIDE,worktypeInt,ierr)
            case('ncv');   call EPSSetDimensions(eps,1,worktypeInt,PETSC_DECIDE,ierr)
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
    
    subroutine qw_cheby(psi0,psi,dt,H,Emin,Emax,rank,N)
        PetscMPIInt, intent(in) :: rank
        PetscInt, intent(in)    :: N
        PetscScalar, intent(in) :: dt, Emin, Emax
        Vec, intent(in)         :: psi0
        Mat, intent(in)         :: H
        
        Vec, intent(out)        :: psi
        
        ! local variables
        PetscInt      :: m, terms, i, j
        PetscReal     :: alpha
        PetscScalar   :: bessj0, bessj1, bessjn
        Vec, pointer  :: work(:)
        
        alpha = PetscRealPart((Emax-Emin)*dt/2.0)
        
        call VecDuplicateVecsF90(psi0,4,work,ierr)
        
        call VecCopy(psi0,work(1),ierr)
        call MatMult(H,work(1),work(2),ierr)
        call VecAXPBY(work(2), (Emax+Emin)/(Emax-Emin),-2.0/(Emax-Emin), work(1),ierr)
        
        bessj0 = dbesjn(0,alpha)
        bessj1 = dbesjn(1,alpha)
        call VecCopy(psi0,work(4),ierr)
        call VecAXPBY(work(4), 2.0*PETSC_i*bessj1,bessj0, work(2),ierr)
        
        terms = 0
        do while (abs(2.d0*dbesjn(terms,alpha)) > 1.d-18)
            terms = terms + 1
        end do
        
        do m = 2, terms
            call MatMult(H,work(2),work(3),ierr)
            call VecAXPBY(work(3), 2.0*(Emax+Emin)/(Emax-Emin),-4.0/(Emax-Emin), work(2),ierr)
            call VecAXPY(work(3),-1.0+0.*PETSC_i,work(1),ierr)
            
            bessjn = dbesjn(m,alpha)
            call VecAXPY(work(4),2.0*(PETSC_i**m)*bessjn,work(3),ierr)

            call VecCopy(work(2),work(1),ierr)
            call VecCopy(work(3),work(2),ierr)
        end do
        
        call VecScale(work(4),exp(-PETSC_i*(Emax+Emin)*dt/2.0),ierr)
        call VecCopy(work(4),psi,ierr)
        
        call VecDestroyVecsF90(4,work,ierr)
    
    end subroutine qw_cheby

end module eigs

!~~~~~~~~~~~~~~~~~~ Todo ~~~~~~~~~~~~~~~~~~~~~~
!  1) Marginal prob subroutine
!  2) Create intial statespace subroutine
!       -this would require a coord function
!  3) Test with f2py
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

program slepc

    use eigs
    
    PetscViewer    :: matlabV
    PetscMPIInt    :: rank
    PetscErrorCode :: ierr
    PetscBool      :: flag
    PetscInt       :: i, j, its, n, d(4)
    PetscScalar    :: Emin, Emax, t, value(3), amp(4)
    PetscReal      :: Emin_error, Emax_error
    Mat            :: H1, H2
    Vec            :: psi0, psi
    
    ! initialize SLEPc and PETSc
    call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
    
    ! get command line arguments
    call PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-n",n,flag,ierr)
    CHKERRQ(ierr)
    if (flag .eqv. .false.) n = 10
    
    ! 1p matrix creater
    d = [-2,3,0,2]
    amp = [-0.50,5.0,1.0,0.50]
    call MatCreate(PETSC_COMM_WORLD,H1,ierr)
    call hamiltonian_1p(H1,d,amp,size(d),n)
    !call MatView(H1,PETSC_VIEWER_STDOUT_WORLD,ierr)
    
    ! 2p matrix creater
    call MatCreate(PETSC_COMM_WORLD,H2,ierr)
    call hamiltonian_2p(H1,H2,n)
    !call MatView(H2,PETSC_VIEWER_STDOUT_WORLD,ierr)
    
    ! Eigenvalue solver
    call min_max_eigs(H2,rank,Emin,Emin_error,'min','krylovschur','mpd',50,0.d0,0,.false.,ierr)
    if (rank==0 .and. ierr==0) write(*,*)'    ',PetscRealPart(Emin),'+-', Emin_error
    call min_max_eigs(H2,rank,Emax,Emax_error,'max','krylovschur','mpd',50,0.d0,0,.false.,ierr)
    if (rank==0 .and. ierr==0) write(*,*)'    ',PetscRealPart(Emax),'+-',Emax_error
    
    ! create vectors
    call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,n**2,psi0,ierr)
    call VecSetFromOptions(psi0,ierr)
    call VecDuplicate(psi0,psi,ierr)
    
    ! create initial state
    value = [1.0,2.0,-0.5]
    call VecSetValues(psi0,3,[0,1,2],value,INSERT_VALUES,ierr)
    call VecAssemblyBegin(psi0,ierr)
    call VecAssemblyEnd(psi0,ierr)
    
    ! matrix exponential
    t = 5.0
    call expm(H2,t,psi0,psi)
    !call VecView(psi,PETSC_VIEWER_STDOUT_WORLD,ierr)
    
    ! QW chebyshev
    call qw_cheby(psi0,psi,t,H2,Emin,Emax,rank,n**2)
    call VecView(psi,PETSC_VIEWER_STDOUT_WORLD,ierr)

    ! destroy matrix/SLEPc
    call MatDestroy(H1,ierr)
    call MatDestroy(H2,ierr)
    call VecDestroy(psi,ierr)
    call VecDestroy(psi0,ierr)
    
    call SlepcFinalize(ierr)

end program slepc
