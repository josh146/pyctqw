!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ctqwMPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~ Todo ~~~~~~~~~~~~~~~~~~~~~~
!  1) Kronecker product
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
program main
    
    use ctqwMPI
    implicit none
#include <finclude/petsc.h>
    
    ! declare variables
    PetscLogStage  :: stage
    PetscLogDouble :: te10, te21,te11, te20, ts0, ts1, tc1, tc0
    PetscMPIInt    :: rank
    PetscErrorCode :: ierr
    PetscBool      :: flag
    PetscInt       :: i, j, its, n, d(2)
    PetscScalar    :: Emin, Emax, t, init_state(2,2), amp(2)
    PetscReal      :: Emin_error, Emax_error
    Mat            :: H, H2
    Vec            :: psi0, psi, psix, psiy
    
    ! initialize SLEPc and PETSc
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
    
    ! get command line arguments
    call PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-n",n,flag,ierr)
    CHKERRQ(ierr)
    if (flag .eqv. .false.) n = 100
    call PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-t",t,flag,ierr)
    CHKERRQ(ierr)
    if (flag .eqv. .false.) t = 10.0
    
    !write(*,*)'hi'
    
    ! create the Hamiltonian
    d = [3,4]
    amp = [2.0,1.5]
    call PetscLogStageRegister('Hamiltonian',stage,ierr)
    call PetscLogStagePush(stage,ierr)
    call MatCreate(PETSC_COMM_WORLD,H,ierr)
    call hamiltonian_1p_line(H,d,amp,size(d),n)
    call PetscBarrier(H,ierr)
    call PetscLogStagePop(ierr)
!    call MatView(H,PETSC_VIEWER_STDOUT_WORLD,ierr)
    
    ! Eigenvalue solver
    call PetscTime(te10,ierr)
    call PetscLogStageRegister('Emax',stage,ierr)
    call PetscLogStagePush(stage,ierr)
    call min_max_eigs(H,rank,Emax,Emax_error,'max','krylovschur','null',35,1.d-2,0,.false.,ierr)
    call PetscLogStagePop(ierr)
    !if (rank==0 .and. ierr==0) write(*,*)'    ',PetscRealPart(Emax),'+-',Emax_error
    call PetscTime(te11,ierr)
    
    call PetscTime(te20,ierr)
    call PetscLogStageRegister('Emin',stage,ierr)
    call PetscLogStagePush(stage,ierr)
    call min_max_eigs(H,rank,Emin,Emin_error,'min','krylovschur','null',35,1.d-2,0,.false.,ierr)
    call PetscLogStagePop(ierr)
    call PetscBarrier(H,ierr)
   ! if (rank==0 .and. ierr==0) write(*,*)'    ',PetscRealPart(Emin),'+-', Emin_error
    call PetscTime(te21,ierr)
    
    ! create vectors
    call VecCreate(PETSC_COMM_WORLD,psi0,ierr)
    call VecSetSizes(psi0,PETSC_DECIDE,n,ierr)
    !call VecSetBlockSize(psi0,n,ierr)
    call VecSetFromOptions(psi0,ierr)
    
    call VecDuplicate(psi0,psi,ierr)
    
    ! create initial state
    call PetscLogStageRegister('InitState',stage,ierr)
    call PetscLogStagePush(stage,ierr)
    init_state(1,:) = [0.,1.0/sqrt(2.0)]
    init_state(2,:) = [1.,1.0/sqrt(2.0)]
    call p1_init(psi0,init_state,2,n)
!    call VecView(psi0,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call PetscBarrier(psi0,ierr)
    call PetscLogStagePop(ierr)
    
    ! matrix exponential
    call PetscTime(ts0,ierr)
    call PetscLogStageRegister('SLEPc expm',stage,ierr)
    call PetscLogStagePush(stage,ierr)
    call expm(H,t,psi0,psi)
    !call VecView(psi,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call PetscBarrier(psi,ierr)
    call PetscLogStagePop(ierr)
    call PetscTime(ts1,ierr)
    
    ! QW chebyshev
    call PetscTime(tc0,ierr)
    call PetscLogStageRegister('Chebyshev',stage,ierr)
    call PetscLogStagePush(stage,ierr)
    call qw_cheby(psi0,psi,t,H,0.*PETSc_i,Emax,rank,n)
    !call VecView(psi,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call PetscBarrier(psi,ierr)
    call PetscLogStagePop(ierr)
    call PetscTime(tc1,ierr)

!    ! get marginal prob
!    call VecCreate(PETSC_COMM_WORLD,psix,ierr)
!    call VecSetSizes(psix,PETSC_DECIDE,n,ierr)
!    call VecSetFromOptions(psix,ierr)
!    call VecDuplicate(psix,psiy,ierr)
!    
!    call PetscLogStageRegister('ProbX',stage,ierr)
!    call PetscLogStagePush(stage,ierr)
!    call marginal(psi,psix,'x',n)
!    call VecView(psix,PETSC_VIEWER_STDOUT_WORLD,ierr)
!    call PetscBarrier(psiX,ierr)
!    call PetscLogStagePop(ierr)
!    
!    call PetscLogStageRegister('ProbY',stage,ierr)
!    call PetscLogStagePush(stage,ierr)
!    call marginal(psi,psiy,'y',n)
!    !call VecView(psiy,PETSC_VIEWER_STDOUT_WORLD,ierr)
!    call PetscBarrier(psiY,ierr)
!    call PetscLogStagePop(ierr)    

    ! destroy matrix/SLEPc
    call MatDestroy(H,ierr)
    call VecDestroy(psi,ierr)
    call VecDestroy(psi0,ierr)
    !call VecDestroy(psix,ierr)
    !call VecDestroy(psiy,ierr)
    
    call PetscFinalize(ierr)
    if (rank==0) then
        open(15,file='petsc20.txt',access='APPEND')
        write(15,'(i10,4E18.10E3)')n,te11-te10,te21-te20,tc1-tc0,ts1-ts0    
        close(15)
    endif

end program main
