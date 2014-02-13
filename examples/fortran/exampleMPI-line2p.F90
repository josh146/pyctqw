program exampleMPI
    
    use ctqwMPI
    implicit none
#include <finclude/petsc.h>
    
    ! declare variables
    PetscLogStage  :: stage
    PetscLogDouble :: te10, te21,te11, te20, ts0, ts1, tc1, tc0, th10, th11
    PetscMPIInt    :: rank
    PetscErrorCode :: ierr
    PetscBool      :: flag, exist
    PetscInt       :: i, j, its, n, d(2), numproc
    PetscScalar    :: Emin, Emax, t, init_state(2,3), amp(2), interaction
    PetscReal      :: Emin_error, Emax_error
    Mat            :: H, H2
    Vec            :: psi0, psi, psix, psiy

    PetscInt,allocatable  :: adjarray(:,:)
    character(len=50)     :: filename, tmpStr
    
    ! initialize SLEPc and PETSc
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
    
    ! get command line arguments
    call PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-n",n,flag,ierr)
    CHKERRQ(ierr)
    if (flag .eqv. .false.) n = 100

    call PetscOptionsGetReal(PETSC_NULL_CHARACTER,"-t",t,flag,ierr)
    CHKERRQ(ierr)
    if (flag .eqv. .false.) then
        t = 5.d0+0.d0*PETSC_i
    else
        t = PetscRealPart(t)+0.d0*PETSC_i
    endif

    call PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-numproc",numproc,flag,ierr)
    
    ! create the Hamiltonian
    d = [0,1]
    amp = [2.0,1.5]
    interaction = 1.0

    call PetscTime(th10,ierr)
    call PetscLogStageRegister('Hamiltonian',stage,ierr)
    call PetscLogStagePush(stage,ierr)
    call MatCreate(PETSC_COMM_WORLD,H,ierr)
    call hamiltonian_p2_line(H,d,amp,interaction,size(d),n)
    call PetscBarrier(H,ierr)
    call PetscLogStagePop(ierr)
    call PetscTime(th11,ierr)
    !call MatView(H,PETSC_VIEWER_STDOUT_WORLD,ierr)
    
    ! Eigenvalue solver
    call PetscTime(te10,ierr)
    call PetscLogStageRegister('Emax',stage,ierr)
    call PetscLogStagePush(stage,ierr)
    call min_max_eigs(H,rank,Emax,Emax_error,'max','krylovschur','null',35,1.d-2,0,.false.,ierr)
    call PetscLogStagePop(ierr)
    !if (rank==0 .and. ierr==0) write(*,*)'    ',PetscRealPart(Emax),'+-',Emax_error
    call PetscTime(te11,ierr)
   
    ! create vectors
    call VecCreate(PETSC_COMM_WORLD,psi0,ierr)
    call VecSetSizes(psi0,PETSC_DECIDE,n**2,ierr)
    call VecSetFromOptions(psi0,ierr)
   
    call VecDuplicate(psi0,psi,ierr)
   
    ! create initial state
    call PetscLogStageRegister('InitState',stage,ierr)
    call PetscLogStagePush(stage,ierr)
    init_state(1,:) = [real(0),real(0),1.0/sqrt(2.0)]
    init_state(2,:) = [real(1),real(1),1.0/sqrt(2.0)]
    call p2_init(psi0,init_state,2,n)
    !call VecView(psi0,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call PetscBarrier(psi0,ierr)
    call PetscLogStagePop(ierr)
  
    ! matrix exponential
    call PetscTime(ts0,ierr)
    call PetscLogStageRegister('SLEPc expm',stage,ierr)
    call PetscLogStagePush(stage,ierr)
    call qw_krylov(H,t,psi0,psi)
    call VecView(psi,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call PetscBarrier(psi,ierr)
    call PetscLogStagePop(ierr)
    call PetscTime(ts1,ierr)
  
    ! QW chebyshev
    call PetscTime(tc0,ierr)
    call PetscLogStageRegister('Chebyshev',stage,ierr)
    call PetscLogStagePush(stage,ierr)
    call qw_cheby(psi0,psi,t,H,0.*PETSc_i,Emax,rank,n)
    call VecView(psi,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call PetscBarrier(psi,ierr)
    call PetscLogStagePop(ierr)
    call PetscTime(tc1,ierr)

    ! destroy matrix/SLEPc
    call MatDestroy(H,ierr)
    call VecDestroy(psi,ierr)
    call VecDestroy(psi0,ierr)
    call PetscFinalize(ierr)
    
    if (rank==0) then
        write(tmpStr,*)numproc
        write(filename,*)int(PetscRealPart(t))
        filename = "out/line-2p-int-t" // trim(adjustl(filename)) &
                    & // "-nodes" // trim(adjustl(tmpStr)) // "_intel.txt"
        inquire(file=filename, exist=exist)
        if (exist) then
            open(12, file=filename, status="old", position="append", action="write")
        else
            open(12, file=filename, status="new", action="write")
        end if
        !write(12,15) n,t,th11-th10,te11-te10,ts1-ts0,tc1-tc0
        write(12,15) n,t,th11-th10,te11-te10,ts1-ts0,tc1-tc0
        15 format(i5,f8.2,5e20.10)
        close(12)
    endif

end program exampleMPI
