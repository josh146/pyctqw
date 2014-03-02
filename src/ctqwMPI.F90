!
!  This file contains the parallelised CTQW propagation and analysis functions
!  ----------------------------------------------------------------------------
!  pyCTQW - Distributed memory CTQW Fortran library and Python module
!  Copyright (C) 2013-2014, Joshua Izaac
!
!  pyCTQW is free software: you can redistribute it and/or modify it under the
!  terms of version 3 of the GNU General Public License as published by
!  the Free Software Foundation.
!
!  pyCTQW is distributed in the hope that it will be useful, but WITHOUT ANY
!  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
!  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
!  more details.
!
!  You should have received a copy of the GNU General Public License
!  along with pyCTQW. If not, see <http://www.gnu.org/licenses/>.
!  ----------------------------------------------------------------------------
!
!
module ctqwMPI

#ifdef __INTEL_COMPILER
    use IFPORT
#endif

    implicit none

    private

    public    :: exportVec, importVec, exportMat, importMat
    public    :: importAdjToH, adjToH

    public    :: identity, kron
    public    :: coord, coord3P

    public    :: marginal1, marginal2, marginal3
    public    :: p1_init, p2_init, p3_init
    public    :: hamiltonian_p1_line, hamiltonian_p2_line, hamiltonian_p3_line
    public    :: min_max_eigs, qw_krylov, qw_cheby

    public    :: partial_trace_array, partial_trace_mat
    public    :: entanglement

    public    :: number_of_edges, getEdgeState
    public    :: getAllEdgeStates, getAllEdgeStates3P
    public    :: GraphISCert
    public    :: refsor

#define PETSC_AVOID_DECLARATIONS
#include <finclude/petsc.h>
#include <finclude/petscvec.h90>
#include <finclude/petscviewer.h90>

#include <finclude/slepcsys.h>
#include <finclude/slepceps.h>
#include <finclude/slepcmfn.h>
#undef PETSC_AVOID_DECLARATIONS

    Integer, parameter   :: kdp = selected_real_kind(15)
    private              :: kdp

    interface refsor
        module procedure d_refsor, r_refsor, i_refsor
    end interface refsor

    contains

    include 'sort.f90'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Vector I/O ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine exportVec(vec,filename,filetype)
        ! export a PETSc vector to a file
        Vec, intent(in)                :: vec ! the PETSc vector to export
        character(len=50), intent(in)  :: filename ! filename *only*, no path, of the output file
        character(len=3), intent(in)   :: filetype ! filetype to export ('txt' or 'bin')
        
        ! local variables
        PetscMPIInt          :: rank
        PetscErrorCode       :: ierr
        PetscViewer          :: binSave
        
        if (filetype == 'bin') then
            call PetscViewerBinaryOpen(PETSC_COMM_WORLD,trim(adjustl(filename)), &
                                         & FILE_MODE_WRITE,binSave,ierr)
            call VecView(vec,binSave,ierr)
            call PetscViewerDestroy(binSave,ierr)
            
        else if (filetype == 'txt') then
            call PetscViewerASCIIOpen(PETSC_COMM_WORLD,trim(adjustl(filename)),binSave,ierr)
            call PetscViewerSetFormat(binSave,PETSC_VIEWER_ASCII_SYMMODU,ierr)
            call VecView(vec,binSave,ierr)
            call PetscViewerDestroy(binSave,ierr)
        endif
    end subroutine exportVec
    
    subroutine importVec(vec,filename,filetype)
        ! import a PETSc vector from a file
        character(len=50), intent(in)  :: filename ! filename *only*, no path, of the input file
        character(len=3), intent(in)   :: filetype ! filetype to export ('txt' or 'bin')
        
        Vec, intent(out)                :: vec ! the PETSc vector to import the file into
        
        ! local variables
        PetscMPIInt          :: rank
        PetscErrorCode       :: ierr
        PetscViewer          :: binLoad
        PetscInt             :: stat, line, i, Istart, Iend
        PetscReal            :: reVec, imVec
        PetscScalar, allocatable :: vecArray(:)
        character(len=100)   :: buffer
        
        call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
        
        if (filetype == 'bin') then
            call PetscViewerBinaryOpen(PETSC_COMM_WORLD,trim(adjustl(filename)), &
                                         & FILE_MODE_READ,binLoad,ierr)
            call VecLoad(vec,binLoad,ierr)
            call PetscViewerDestroy(binLoad,ierr)
            
        else if (filetype == 'txt') then
            open(15,file=trim(adjustl(filename)),status='old', action='read')
            call PetscBarrier(vec,ierr)
            line = 0
            do
               read(15,*,iostat=stat)buffer
               if (stat /= 0) exit
               line = line + 1
            end do
                
            allocate(vecArray(line))
        
            rewind(unit=15)
            do i=1,line
               read(15,*,iostat=stat) reVec, imVec
               if (stat /= 0) exit
               vecArray(i) = reVec + PETSc_i*imVec
            end do
        
            close(15)
            
            call PetscBarrier(vec,ierr)
            
            call VecGetOwnershipRange(vec,Istart,Iend,ierr)
            call VecSetValues(vec,Iend-Istart,[(i,i=Istart,Iend-1)], &
                              & vecArray(Istart+1:Iend),INSERT_VALUES,ierr)
            call VecAssemblyBegin(vec,ierr)
            call VecAssemblyEnd(vec,ierr)
            
        endif
    end subroutine importVec
    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Matrix I/O ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine exportMat(mat,filename,filetype)
        ! export a PETSc matrix to a file
        Mat, intent(in)                :: mat ! the input PETSc matrix
        character(len=50), intent(in)  :: filename ! the filename of the exported file
        character(len=3), intent(in)   :: filetype ! filetype to export ('txt' or 'bin')      
        
        ! local variables
        PetscErrorCode :: ierr
        PetscViewer    :: binSave
        
        if (filetype == 'bin') then
            call PetscViewerBinaryOpen(PETSC_COMM_WORLD,trim(adjustl(filename)), &
                                         & FILE_MODE_WRITE,binSave,ierr)
            call MatView(mat,binSave,ierr)
            call PetscViewerDestroy(binSave,ierr)
            
        else if (filetype == 'txt') then
            call PetscViewerASCIIOpen(PETSC_COMM_WORLD,trim(adjustl(filename)),binSave,ierr)
            call PetscViewerSetFormat(binSave,PETSC_VIEWER_ASCII_DENSE,ierr)
            call MatView(mat,binSave,ierr)
            call PetscViewerDestroy(binSave,ierr)
            
        else if (filetype == 'sparse') then
            call PetscViewerASCIIOpen(PETSC_COMM_WORLD,trim(adjustl(filename)),binSave,ierr)
            call PetscViewerSetFormat(binSave,PETSC_VIEWER_NATIVE,ierr)
            call MatView(mat,binSave,ierr)
            call PetscViewerDestroy(binSave,ierr)
            
        else if (filetype == 'mat') then
            call PetscViewerASCIIOpen(PETSC_COMM_WORLD,trim(adjustl(filename)),binSave,ierr)
            call PetscViewerSetFormat(binSave,PETSC_VIEWER_ASCII_MATLAB,ierr)
            call MatView(mat,binSave,ierr)
            call PetscViewerDestroy(binSave,ierr)
            
        endif
    end subroutine exportMat
    
    subroutine importMat(mat,filename,filetype)
        ! import a PETSc matrix from a file
        character(len=50), intent(in)  :: filename ! filename *only*, no path, of the input file
        character(len=3), intent(in)   :: filetype ! filetype to export ('txt' or 'bin')
        
        Mat, intent(out)                :: mat ! the PETSc matrix to import the file into
        
        ! local variables
        PetscMPIInt          :: rank
        PetscErrorCode       :: ierr
        PetscViewer          :: binLoad
        PetscInt             :: stat, line, i, Istart, Iend, j
        PetscReal            :: reVec, imVec
        PetscReal, allocatable :: matArray(:,:)
        PetscInt, allocatable :: adjArray(:,:)
        character(len=100)   :: buffer
        
        call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
        
        if (filetype == 'bin') then
            call PetscViewerBinaryOpen(PETSC_COMM_WORLD,trim(adjustl(filename)), &
                                         & FILE_MODE_READ,binLoad,ierr)
            call MatLoad(mat,binLoad,ierr)
            call PetscViewerDestroy(binLoad,ierr)
            
        else if (filetype == 'txt') then
            open(15,file=trim(adjustl(filename)),status='old', action='read')
            call PetscBarrier(mat,ierr)
            line = 0
            do
               read(15,*,iostat=stat)buffer
               if (stat /= 0) exit
               line = line + 1
            end do
                
            allocate(matArray(line,line))
        
            rewind(unit=15)
            do i=1,line
               read(15,*,iostat=stat) matArray(i,:)
               if (stat /= 0) exit
            end do
        
            close(15)
            
            call PetscBarrier(mat,ierr)
            
            call MatSetSizes(mat,PETSC_DECIDE,PETSC_DECIDE,line,line,ierr)
            call MatSetFromOptions(mat,ierr)
            call MatSetUp(mat,ierr)
            call MatGetOwnershipRange(mat,Istart,Iend,ierr)
            
            do i=Istart, Iend-1
                do j=1,line
                    if (matArray(i+1,j) .ne. 0) then
                        call MatSetValue(mat,i,j-1,matArray(i+1,j)+0*PETSc_i,INSERT_VALUES,ierr)
                    endif
                enddo
            enddo
            
            call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
            call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)
            
            deallocate(matArray)
        endif
    end subroutine importMat

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~   Some Matrix Array Fuctions  ~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine identity(mat,n)
        ! creates a :math:`n\times n` identity matrix
        PetscInt, intent(in)     :: n ! dimension of the identity matrix
        PetscScalar, intent(out) :: mat(n,n) ! the :math:`n\times n` identity matrix output
        PetscInt                 :: i

        mat = 0.d0
        do i = 1, n
            mat(i,i)=1.d0
        end do
    end subroutine identity
    
    subroutine kron(M1,r1,c1, M2,r2,c2, kronProd)
        ! compute the Kronecker product of two arrays
        PetscInt, intent(in)    :: r1 ! number of rows in matrix 1
        PetscInt, intent(in)    :: c1 ! number of columns in matrix 1
        PetscInt, intent(in)    :: r2 ! number of rows in matrix 2
        PetscInt, intent(in)    :: c2 ! number of columns in matrix 2
        PetscScalar, intent(in) :: M1(r1,c1) ! matrix 1
        PetscScalar, intent(in) :: M2(r2,c2) ! matrix 2
        PetscScalar, intent(out):: kronProd(r1*r2,c1*c2) ! kronecker product output

        ! local variables
        PetscInt :: i, j
        
        kronProd = 0.d0
        forall (i=1:r1, j=1:c1)
            kronProd(r2*(i-1)+1:r2*i, c2*(j-1)+1:c2*j) = M1(i,j)*M2
        end forall
    end subroutine kron

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~Import and convert Adjacency matrix to Hamiltonian ~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine importAdjToH(mat,filename,p,d,amp,interaction,nd)
        ! Import an adjacency matrix from a file, and create a PETSc Hamiltonian matrix.
        character, intent(in)          :: p ! number of particles in the system ('1', '2', or '3')
        character(len=50), intent(in)  :: filename ! dense adjacency matrix file in text format
        PetscInt, intent(in)           :: nd ! number of defects
        PetscInt, intent(in)           :: d(nd) ! nodes to place defects on
        PetscScalar, intent(in)        :: amp(nd) ! amplitude of defects
        PetscScalar, intent(in)        :: interaction ! interaction amplitude
        
        Mat, intent(out)                :: mat ! output Hamiltonian matrix
        
        ! local variables
        PetscMPIInt              :: rank
        PetscErrorCode           :: ierr
        PetscInt                 :: stat, N, i
        PetscInt, allocatable    :: adjArray(:,:)
        PetscScalar, allocatable :: HArray(:,:)
        character(len=100)       :: buffer
        
        call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
        
        open(15,file=trim(adjustl(filename)),status='old', action='read')
        call PetscBarrier(mat,ierr)
        N = 0
        do
           read(15,*,iostat=stat)buffer
           if (stat /= 0) exit
           N = N + 1
        end do

        N = N+1
        
        allocate(adjArray(N,N),HArray(N,N))

        rewind(unit=15)
        do i=1,N
           read(15,*,iostat=stat) adjArray(i,:)
           if (stat /= 0) exit
        end do

        close(15)

        call adjToH(mat,adjArray,p,d,amp,interaction,nd,N)
        deallocate(adjArray,HArray)

    end subroutine importAdjToH

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~ Convert Adjacency array to Hamiltonian ~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine adjToH(mat,adjArray,p,d,amp,interaction,nd,N)
        ! convert an adjacency array to a PETSc Hamiltonian matrix
        character, intent(in)          :: p ! number of particles in the system ('1', '2', or '3')
        PetscInt, intent(in)           :: nd ! number of defects
        PetscInt, intent(in)           :: d(nd) ! nodes to place defects on
        PetscInt, intent(in)           :: N ! number of vertices in the graph
        PetscInt, intent(in)           :: adjArray(N,N) ! nxn array containing adjacency matrix
        PetscScalar, intent(in)        :: amp(nd) ! amplitude of defects
        PetscScalar, intent(in)        :: interaction ! interaction amplitude
        
        Mat, intent(out)                :: mat ! output Hamiltonian matrix
        
        ! local variables
        PetscMPIInt              :: rank
        PetscErrorCode           :: ierr
        PetscInt                 :: stat, i, Istart, Iend, j, k, Nsq, Ncb
        PetscScalar              :: HArray(N,N)
        PetscInt, allocatable    :: intTestArray(:)
        PetscScalar, allocatable :: ident(:,:),ident2(:,:), temp(:,:),temp3(:,:), H2Array(:,:)
        
        call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
        
        HArray = 0.d0
        do i=1, N            
            do j=1,N
               if (i==j) then
                   do k=1,nd
                       if (d(k)==i-1) then
                           HArray(i,j) = amp(k)+sum(adjArray(i,:))-adjArray(i,j)
                           exit
                       else
                           HArray(i,j) = sum(adjArray(i,:))-adjArray(i,j)
                       endif
                   enddo
                else
                    HArray(i,j) = -adjArray(i,j)
                endif
            enddo
        enddo
        
        if (p == '2') then
            Nsq = N*N
            allocate(H2Array(Nsq,Nsq), ident(N,N), temp(Nsq,Nsq))
        
            call identity(ident,N)
            call kron(HArray,N,N,ident,N,N,H2Array)
            call kron(ident,N,N,HArray,N,N,temp)

            H2Array = H2Array + temp
            deallocate(ident,temp)
            
            
            call PetscBarrier(mat,ierr)
        
            call MatSetSizes(mat,PETSC_DECIDE,PETSC_DECIDE,Nsq,Nsq,ierr)
            call MatSetFromOptions(mat,ierr)
            call MatSetUp(mat,ierr)
            call MatGetOwnershipRange(mat,Istart,Iend,ierr)
        
            do i=Istart, Iend-1            
                do j=1,Nsq
                    if (H2Array(i+1,j) .ne. 0) then
                        if (i==j-1 .and. mod(j,n+1)==1) then
                            call MatSetValue(mat,i,j-1,H2Array(i+1,j)+interaction,INSERT_VALUES,ierr)
                        else
                            call MatSetValue(mat,i,j-1,H2Array(i+1,j),INSERT_VALUES,ierr)
                        endif
                    endif
                enddo
            enddo
        
            call PetscBarrier(mat,ierr)

            call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
            CHKERRQ(ierr)
            call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)
            CHKERRQ(ierr)
        
            deallocate(H2Array)

        elseif (p == '3') then
            Nsq = N*N
            Ncb = N*Nsq

            allocate(H2Array(Ncb,Ncb), ident(N,N), ident2(Nsq,Nsq), &
                       & intTestArray(N),temp(Ncb,Ncb),temp3(Ncb,Ncb))
        
            call identity(ident,N)
            call identity(ident2,Nsq)

            call kron(ident2,Nsq,Nsq,HArray,N,N,temp)
            call kron(HArray,N,N,ident2,Nsq,Nsq,temp3)

            call kron(ident,N,N,HArray,N,N,ident2)
            call kron(ident2,Nsq,Nsq,ident,N,N,H2Array)

            H2Array = temp + H2Array + temp3
            deallocate(ident,ident2,temp,temp3)
            
            call PetscBarrier(mat,ierr)
        
            call MatSetSizes(mat,PETSC_DECIDE,PETSC_DECIDE,Ncb,Ncb,ierr)
            call MatSetFromOptions(mat,ierr)
            call MatSetUp(mat,ierr)
            call MatGetOwnershipRange(mat,Istart,Iend,ierr)

            intTestArray = [(i, i=1, Nsq, n)]

            do i=Istart, Iend-1
                do j=1,Ncb
                    if (H2Array(i+1,j) .ne. 0) then 
                        if (i==j-1) then
                            if ( mod(j,Nsq+n+1) == 1 ) then
                                 call MatSetValue(mat,i,j-1,H2Array(i+1,j)+2*interaction,INSERT_VALUES,ierr)
                            elseif ((mod(j,Nsq+n) .le. n) .and. (mod(j,Nsq+1) .ge. 1)) then
                                 call MatSetValue(mat,i,j-1,H2Array(i+1,j)+interaction,INSERT_VALUES,ierr)
                            elseif ( mod(j,n+1) == int(ceiling(real(j)/real(Nsq))) ) then
                                 call MatSetValue(mat,i,j-1,H2Array(i+1,j)+interaction,INSERT_VALUES,ierr)
                            elseif ( any(intTestArray .eq. mod(j,Nsq+1)) ) then
                                 call MatSetValue(mat,i,j-1,H2Array(i+1,j)+interaction,INSERT_VALUES,ierr)
                            else
                                 call MatSetValue(mat,i,j-1,H2Array(i+1,j),INSERT_VALUES,ierr)
                            endif
                        else
                             call MatSetValue(mat,i,j-1,H2Array(i+1,j),INSERT_VALUES,ierr)
                        endif
                    endif
                enddo
            enddo

            call PetscBarrier(mat,ierr)
        
            call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
            CHKERRQ(ierr)
            call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)
        
            deallocate(H2Array,intTestArray)
        
        else
            call PetscBarrier(mat,ierr)
        
            call MatSetSizes(mat,PETSC_DECIDE,PETSC_DECIDE,N,N,ierr)
            call MatSetFromOptions(mat,ierr)
            call MatSetUp(mat,ierr)
            call MatGetOwnershipRange(mat,Istart,Iend,ierr)
        
            do i=Istart, Iend-1            
                do j=1,N
                    if (HArray(i+1,j) .ne. 0) then
                        call MatSetValue(mat,i,j-1,HArray(i+1,j),INSERT_VALUES,ierr)
                    endif
                enddo
            enddo
        
            call PetscBarrier(mat,ierr)
            call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
            call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)
        endif

    end subroutine adjToH 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~ Convert from 2D to 1D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~ statespace for 2 particles ~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function coord(x,y,n)
        ! convert from a 2D to 1D statespace for 2 particles, using
        ! :math:`coord = n(x + n/2 - 1) + y + n/2 - 1`
        integer, intent(in)    :: n ! number of nodes in the system
        integer, intent(in)    :: x ! vertex location of particle 1
        integer, intent(in)    :: y ! vertex location of particle 2
        integer                :: coord   ! output 2P statepace coordinate
    
        coord = n*(x + n/2 - 1) + y + n/2 - 1
    end function coord

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~ Convert from 3D to 1D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~ statespace for 3 particles ~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function coord3P(x,y,z,n)
        ! convert from a 2D to 1D statespace for 2 particles, using
        ! :math:`coord2P = xn^2 + ny + z`
        integer, intent(in)    :: n ! number of nodes in the system
        integer, intent(in)    :: x ! vertex location of particle 1
        integer, intent(in)    :: y ! vertex location of particle 2
        integer, intent(in)    :: z ! vertex location of particle 3
        integer                :: coord3P ! output 3P statepace coordinate
    
        coord3P = (n*n)*x + n*y + z
    end function coord3P
    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~ 1P  probabilities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine marginal1(psi,prob,n)
        ! Calculates the marginal probability of 1 particle state psi; i.e. 
        ! :math:`|\psi|^2`
        PetscInt, intent(in)      :: n ! length of vector psi
        Vec, intent(in)           :: psi ! input statespace PETSc vector
        
        Vec, intent(out)          :: prob ! output PETSc vector containing the probabilities
        
        ! local variables
        PetscErrorCode :: ierr
        PetscScalar, pointer :: probArray(:)
        
        ! copy psi to prob
        call VecCopy(psi,prob,ierr)
        call VecAbs(prob,ierr)
        
        call VecGetArrayF90(prob,probArray,ierr)
        probArray = probArray**2.d0
        call VecRestoreArrayF90(prob,probArray,ierr)
    end subroutine marginal1

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~ marginal probabilities~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine marginal2(psi,psiM,p,n)
        ! Calculates the marginal probability of particle number :f:var:`p`
        character, intent(in)     :: p ! the particle to calculate the marginal probability ('1', '2', '3')
        PetscInt, intent(in)      :: n ! number of vertices in the graph
        Vec, intent(in)           :: psi ! input PETSc statespace vector of length :math:`n^p`
        
        Vec, intent(out)          :: psiM ! output PETSc vector containing the probabilities
        
        ! local variables
        PetscErrorCode :: ierr
        PetscMPIInt    :: rank
        PetscInt    :: i, j, NN, Istart, Iend
        PetscInt, allocatable :: ind(:)
        PetscScalar, allocatable :: temp(:)
        PetscScalar, pointer :: workArray(:)
        Vec         :: work
        VecScatter  :: ctx
        
        NN = N*N
        
        ! create work vector
        call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,NN,work,ierr)
        call VecSetFromOptions(work,ierr)
        
        ! Vector scatter to all processes
        call VecScatterCreateToAll(psi,ctx,work,ierr)
        call VecScatterBegin(ctx,psi,work,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(ctx,psi,work,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterDestroy(ctx,ierr)
        
        call VecGetArrayF90(work,workArray,ierr)
        
        call VecGetOwnershipRange(psiM,Istart,Iend,ierr)
        allocate(ind(Iend-Istart),temp(Iend-Istart))
        ind = [(i,i=Istart,Iend-1)]
        
        call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
        
        if (p=='x') then
            do i=Istart,Iend-1
                temp(i+1-Istart) = sum(abs(workArray(1+(i+1-1)*n:(i+1)*n))**2.d0)
            end do
        elseif (p=='y') then
            do i=Istart,Iend-1
                temp(i+1-Istart) = sum(abs([(workArray(i+1+j*n), j=0, n-1)])**2.d0)
            end do
        endif
        
        call VecRestoreArrayF90(work,workArray,ierr)
        call VecDestroy(work,ierr)
        
        call VecSetValues(psiM,size(ind),ind,temp,INSERT_VALUES,ierr)
        call VecAssemblyBegin(psiM,ierr)
        call VecAssemblyEnd(psiM,ierr)
        deallocate(ind,temp)
    end subroutine marginal2

    subroutine marginal3(psi,psiM,p,n)
        character, intent(in)     :: p
        PetscInt, intent(in)      :: n
        Vec, intent(in)           :: psi
        
        Vec, intent(out)          :: psiM
        
        ! local variables
        PetscErrorCode           :: ierr
        PetscMPIInt              :: rank
        PetscInt                 :: i, j, NN, Nsq, Istart, Iend
        PetscInt, allocatable    :: ind(:)
        PetscScalar, allocatable :: temp(:)
        PetscScalar, pointer     :: workArray(:)
        Vec                      :: work
        VecScatter               :: ctx
        
        Nsq = N*N
        NN = N*Nsq
        
        ! create work vector
        call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,NN,work,ierr)
        call VecSetFromOptions(work,ierr)
        
        ! Vector scatter to all processes
        call VecScatterCreateToAll(psi,ctx,work,ierr)
        call VecScatterBegin(ctx,psi,work,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(ctx,psi,work,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterDestroy(ctx,ierr)
        
        call VecGetArrayF90(work,workArray,ierr)
        
        call VecGetOwnershipRange(psiM,Istart,Iend,ierr)
        allocate(ind(Iend-Istart),temp(Iend-Istart))
        ind = [(i,i=Istart,Iend-1)]
        
        call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
        
        if (p=='x') then
            do i=Istart,Iend-1
                temp(i+1-Istart) = sum(abs(workArray(1+i*(Nsq):(Nsq)*(1+i)))**2.d0)
            end do
        elseif (p=='y') then
            do i=Istart,Iend-1
                temp(i+1-Istart) = sum(abs([(workArray(1+i*n+j*(Nsq):n+i*n+j*(Nsq)), j=0, n-1)])**2.d0)
            end do
        elseif (p=='z') then
            do i=Istart,Iend-1
                temp(i+1-Istart) = sum(abs(workArray(1+i::n))**2.d0)
            end do
        endif
        
        call VecRestoreArrayF90(work,workArray,ierr)
        call VecDestroy(work,ierr)
        
        call VecSetValues(psiM,size(ind),ind,temp,INSERT_VALUES,ierr)
        call VecAssemblyBegin(psiM,ierr)
        call VecAssemblyEnd(psiM,ierr)
        deallocate(ind,temp)
    end subroutine marginal3
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
!~~~~~~~~~~~~~~~~~~~~~ create 1P line state space ~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    subroutine p1_init(psi0,init_state,num,N)
        PetscInt, intent(in)     :: num, N
        PetscScalar, intent(in)  :: init_state(num,2)
        
        Vec, intent(out)         :: psi0
        
        ! local variables
        PetscErrorCode :: ierr
        PetscInt       :: i, ind(num)
        PetscScalar    :: val(num)
        
        do i=1, num
            ind(i) = int(init_state(i,1))+N/2-1
            val(i) = init_state(i,2)
        end do
        
        call VecSetValues(psi0,num,ind,val,INSERT_VALUES,ierr)
        call VecAssemblyBegin(psi0,ierr)
        call VecAssemblyEnd(psi0,ierr)    
    end subroutine p1_init
    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
!~~~~~~~~~~~~~~~~~~~~~ create 2P line state space ~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    subroutine p2_init(psi0,init_state,num,N)
        PetscInt, intent(in)     :: num, N
        PetscScalar, intent(in)  :: init_state(num,3)
        
        Vec, intent(out)         :: psi0
        
        ! local variables
        PetscErrorCode :: ierr
        PetscInt    :: i, ind(num)
        PetscScalar :: val(num)
        
        do i=1, num
            ind(i) = coord(int(init_state(i,1)),int(init_state(i,2)),N)
            val(i) = init_state(i,3)
        end do
        
        call VecSetValues(psi0,num,ind,val,INSERT_VALUES,ierr)
        call VecAssemblyBegin(psi0,ierr)
        call VecAssemblyEnd(psi0,ierr)    
    end subroutine p2_init

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
!~~~~~~~~~~~~~~~~~~~~~ create 3P line state space ~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    subroutine p3_init(psi0,init_state,num,N)
        PetscInt, intent(in)     :: num, N
        PetscScalar, intent(in)  :: init_state(num,4)
        
        Vec, intent(out)         :: psi0
        
        ! local variables
        PetscErrorCode :: ierr
        PetscInt    :: i, ind(num)
        PetscScalar :: val(num)
        
        do i=1, num
            ind(i) = coord3P(int(init_state(i,1)),int(init_state(i,2)),int(init_state(i,3)),N)
            val(i) = init_state(i,4)
        end do
        
        call VecSetValues(psi0,num,ind,val,INSERT_VALUES,ierr)
        call VecAssemblyBegin(psi0,ierr)
        call VecAssemblyEnd(psi0,ierr)    
    end subroutine p3_init
    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
!~~~~~~~~~~~~~~ create a sparse hamiltonian matrix of size n ~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~ using PETSc's sparse matrix routines ~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine hamiltonian_p1_line(A,d,amp,nd,n)
        PetscInt, intent(in)    :: nd, n, d(nd)
        PetscScalar, intent(in) :: amp(nd)
        
        Mat, intent(out)        :: A
        
        ! local variables
        PetscErrorCode :: ierr
        PetscInt       :: i, j, its, Istart, Iend, col(3)
        PetscScalar    :: value(3)
        
        call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
        call MatSetFromOptions(A,ierr)
        call MatSetUp(A,ierr)
        
        call MatGetOwnershipRange(A,Istart,Iend,ierr)
        
        if (Istart == 0) then
            col(1:2) = [0,1]
            value(1:2) = [2.d0,-1.d0]
            
            do j = 1, nd
                if (int(d(j)+n/2.0) == 1) then
                    value(1) = 2.d0+amp(j)
                    exit
                endif
            end do
            
            call MatSetValues(A,1,0,2,col,value,INSERT_VALUES,ierr)
            Istart = Istart + 1
        end if
        
        if (Iend == n) then
            col(1:2) = [n-2,n-1]
            value(1:2) = [-1.d0,2.d0]
            
            do j = 1, nd
                if (int(d(j)+n/2.0) == n) then
                    value(2) = 2.d0+amp(j)
                    exit
                endif
            end do
            
            call MatSetValues(A,1,n-1,2,col,value,INSERT_VALUES,ierr)
            Iend = Iend - 1
        endif
        
        value(1:3) = [-1.d0,2.d0,-1.d0]
        do i=Istart, Iend-1
            col = [i-1,i,i+1]
            
            value = [-1.d0,2.d0,-1.d0]
            do j = 1, nd
                if (int(d(j)+n/2.0) == i+1) then
                    value(2) = 2.d0+amp(j)
                    exit
                endif
            end do
            
            call MatSetValues(A,1,i,3,col,value,INSERT_VALUES,ierr)
        enddo
        
        call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    
    end subroutine hamiltonian_p1_line

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!~~~~~~~~~ create a sparse 2P hamiltonian matrix of size n^2 x n^2 ~~~~~~~
!~~~~~~~~~~~~~~~~~~ using PETSc's sparse matrix routines ~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine hamiltonian_p2_line(H2,d,amp,interaction,nd,n)
        PetscInt, intent(in)    :: n, nd, d(nd)
        PetscScalar, intent(in) :: amp(nd), interaction
        
        Mat, intent(out)        :: H2
        
        ! local variables
        PetscErrorCode :: ierr
        PetscInt       :: NN, i, j, its, Istart, Iend, col(3)
        PetscScalar    :: value(3), alpha
        PetscScalar    :: diagArray(n)
        
        NN = N*N
        
        ! create array of 1D diagonal entries
        diagArray = 2.d0
        do i = 1, n
            do j = 1, nd
                if (int(d(j)+n/2.0) == i) then
                    diagArray(i) = 2.d0+amp(j)
                    exit
                endif
            end do
        end do
        
        ! create H2 matrix
        call MatSetSizes(H2,PETSC_DECIDE,PETSC_DECIDE,NN,NN,ierr)
        call MatSetFromOptions(H2,ierr)
        call MatSetUp(H2,ierr)
        
        ! create KronProd(H1,I) matrix and store in H2
        call MatGetOwnershipRange(H2,Istart,Iend,ierr)
        
        if (Istart == 0) then
            col(1:2) = [0,1]
            value(1:2) = [diagArray(1)+interaction,-1.d0+0*PETSC_i]            
            call MatSetValues(H2,1,0,2,col,value,INSERT_VALUES,ierr)
            Istart = Istart + 1
        end if
        
        if (Iend == NN) then
            col(1:2) = [NN-2,NN-1]
            value(1:2) = [-1.d0+0*PETSC_i,diagArray(n)+interaction]            
            call MatSetValues(H2,1,NN-1,2,col,value,INSERT_VALUES,ierr)
            Iend = Iend - 1
        endif
        
        do i=Istart, Iend-1
            if (mod(i,n)+1==1) then

                if (mod(i,n+1)==0) then
                    alpha = interaction
                else
                    alpha = 0.d0
                endif
            
                col(1:2) = [i,i+1]
                value(1:2) = [diagArray(mod(i,n)+1)+alpha,-1.d0+0*PETSC_i]            
                call MatSetValues(H2,1,i,2,col,value,INSERT_VALUES,ierr)
            
            elseif (mod(i,n)+1==n) then

                if (mod(i,n+1)==0) then
                    alpha = interaction
                else
                    alpha = 0.d0
                endif
            
                col(1:2) = [i-1,i]
                value(1:2) = [-1.d0+0*PETSC_i,diagArray(mod(i,n)+1)+alpha]            
                call MatSetValues(H2,1,i,2,col,value,INSERT_VALUES,ierr)
            
            else

                if (mod(i,n+1)==0) then
                    alpha = interaction
                else
                    alpha = 0.d0
                endif

                col = [i-1,i,i+1]
                value = [-1.d0+0*PETSC_i,diagArray(mod(i,n)+1)+alpha,-1.d0+0*PETSC_i]            
                call MatSetValues(H2,1,i,3,col,value,INSERT_VALUES,ierr)
            endif
        enddo
        
        call MatAssemblyBegin(H2,MAT_FLUSH_ASSEMBLY,ierr)
        call MatAssemblyEnd(H2,MAT_FLUSH_ASSEMBLY,ierr)
        
        ! create KronProd(I,H1) matrix and add it to H2
        call MatGetOwnershipRange(H2,Istart,Iend,ierr)
        
        do while (Istart<n)
            col(1:2) = [Istart,Istart+n]
            value(1:2) = [diagArray(1),-1.d0+0*PETSC_i]            
            call MatSetValues(H2,1,Istart,2,col,value,ADD_VALUES,ierr)
            Istart = Istart + 1
        end do
        
        do while (Iend>NN-n)
            col(1:2) = [Iend-1-n,Iend-1]
            value(1:2) = [-1.d0+0*PETSC_i,diagArray(n)]            
            call MatSetValues(H2,1,Iend-1,2,col,value,ADD_VALUES,ierr)
            Iend = Iend - 1
        end do
        
        do i=Istart, Iend-1
            col = [i-n,i,i+n]
            value = [-1.d0+0*PETSC_i,diagArray(int(real(i)/n)+1),-1.d0+0*PETSC_i]
            call MatSetValues(H2,1,i,3,col,value,ADD_VALUES,ierr)
        enddo
        
        call MatAssemblyBegin(H2,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(H2,MAT_FINAL_ASSEMBLY,ierr)
    
    end subroutine hamiltonian_p2_line

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!~~~~~~~~~ create a sparse 3P hamiltonian matrix of size n^3 x n^3 ~~~~~~~
!~~~~~~~~~~~~~~~~~~ using PETSc's sparse matrix routines ~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine hamiltonian_p3_line(H3,d,amp,interaction,nd,n)
        PetscInt, intent(in)    :: n, nd, d(nd)
        PetscScalar, intent(in) :: amp(nd), interaction
        
        Mat, intent(out)        :: H3
        
        ! local variables
        PetscErrorCode :: ierr
        PetscInt       :: i, adjArray(n,n)
        
        ! create adjacency matrix
        adjArray = 0
        adjArray(1,2) = 1
        adjArray(n,n-1) = 1
        do i = 2, n-1
            adjArray(i,i-1) = 1
            adjArray(i,i+1) = 1
        enddo
                
        call adjToH(H3,adjArray,'3',d,amp,interaction,nd,N) 
    
    end subroutine hamiltonian_p3_line
    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
!~~~~~~~~~~~~~~~~~~~~~ calculate y=e^(-iHt).v using SLEPc ~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine qw_krylov(A,t,v,y)
        Vec, intent(in)         :: v
        Mat, intent(in)         :: A
        PetscScalar, intent(in) :: t
        
        Vec, intent(out)        :: y
        
        ! local variables
        character(len=40)       :: vecType
        PetscErrorCode          :: ierr
        PetscScalar             :: alpha
        MFN                     :: mfn

        ! check if vector is sequential or MPI
        call VecGetType(v,vecType,ierr)
        vecType = trim(adjustl(vecType))

        call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)

        if (vecType == "seq") then
            call MFNCreate(MPI_COMM_SELF, mfn, ierr)
        else
            call MFNCreate(PETSC_COMM_WORLD, mfn, ierr)
        endif

        call MFNSetOperator(mfn,A,ierr)
        call MFNSetFunction(mfn,SLEPC_FUNCTION_EXP,ierr)
        call MFNSetFromOptions(mfn,ierr)
        
        alpha = -t*PETSC_i
        call MatScale(A,alpha,ierr)
        !call MFNSetScaleFactor(mfn, alpha,ierr)
        call MFNSolve(mfn,v,y,ierr)
        call MatScale(A,1.d0/alpha,ierr)
        
        call MFNDestroy(mfn,ierr)
        call SlepcFinalize(ierr)
    
    end subroutine qw_krylov

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
!~~~~~~~~~~~~ calculate the min or max eiganvalue using SLEPc ~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

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
        character(len=40) :: matType
        PetscErrorCode    :: ierr
        PetscInt          :: i, j, its, nev, maxit, nconv
        PetscScalar       :: kr, ki
        PetscReal         :: tol
        EPS               :: eps
        EPSType           :: tname
        ST                :: st
        character(len=12) :: arg
        Vec               :: vec
        
        ! check if vector is sequential or MPI
        call MatGetType(A,matType,ierr)
        matType = trim(adjustl(matType))

        ! initialise SLEPc and create the EPS
        call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
        if (matType == "seqaij") then
            call EPSCreate(MPI_COMM_SELF,eps,ierr)
        else
            call EPSCreate(PETSC_COMM_WORLD,eps,ierr)
        endif
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
            case('power');     call EPSSetType(eps,EPSPOWER,ierr)
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
                !if (rank==0) write(*,*)'Calculating Emax...'
                call EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL,ierr)
            case default
                !if (rank==0) write(*,*)'Calculating Emin...'
                if (arg=='power') then
                    call EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE,ierr)
                else
                    call EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL,ierr)
                endif
                
                call EPSGetST(eps,st,ierr)
                call STSetType(st,STSINVERT,ierr)
                call STSetShift(st,(0.d0,0.d0),ierr)
                
                !call MatGetVecs(A,vec,PETSC_NULL_OBJECT,ierr)
                !call VecSet(vec,1.0,ierr)
                !call EPSSetDeflationSpace(eps,1,vec,ierr)
                !call VecDestroy(vec,ierr)
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
        call SlepcFinalize(ierr)
        
    end subroutine min_max_eigs

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Chebyshev method ~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    subroutine qw_cheby(psi0,psi,dt,H,Emin,Emax,rank,N)
        PetscMPIInt, intent(in) :: rank
        PetscInt, intent(in)    :: N
        PetscScalar, intent(in) :: dt, Emin, Emax
        Vec, intent(in)         :: psi0
        Mat, intent(in)         :: H
        
        Vec, intent(out)        :: psi
        
        ! local variables
        PetscErrorCode :: ierr
        PetscInt      :: m, terms, i, j
        PetscReal     :: alpha
        PetscScalar   :: bessj0, bessj1, bessjn, &
                            EmEm, d2EmEm, imagM, neg1
        Vec, pointer  :: work(:)

        EmEm = (Emax+Emin)/(Emax-Emin)
        d2EmEm = -2.d0/(Emax-Emin)
        alpha = PetscRealPart((Emax-Emin)*dt/2.d0)
        
        call VecDuplicateVecsF90(psi0,4,work,ierr)
        
        call VecCopy(psi0,work(1),ierr)
        call MatMult(H,work(1),work(2),ierr)
        call VecAXPBY(work(2), EmEm,d2EmEm, work(1),ierr)
        
        bessj0 = dbesjn(0,alpha)
        bessj1 = dbesjn(1,alpha)
        call VecCopy(psi0,work(4),ierr)
        call VecAXPBY(work(4), 2.d0*PETSC_i*bessj1,bessj0, work(2),ierr)
        
        terms = 0
        do while (abs(2.d0*dbesjn(terms,alpha)) > 1.d-18)
            terms = terms + 1
        end do
        
        EmEm = 2.d0*EmEm
        d2EmEm = 2.d0*d2EmEm
        imagM = 2.d0*(PETSC_i*PETSC_i)
        neg1 = -1.d0

        do m = 2, terms
            call MatMult(H,work(2),work(3),ierr)
            call VecAXPBY(work(3), EmEm, d2EmEm, work(2),ierr)
            call VecAXPY(work(3),neg1,work(1),ierr)
            
            bessjn = dbesjn(m,alpha)
            call VecAXPY(work(4),imagM*bessjn,work(3),ierr)
            imagM = imagM*PETSC_i

            call VecCopy(work(2),work(1),ierr)
            call VecCopy(work(3),work(2),ierr)
        end do
        
        call VecScale(work(4),exp(-PETSC_i*(Emax+Emin)*dt/2.d0),ierr)
        call VecCopy(work(4),psi,ierr)
        
        call VecDestroyVecsF90(4,work,ierr)
    
    end subroutine qw_cheby

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2P Entanglement ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    subroutine partial_trace_array(psi,rhoX,n)
        Vec, intent(in)          :: psi
        PetscInt, intent(in)     :: n

        PetscScalar, intent(out) :: rhoX(0:n-1,0:n-1)

        ! local variables
        character(len=40)        :: vecType
        PetscErrorCode           :: ierr
        PetscMPIInt              :: rank
        PetscInt                 :: v, k, i, rhoLength
        PetscScalar, pointer     :: workArray(:)
        Vec                      :: work
        VecScatter               :: ctx

        call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)

        ! check if vector is sequential or MPI
        call VecGetType(psi,vecType,ierr)
        vecType = trim(adjustl(vecType))
        
        call VecGetSize(psi,rhoLength,ierr)

        if (vecType == "seq") then
            call VecGetArrayF90(psi,workArray,ierr)
        else
            ! create work vector
            call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,rhoLength,work,ierr)
            call VecSetFromOptions(work,ierr)
            
            ! Vector scatter to all processes
            call VecScatterCreateToAll(psi,ctx,work,ierr)
            call VecScatterBegin(ctx,psi,work,INSERT_VALUES,SCATTER_FORWARD,ierr)
            call VecScatterEnd(ctx,psi,work,INSERT_VALUES,SCATTER_FORWARD,ierr)
            call VecScatterDestroy(ctx,ierr)
            
            call VecGetArrayF90(work,workArray,ierr)
        endif

        ! calculate the partial trace of *last* particle
        rhoX = 0.d0
        do v=0, n-1
            do k=v, n-1
                do i=0, rhoLength-1-k+v
                    if (mod(i,n) == v) then
                        rhoX(v,k) = rhoX(v,k) + workArray(i)*conjg(workArray(i+k-v))

                        if (k.ne.v) rhoX(k,v) = rhoX(k,v) + conjg(workArray(i))*workArray(i+k-v)
                    endif
                enddo
            enddo
        enddo
        
        if (vecType == "seq") then
            call VecRestoreArrayF90(psi,workArray,ierr)
        else
            call VecRestoreArrayF90(work,workArray,ierr)
            call VecDestroy(work,ierr)
        endif

    end subroutine partial_trace_array

    subroutine partial_trace_mat(psi,rhoX,n)
        Vec, intent(in)          :: psi
        PetscInt, intent(in)     :: n

        Mat, intent(out)         :: rhoX

        ! local variables
        character(len=40)        :: vecType
        PetscErrorCode           :: ierr
        PetscMPIInt              :: rank
        PetscInt                 :: v, k, i, Istart, Iend, rhoLength
        PetscScalar              :: temp(0:n-1)
        PetscScalar, pointer     :: workArray(:)
        Vec                      :: work
        VecScatter               :: ctx

        call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)

        ! check if vector is sequential or MPI
        call VecGetType(psi,vecType,ierr)
        vecType = trim(adjustl(vecType))

        call VecGetSize(psi,rhoLength,ierr)

        if (vecType == "seq") then
            call VecGetArrayF90(psi,workArray,ierr)
        else
            ! create work vector
            call VecGetSize(psi,rhoLength,ierr)
            call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,rhoLength,work,ierr)
            call VecSetFromOptions(work,ierr)
            
            ! Vector scatter to all processes
            call VecScatterCreateToAll(psi,ctx,work,ierr)
            call VecScatterBegin(ctx,psi,work,INSERT_VALUES,SCATTER_FORWARD,ierr)
            call VecScatterEnd(ctx,psi,work,INSERT_VALUES,SCATTER_FORWARD,ierr)
            call VecScatterDestroy(ctx,ierr)
            
            call VecGetArrayF90(work,workArray,ierr)
        endif

        ! set up partial trace matrix
        call MatSetSizes(rhoX,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
        call MatSetOption(rhoX,MAT_HERMITIAN,PETSC_TRUE,ierr)
        call MatSetFromOptions(rhoX,ierr)
        call MatSetUp(rhoX,ierr)
        call MatGetOwnershipRange(rhoX,Istart,Iend,ierr)

        ! calculate the partial trace of the *last* particle
        do v=Istart, Iend-1
            temp = 0.d0

            do k=v, n-1
                do i=0, rhoLength-1-k+v
                    if (mod(i,n) == v) then
                        temp(k) = temp(k) + workArray(i)*conjg(workArray(i+k-v))
                    endif
                enddo
            enddo

            call MatSetValues(rhoX,1,v,n-v,[(k,k=v,n-1)],temp(v:n-1),INSERT_VALUES,ierr)
            call MatSetValues(rhoX,n-v-1,[(k,k=v+1,n-1)],1,v,conjg(temp(v+1:n-1)),INSERT_VALUES,ierr)
        enddo

        if (vecType == "seq") then
            call PetscBarrier(PETSC_COMM_SELF,ierr)
        else
            call PetscBarrier(PETSC_COMM_WORLD,ierr)
        endif

        call MatAssemblyBegin(rhoX,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(rhoX,MAT_FINAL_ASSEMBLY,ierr)
        
        if (vecType == "seq") then
            call VecRestoreArrayF90(psi,workArray,ierr)
        else
            call VecRestoreArrayF90(work,workArray,ierr)
            call VecDestroy(work,ierr)
        endif

    end subroutine partial_trace_mat

    subroutine entanglement(psi,n,vNE,eig_solver,worktype,worktypeInt,tolIn,max_it,verbose,error)
        Vec, intent(in)              :: psi
        character(len=*),intent(in)  :: eig_solver
        PetscInt, intent(in)         :: worktypeInt, max_it, n
        PetscReal, intent(in)        :: tolIn
        character(len=3), intent(in) :: worktype
        PetscBool, intent(in)        :: verbose

        PetscReal, intent(out)       :: vNE
        PetscInt, intent(out)        :: error

        ! local variables
        character(len=40)        :: vecType
        PetscScalar              :: lambda(n), kr, ki
        PetscInt                 :: i, its, nev, maxit, nconv
        PetscReal                :: tol
        PetscErrorCode           :: ierr
        PetscMPIInt              :: rank
        Mat                      :: rhoX
        EPS                      :: eps
        EPSType                  :: tname
        character(len=12)        :: arg

        call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)

        ! check if vector is sequential or MPI
        call VecGetType(psi,vecType,ierr)
        vecType = trim(adjustl(vecType))

        ! initialise slepc
        call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)

        if (vecType == "seq") then
            call MatCreate(MPI_COMM_SELF,rhoX,ierr)
            call EPSCreate(MPI_COMM_SELF,eps,ierr)
        else
            call MatCreate(PETSC_COMM_WORLD,rhoX,ierr)
            call EPSCreate(PETSC_COMM_WORLD,eps,ierr)
        endif

        ! get partial trace
        call partial_trace_mat(psi,rhoX,n)

        ! set the problem type as EPS_HEP: Hermitian matrix
        call EPSSetOperators(eps,rhoX,PETSC_NULL_OBJECT,ierr)
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
            case('krylovschur'); call EPSSetType(eps,EPSKRYLOVSCHUR,ierr)
            case('lapack');      call EPSSetType(eps,EPSLAPACK,ierr)
            case default;        continue ! SLEPC default is Krylov-Schur
        end select   

        ! select all eigenvalues
        call EPSSetWhichEigenpairs(eps,EPS_ALL,ierr)
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

        lambda = 0.d0


        ! determine if convergence occurs
        call EPSGetConverged(eps,nconv,ierr)
        if (nconv==n) then
            do i=1,n
                call EPSGetEigenpair(eps,i-1,kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
                lambda(i) = PetscRealPart(kr)
            enddo
            error = 0
        else
            error = 1
            if (rank==0) write(*,*)'Not all eigenvalues converged. Reduce tolerance, increase &
                        & number of interations, or increase work vector size'
        endif

        ! calculate the von Neumann Entropy
        vNE = 0.d0
        if (error == 0) then
            do i=1, n
                if (lambda(i) .ne. 0.d0) then
                    vNE = vNE - lambda(i)*log(lambda(i))/log(2.d0)
                endif
            enddo
        endif

        ! clean up
        call EPSDestroy(eps,ierr)
        call MatDestroy(rhoX,ierr)
        call SlepcFinalize(ierr)

    end subroutine entanglement

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Graph Isomorphism ~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function number_of_edges(adjArray, n)
        PetscInt, intent(in)    :: n, adjArray(n,n)
        PetscInt                :: number_of_edges

        ! local variables
        PetscInt                :: i, j

        number_of_edges = 0

        do i=1,n-1
            do j=i+1,n
                if (adjArray(i,j) == 1) number_of_edges = number_of_edges + 1
            enddo
        enddo

    end function number_of_edges

    function getEdgeState(edgeNum,adjArray,n)
        PetscInt, intent(in)    :: n, edgeNum, adjArray(n,n)
        PetscScalar             :: getEdgeState(2,3)

        ! local variables
        PetscInt                :: i, j, number_of_edges

        number_of_edges = 0

        do i=1,n-1
            do j=i+1,n
                if (adjArray(i,j) == 1) number_of_edges = number_of_edges + 1

                if (number_of_edges == edgeNum) then
                    getEdgeState(1,:) = [real(i-1,8),real(i-1,8),1.d0/sqrt(2.d0)]
                    getEdgeState(2,:) = [real(j-1,8),real(j-1,8),1.d0/sqrt(2.d0)]
                    return
                endif
            enddo
        enddo
    end function getEdgeState

    subroutine getAllEdgeStates(init_states,localStateNum,adjArray,n)
        PetscInt, intent(in)    :: n, adjArray(n,n)
        PetscInt, intent(out)   :: localStateNum
        PetscScalar, intent(out):: init_states(n*(n-1)/2,2,3)

        ! local variables
        PetscMPIInt    :: rank, size
        PetscInt       :: edgeNum, i ,j
        PetscErrorCode :: ierr
        
        edgeNum = number_of_edges(adjArray,n)

        call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
        call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr)

        if (edgeNum > size) then
            localStateNum = edgeNum/size

            if (rank < mod(edgeNum,size)) then
                localStateNum = localStateNum + 1
                do i=1, localStateNum
                    init_states(i,:,:) = getEdgeState(localStateNum*rank+i,adjArray,n)
                end do
            else
                do i=1, localStateNum
                    init_states(i,:,:) = getEdgeState(localStateNum*rank+mod(edgeNum,size)+i,adjArray,n)
                end do
            endif

        else
            if (rank < edgeNum) then
                localStateNum = 1
                init_states(1,:,:) = getEdgeState(rank+1,adjArray,n)
            else
                localStateNum = 0
            endif
        endif
    end subroutine getAllEdgeStates

    subroutine getAllEdgeStates3P(init_states,localStateNum,adjArray,n)
        PetscInt, intent(in)    :: n, adjArray(n,n)
        PetscInt, intent(out)   :: localStateNum
        PetscScalar, intent(out):: init_states(n*(n-1)/2,2,4)

        ! local variables
        PetscMPIInt    :: rank, size
        PetscInt       :: edgeNum, i ,j
        PetscReal      :: edgeStates(2,3), vertex(2)
        PetscErrorCode :: ierr
        
        edgeNum = number_of_edges(adjArray,n)

        call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
        call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr)

        if (edgeNum > size) then
            localStateNum = edgeNum/size

            if (rank < mod(edgeNum,size)) then
                localStateNum = localStateNum + 1
                do i=1, localStateNum
                    edgeStates = getEdgeState(localStateNum*rank+i,adjArray,n)
                    vertex = edgeStates(:,1)
                    init_states(i,1,:) = [vertex(1),vertex(1),vertex(1),1.d0/sqrt(2.d0)]
                    init_states(i,2,:) = [vertex(2),vertex(2),vertex(2),1.d0/sqrt(2.d0)]
                end do
            else
                do i=1, localStateNum
                    edgeStates = getEdgeState(localStateNum*rank+mod(edgeNum,size)+i,adjArray,n)
                    vertex = edgeStates(:,1)
                    init_states(i,1,:) = [vertex(1),vertex(1),vertex(1),1.d0/sqrt(2.d0)]
                    init_states(i,2,:) = [vertex(2),vertex(2),vertex(2),1.d0/sqrt(2.d0)]
                end do
            endif

        else
            if (rank < edgeNum) then
                localStateNum = 1
                edgeStates = getEdgeState(rank+1,adjArray,n)
                vertex = edgeStates(:,1)
                init_states(1,1,:) = [vertex(1),vertex(1),vertex(1),1.d0/sqrt(2.d0)]
                init_states(1,2,:) = [vertex(2),vertex(2),vertex(2),1.d0/sqrt(2.d0)]
            else
                localStateNum = 0
            endif
        endif
    end subroutine getAllEdgeStates3P

    subroutine GraphISCert(cert,certLength,adjArray,p,tol,expm_method, & 
        eig_solver,Emax_estimate,worktype,worktypeInt,tolIn,max_it,verbose,N)

        PetscInt, intent(in)         :: N, adjArray(N,N), worktypeInt, max_it, p
        PetscReal, intent(in)        :: tolIn, Emax_estimate, tol
        character(len=*),intent(in)  :: expm_method, eig_solver
        character(len=3), intent(in) :: worktype
        PetscBool, intent(in)        :: verbose

        PetscReal, intent(out)       :: cert(2,(N**(p+1))*(N-1)/2)
        PetscInt, intent(out)        :: certLength

        ! local variables
        PetscMPIInt              :: rank, comm_size
        PetscErrorCode           :: ierr
        Mat                      :: H
        Vec                      :: psi0, psi, temp, cert0
        VecScatter               :: ctx
        PetscInt                 :: i, j, NN, localStateNum, certSize, &
                                    Istart, Iend, counter, certPos
        PetscScalar              :: Emin, Emax, t
        PetscReal                :: Emax_error, currentProb
        PetscReal, allocatable   :: certArrayReal(:)
        PetscScalar, allocatable :: localCert(:), init_states(:,:,:)
        PetscScalar, pointer     :: certArray(:)
        character(len=12)        :: expArg, eigArg, wtArg

        call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
        call MPI_Comm_size(PETSC_COMM_WORLD,comm_size,ierr)

        ! set time
        t = 2*N

        ! create local Hamiltonian
        call MatCreate(MPI_COMM_SELF,H,ierr)
        call MatSetType(H, MATSEQAIJ,ierr)
        if (p==3) then
            call adjToH(H,adjArray,'3',[0],[0.d0*PETSC_i],1.d0+0.d0*PETSC_i,1,N)
        else
            call adjToH(H,adjArray,'2',[0],[0.d0*PETSC_i],1.d0+0.d0*PETSC_i,1,N)
        endif

        NN = N**p

        allocate(init_states(N*(N-1)/2,2,p+1))
        ! get local initial states
        if (p==3) then
            call getAllEdgeStates3P(init_states,localStateNum,adjArray,N)
        else
            call getAllEdgeStates(init_states,localStateNum,adjArray,N)
        endif
        allocate(localCert(NN*localStateNum))

        call VecCreateSeq(MPI_COMM_SELF,NN,psi0,ierr)
        call VecCreateSeq(MPI_COMM_SELF,NN,psi,ierr)

        ! if using the Chebyshev method, calculate the eigenvalues
        expArg = trim(adjustl(expm_method))
        if (expArg == 'chebyshev') then
            if (Emax_estimate == 0.) then
                eigArg = trim(adjustl(eig_solver))
                wtArg = trim(adjustl(worktype))
                call min_max_eigs(H,rank,Emax,Emax_error,'max',&
                    eigArg,wtArg,worktypeInt,tolIn,max_it,verbose,ierr)
            else
                Emax = Emax_estimate
            endif

            Emin = 0.d0
        endif

        ! loop over all local initial states
        do i=1, localStateNum

            if (verbose) write(*,*)'rank = ',rank,' Calculating init state = ', &
                int(init_states(i,1,1)), int(init_states(i,2,1))

            ! create the initial state vector
            if (p==3) then 
                call p3_init(psi0,init_states(i,:,:),2,N)
            else
                init_states(i,:,1:2) = init_states(i,:,1:2)-N/2+1
                call p2_init(psi0,init_states(i,:,:),2,N)
            endif

            ! propagate the CTQW
            if (expArg == 'chebyshev') then
                call qw_cheby(psi0,psi,t,H,Emin,Emax,rank,N)
            else
                call qw_krylov(H,t,psi0,psi)
            endif

            ! find the sqrt(|psi|) of the vector
             call VecSqrtAbs(psi,ierr)

            ! store the probabilities in a local certificate array
            call VecGetValues(psi,NN,[(j,j=0,NN-1)],localCert((i-1)*NN+1:i*NN),ierr)

            ! reset initial state vector
            call VecSet(psi0,0.d0*PETSc_i,ierr)
        enddo

        if (comm_size>1) then
            ! create a global vector
            call VecCreateMPI(PETSC_COMM_WORLD,NN*localStateNum,PETSC_DECIDE,temp,ierr)
            call VecGetOwnershipRange(temp,Istart,Iend,ierr)

            ! each node places it's local certificate in the global vector
            call VecSetValues(temp,NN*localStateNum,[(j,j=Istart,Iend-1)],localCert,INSERT_VALUES,ierr) 

            ! get total size of certificate, and create
            ! a sequential vector of this length on node 0
            if (rank == 0) then
                call VecGetSize(temp,certSize,ierr)
            else
                certSize = 0
            endif
            
            ! create vector on node 0
            call VecCreate(MPI_COMM_SELF,cert0,ierr)
            call VecSetType(cert0,VECSEQ,ierr)
            call VecSetSizes(cert0,certSize,PETSC_DECIDE,ierr)

            ! Vector scatter the global certificate to node 0
            call VecScatterCreateToZero(temp,ctx,cert0,ierr)
            call VecScatterBegin(ctx,temp,cert0,INSERT_VALUES,SCATTER_FORWARD,ierr)
            call VecScatterEnd(ctx,temp,cert0,INSERT_VALUES,SCATTER_FORWARD,ierr)
            call VecScatterDestroy(ctx,ierr)

            call VecDestroy(temp,ierr)
        endif

        ! clean up
        call MatDestroy(H,ierr)
        call VecDestroy(psi0,ierr)
        call VecDestroy(psi,ierr)
        deallocate(init_states)

        !call VecView(cert0,PETSC_VIEWER_STDOUT_WORLD,ierr)
        if (rank==0) then 
            if (comm_size>1) then
                call VecGetArrayF90(cert0,certArray,ierr)
                
                allocate(certArrayReal(certSize))
                certArrayReal = certArray
            else
                certSize = NN*localStateNum
                allocate(certArrayReal(certSize))
                certArrayReal = localCert
            endif


            call d_refsor(certArrayReal)

            currentProb = certArrayReal(1)
            counter = 1
            certPos = 1

            do i=1, certSize
                if (certArrayReal(i) - currentProb > tol) then
                    cert(:,certPos) = [currentProb,real(counter,8)]
                    currentProb = certArrayReal(i)
                    counter = 1
                    certPos = certPos + 1
                else
                    counter = counter + 1
                endif
            enddo

            deallocate(certArrayReal,localCert)

            certLength = certPos-1

            if (comm_size>1) then
                call VecRestoreArrayF90(cert0,certArray,ierr)
                call VecDestroy(cert0,ierr)
            endif

        endif

    end subroutine GraphISCert


end module ctqwMPI
