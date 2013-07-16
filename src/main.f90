program CtQW
	use iso_fortran_env
	use myFunctions
	use fileOps
	implicit none

	integer			::	i, j, k, N, steps, xv01, xv02
	real(8)			::	t, V01, V02
	real(8), allocatable	:: 	H(:,:), H2(:,:)
	complex(8), allocatable	::	psi(:),psi0(:)
	logical			::	writecoeff
	character(len=32)	::	arg

	! set the variables
	N = 20
	t = 1.d0;
	V01 = 1.d0;	xv01 = 3
	V02 = 1.d0;	xv02 = 4

!	! get the first argument provided to the program, and store it in xv02
!	call getarg(1, arg)
!	read(arg,*,iostat=i)width
!	xv02 = xv02 + width

!	if (i==5010) then
!		! if the argument is text, provide syntax help and exit
!		write(*,*)"Usage: CtQW [POTENTIAL WIDTH]"
!		stop
!	elseif (i==-1 .OR. i>0) then
!		write(*,*)"ERROR: No argument provided"
!		write(*,*)"Usage: CtQW [POTENTIAL WIDTH]"
!		stop
!	endif

	! allocate the vectors
	allocate(H(N,N),H2(N**2,N**2),psi(N**2),psi0(N**2))
	
	! set the initial state
	psi0 = 0.d0
	psi0(coord(0,0,N)) = 1.d0/sqrt(2.d0)
	psi0(coord(1,1,N)) = 1.d0/sqrt(2.d0)
	
	call writestate(marginal(psi0,'x'),0,'x')
	call writestate(marginal(psi0,'y'),0,'y')

	write(*,*)"Propagating the continuous time quantum walk..."
	writecoeff = .false.
	
	! create the Hamiltonian matrix
	H = 0.d0
	do j = 1, N
		if (floor(j-N/2.d0)==real(xv01,8)) then
			H(j,j) = 2.d0 + V01
		elseif (floor(j-N/2.d0)==real(xv02,8) .and. floor(j-N/2.d0)/=real(xv01,8)) then
			H(j,j) = 2.d0 + V02
		else
			H(j,j) = 2.d0
		end if

		do k = 1, N
			if (abs(k-j)==1) H(j,k) = -1.d0
		end do
	end do
	
	H2 = kron(H,identity(N)) + kron(identity(N),H)
	call quantumWalk2(psi0,psi,t,H2)
	!psi = quantumWalk(psi0,t,H2,writecoeff)
	
	call writemarginal(marginal(psi,'x'),1,'x')
	call writemarginal(marginal(psi,'y'),1,'y')

	deallocate(psi,psi0,H,H2)
end program CtQW
