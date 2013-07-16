program CtQW
	use iso_fortran_env
	use myFunctions
	use fileOps
	implicit none

	integer			::	narg, i, j, k, N, steps, d1, d2
	real(8)			::	t, a1, a2
	real(8), allocatable	:: 	H(:,:), H2(:,:)
	complex(8), allocatable	::	psi(:),psi0(:)
	logical			::	writecoeff, burkadt=.false.
	character(len=32)	::	arg

	! set the variables
	N = 20
	t = 1.d0;
	a1 = 1.d0;	d1 = 3
	a2 = 1.d0;	d2 = 4

	!Check if arguments are found
	narg = command_argument_count()
	
	if(narg>0)then
		do i=1, narg
			call get_command_argument(i,arg)
			select case(adjustl(arg))
				case("--burkadt"); burkadt=.true.
			end select
		end do
	end if

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
	call hamiltonian_1p(H,d1,a1,d2,a2,N)
	call hamiltonian_2p_noint(H,H2,N)
	
	select case(burkadt)
		case(.true.); call quantumWalk_Burkadt(psi0,psi,t,H2)
		case default; psi = quantumWalk(psi0,t,H2,writecoeff)
	end select
	
	call writemarginal(marginal(psi,'x'),1,'x')
	call writemarginal(marginal(psi,'y'),1,'y')

	deallocate(psi,psi0,H,H2)
end program CtQW
