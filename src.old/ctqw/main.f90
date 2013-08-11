program main
	use iso_fortran_env
	use ctqw
	use fileOps
	implicit none

	integer			:: narg, i, j, k, N, steps, d(2)
	real(8)			:: t, a(2), Emin, Emax
	real(8), allocatable	:: H(:,:), H2(:,:)
	complex(8), allocatable	:: psi(:),psi0(:), psiX(:)
	logical			:: burkadt=.false.
	character(len=32)	:: arg
	integer			:: t1,t2,rate

	! set the variables
	N = 50
	t = 1.d0;
	d = [3, 4]
	a = [1.d0, 1.d0]

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
	allocate(H(N,N),H2(N**2,N**2),psi(N**2),psi0(N**2),psiX(N))
	
	! set the initial state
	psi0 = 0.d0
	psi0(coord(0,0,N)) = 1.d0/sqrt(2.d0)
	psi0(coord(1,1,N)) = 1.d0/sqrt(2.d0)
	
	call writestate(marginal(psi0,'x'),0,'x')
	call writestate(marginal(psi0,'y'),0,'y')

	write(*,*)"Propagating the continuous time quantum walk..."
	
	! create the Hamiltonian matrix
	call system_clock(t1, rate)
	call hamiltonian_1p(H,d,a,2,N)
	call hamiltonian_2p_noint(H,H2,N)
	call system_clock(t2, rate)
	write(*,*) "Hamiltonian creation time: ", real(t2 - t1) / real(rate)
	
	call system_clock(t1, rate)
	select case(burkadt)
		case(.true.); call qw_Burkadt(psi0,psi,t,H2,N**2)
		case default
			call sband_extremeEv(H2,Emin,Emax)
			call qw_cheby(psi0,psi,t,H2, Emin, Emax,N**2)
	end select
	call system_clock(t2, rate)
	write(*,*)
	write(*,*) "QW propagation time: ", real(t2 - t1) / real(rate)
	
	call writemarginal(marginal(psi,'x'),1,'x')
	call writemarginal(marginal(psi,'y'),1,'y')

	deallocate(psi,psi0,H,H2,psiX)
end program main
