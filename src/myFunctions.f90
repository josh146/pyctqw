module myFunctions
	use IFPORT
	implicit none
	complex(8), parameter	::	ii = (0.d0,1.d0)
	real(8), parameter		::	pi = 4.d0*atan(1.d0)

	interface
		subroutine zgeev(JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, &
		& LWORK, RWORK, INFO)
			character 	:: JOBVL, JOBVR
			integer 	:: INFO, LDA, LDVL, LDVR, LWORK, N
			real(8) 	:: RWORK(*)
			complex(8) 	:: A(*), VL(*), VR(*), W(*), WORK(*)
		end subroutine zgeev
	end interface

	contains

	function coord(x,y,n)
		integer, intent(in)	:: n, x, y
		integer			:: coord
		
		coord = n*(x + n/2 - 1) + y + n/2
	end function coord
	
	function marginal(psi,p)
		complex(8), intent(in)	:: psi(:)
		character, intent(in)	:: p

		! local variables
		integer			:: i, j, N
		complex(8), allocatable	:: marginal(:)

		N = sqrt(real(size(psi)))
		allocate(marginal(N))
		
		if (p=='x') then
			do i=1,N
				marginal(i) = sum(abs(psi(1+(i-1)*N:i*N))**2.d0)
			end do
		elseif (p=='y') then
			do i=1,N
				marginal(i) = sum(abs([(psi(i+j*N), j=0, N-1)])**2.d0)
			end do	
		endif		
	end function marginal

	function identity(n)
		integer, intent(in)	::	n
		integer			::	i
		real(8)		::	identity(n,n)

		identity = 0.d0
		do i = 1, n
			identity(i,i)=1.d0
		end do
	end function identity
	
	function kron(M1, M2)
		real(8), intent(in)	:: M1(:,:), M2(:,:)
		real(8), allocatable	:: kron(:,:)

		! local variables
		integer :: i, j, dim1(2), dim2(2), r1, c1, r2, c2
		dim1 = shape(M1); r1=dim1(1); c1=dim1(2)
		dim2 = shape(M2); r2=dim2(1); c2=dim2(2)
		
		allocate(kron(r1*r2,c1*c2))
		kron = (0.d0,0.d0)

		forall (i=1:r1, j=1:c1)
			kron(r2*(i-1)+1:r2*i, c2*(j-1)+1:c2*j) = M1(i,j)*M2
		end forall
	end function kron

	subroutine extremeEv(H,Emin,Emax)
		real(8), intent(in)	::	H(:,:)
		real(8), intent(out)	::	Emin, Emax

		! local variables
		integer			::	N, info, i, j
		real(8), allocatable	::	rwork(:)
		complex(8), allocatable	::	UpperArray(:,:), leftvectors(:,:), &
						& array(:,:), rightvectors(:,:), eigenvalues(:), work(:)

		N = size(H,1)

		! allocate the arrays required for LAPACK
		allocate(eigenvalues(N), work(2*N), rwork(2*N), leftvectors(N, N), rightvectors(N, N),array(N,N))
		array = H

		! find the eigenvalues of H
		call zgeev('N', 'N', N, array, N, eigenvalues, leftvectors, &
				& N, rightvectors, N, work, 2*N, rwork, info)

		Emin = minval(real(eigenvalues))
		Emax = maxval(real(eigenvalues))
		deallocate(rwork,work,eigenvalues,leftvectors,rightvectors)
	end subroutine extremeEv

	function quantumWalk(psi,dt,H,writecoeff)
		complex(8), intent(in)	::	psi(:)
		real(8), intent(in)	::	dt, H(:,:)
		logical, intent(in)	::	writecoeff

		! local variables
		integer			::	N, m, terms, i, j
		real(8)			::	alpha, Emax, Emin
		complex(8), allocatable	::	phi0(:), phi1(:), phi2(:), U(:), quantumWalk(:)

		N = size(psi)
		allocate(phi0(N),phi1(N),phi2(N),U(N),quantumWalk(N))

		! set Chebyshev variables
		call extremeEv(H,Emin,Emax)
		alpha = (Emax-Emin)*dt/2.d0

		phi0 = 1.d0*psi
		phi1 = -(2.d0*matmul(H,psi)-(Emax+Emin)*psi)/(Emax-Emin)
		U = dbesjn(0,alpha)*phi0 + 2.d0*ii*dbesjn(1,alpha)*phi1

		terms = 0
		do while (abs(2.d0*dbesjn(terms,alpha)) > 1.d-45)
			terms = terms + 1
		end do

		open(26,file="coeff.txt",status='replace')
		do m = 2, terms
			call progressbar(m,terms)
			phi2 = -2.d0*(2.d0*matmul(H,phi1)-(Emax+Emin)*phi1)/(Emax-Emin) - phi0
			U = U + 2.d0*(ii**m)*dbesjn(m,alpha)*phi2

			if (writecoeff) write(26,"(E24.15E3)")abs(2.d0*dbesjn(m,alpha))

			phi0 = phi1
			phi1 = phi2
		end do
		close(26)
		
		write(*,*)"\n"
		quantumWalk = exp(-ii*(Emax+Emin)*dt/2.d0)*U
		deallocate(phi0,phi1,phi2,U)
	end function quantumWalk

	! finds the matrix exponential exp(-iHt)
	function matrixExp(H,t)
		real(8), intent(in)	::	H(:,:)
		real(8), intent(in)	::	t

		! local variables
		integer			::	N, m, terms, i, j
		real(8)			::	alpha, Emax, Emin
		complex(8), allocatable	::	phi0(:,:), phi1(:,:), phi2(:,:), U(:,:), matrixExp(:,:)

		N = size(H,1)
		allocate(phi0(N,N),phi1(N,N),phi2(N,N),U(N,N),matrixExp(N,N))

		! set Chebyshev variables
		call extremeEv(H,Emin,Emax)
		alpha = (Emax-Emin)*t/2.d0

		phi0 = identity(N)
		phi1 = -(2.d0*H-(Emax+Emin)*identity(N))/(Emax-Emin)
		U = dbesjn(0,alpha)*phi0 + 2.d0*ii*dbesjn(1,alpha)*phi1

		terms = 0
		do while (abs(2.d0*dbesjn(terms,alpha)) > 1.d-45)
			terms = terms + 1
		end do

		do m = 2, terms
			phi2 = -2.d0*(2.d0*matmul(H,phi1)-(Emax+Emin)*phi1)/(Emax-Emin) - phi0
			U = U + 2.d0*(ii**m)*dbesjn(m,alpha)*phi2

			phi0 = phi1
			phi1 = phi2
		end do

		matrixExp = exp(-ii*(Emax+Emin)*t/2.d0)*U
		deallocate(phi0,phi1,phi2,U)
	end function matrixExp
	
	subroutine quantumWalk2(psi0,psi,dt,H)
		complex(8), intent(in)	:: psi0(:)
		real(8), intent(in)	:: dt, H(:,:)
		complex(8), intent(out)	:: psi(:)

		! local variables
		integer			:: N
		complex(8), allocatable	:: U(:,:)

		N = size(psi0)
		allocate(U(N,N))
		
		call c8mat_expm1(N,-ii*H*dt,U)
		psi = matmul(U,psi0)
		
		deallocate(U)
	end subroutine quantumWalk2

	function ToKspace(psi,x)
		complex(8), intent(in)	::	psi(:)
		real(8), intent(in)	::	x(:)
		complex(8), allocatable	::	toKspace(:)

		! local variables
		integer			::	j, N
		real(8), allocatable	::	kgrid(:)
				
		N = size(psi)
		allocate(toKspace(N),kgrid(N))

		kgrid = [(-pi+2.d0*pi*j/N, j=0, N-1)]
		toKspace = 0.d0
		do j=1, N
			toKspace(j) = sum(psi*exp(-ii*kgrid(j)*x)) / sqrt(real(N,8))
		end do

		deallocate(kgrid)
	end function toKspace

	function tCoeff(psi0,psiT,x)
		complex(8), intent(in)	::	psi0(:), psiT(:)
		real(8), intent(in)	::	x(:)

		! local variables
		integer			::	j, m, N
		integer, allocatable	::	T0pnts(:)
		real(8)			::	kpsi0Max
		real(8), allocatable	::	tCoeff(:,:), kgrid(:)
		complex(8), allocatable	::	kpsi0(:), kpsiT(:)

		N = size(psi0)
		allocate(kpsi0(N),kpsiT(N),kgrid(N),T0pnts(N))

		! Transform to k-space
		kpsi0 = ToKspace(psi0,x)
		kpsiT = ToKspace(psiT,x)
		kgrid = [(-pi+2.d0*pi*j/N, j=0, N-1)]

		! find the max values of kpsi0
		m = 0
		kpsi0Max = maxval(abs(kpsi0)**2.d0)
		do j = 1, N
			if (abs(kpsi0(j))**2.d0 .ge. 0.001d0*kpsi0Max) then
				m = m+1
				T0pnts(m) = j
			end if
		end do

		if (m>0) then
			allocate(tCoeff(2,m))
			do j = 1, m
				tCoeff(:,j) = [kgrid(T0pnts(j)), &
						& abs(kpsiT(T0pnts(j))/kpsi0(T0pnts(j)))**2.d0]
			end do
		endif

		deallocate(kpsi0,kpsiT)
	end function tCoeff

	subroutine progressbar(i,NN)
		integer, intent(in)	::	i, NN
		integer			::	k
		character(len=39)	::	bar=" \r???%|                              |"

		! updates the fraction of calculation done
		write(unit=bar(3:5),fmt="(i3)") (100*i)/NN
		do k = 1, (i*30)/NN
			bar(7+k:7+k)="="
		enddo

		open (unit=6, carriagecontrol='fortran')
		write(6,'(3a)',advance='no')'+',CHAR(13),bar
		flush 6

		! once the progress bar is full, clear it for further use
		if ((100*i)/NN == 100) bar=" \r???%|                              |"
	end subroutine progressbar

end module myFunctions
