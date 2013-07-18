module ctqw
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
	
	interface
		subroutine dsbevd(JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, LWORK, &
		& IWORK, LIWORK, INFO)
			integer		:: N, KD, LDAB, LDZ, LWORK, IWORK(*), LIWORK, INFO
			real(8)		:: AB(LDAB,*), W(*), Z(LDZ,*), WORK(*)
			character	:: JOBZ, UPLO
		end subroutine dsbevd
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
	
	subroutine pymarginal(psi,NN,p,psip,N)
		character, intent(in)	:: p
		integer, intent(in)	:: N, NN
		complex(8), intent(in)	:: psi(NN)
		complex(8), intent(out)	:: psip(N)
		
		psip = marginal(psi,p)	
	end subroutine pymarginal

	subroutine identity(mat,n)
		integer, intent(in) :: n
		real(8), intent(out) :: mat(n,n)
		
		integer :: i

		mat = 0.d0
		do i = 1, n
			mat(i,i)=1.d0
		end do
	end subroutine identity
	
	subroutine kron(M1,r1,c1, M2,r2,c2, kronProd)
		integer, intent(in)	:: r1, c1, r2, c2
		real(8), intent(in)	:: M1(r1,c1), M2(r2,c2)
		real(8), intent(out)	:: kronProd(r1*r2,c1*c2)

		! local variables
		integer :: i, j
		
		kronProd = 0.d0
		forall (i=1:r1, j=1:c1)
			kronProd(r2*(i-1)+1:r2*i, c2*(j-1)+1:c2*j) = M1(i,j)*M2
		end forall
	end subroutine kron
	
	subroutine hamiltonian_1p(H,d,a,nd,N)
		integer, intent(in)	:: nd, N
		real(8), intent(in)	:: a(nd)
		integer, intent(in)	:: d(nd)
		real(8), intent(out)	:: H(N,N)
		
		! local variables
		integer	:: i, j, k
		
		H = 0.d0
		do j = 1, N
			H(j,j) = 2.d0
			do i=1, nd
				H(int(d(i)+N/2.d0),int(d(i)+N/2.d0)) = 2.d0 + a(i)
			end do

			do k = 1, N
				if (abs(k-j)==1) H(j,k) = -1.d0
			end do
		end do
	end subroutine hamiltonian_1p
	
	subroutine hamiltonian_2p_noint(H1,H2,N)
		integer, intent(in)	:: N
		real(8), intent(in)	:: H1(N,N)
		real(8), intent(out)	:: H2(N*N,N*N)
		
		real(8), allocatable	:: ident(:,:), temp(:,:)
		
		allocate(ident(N,N), temp(N**2,N**2))
		
		call identity(ident,N)
		call kron(H1,N,N,ident,N,N,H2)
		call kron(ident,N,N,H1,N,N,temp)

		H2 = H2 + temp
		deallocate(ident, temp)
	end subroutine hamiltonian_2p_noint
	
	subroutine sym_band_storage(A,AB,UPLO,kd,N)
		integer, intent(in)		:: N, kd
		character(len=1), intent(in)	:: UPLO
		real(8), intent(in)		:: A(N,N)
		real(8), intent(out)		:: AB(kd+1,N)
		
		integer	:: i, j
		
		AB = 0.d0
		select case(UPLO)
			case('U')
				do j=1, N
					AB(kd+1+max(1,j-kd)-j:kd+1,j) = A(max(1,j-kd):j,j)
					!(AB(kd+1+i-j,j) = A(i,j), i=max(1,j-kd), j)
				end do
			case('L')
				do j=1, N
					AB(1:1+min(N,j+kd)-j,j) = A(j:min(N,j+kd),j)
					!(AB(1+i-j,j) = A(i,j), i=j, min(N,j+kd))
				end do
		end select
	end subroutine sym_band_storage

	subroutine sband_extremeEv(H,Emin,Emax)
		real(8), intent(in)	::	H(:,:)
		real(8), intent(out)	::	Emin, Emax

		! local variables
		integer			::	N, i, j, kd, IWORK(1), INFO
		real(8)			::	Z(1,1)
		real(8), allocatable	::	eigenvalues(:), AB(:,:), WORK(:)

		N = size(H,1)
		kd = sqrt(real(N))

		! allocate the arrays required for LAPACK
		allocate(AB(kd+1,N), eigenvalues(N), WORK(2*N))
				
		call sym_band_storage(H,AB,'L',kd,N)

		! find the eigenvalues of H
		call dsbevd('N', 'L', N, kd, AB, 1+kd, eigenvalues, Z, 1, WORK, 2*N, IWORK, 1, INFO)
		!call dsbevd(AB, eigenvalues, 'L')

		Emin = minval(real(eigenvalues))
		Emax = maxval(real(eigenvalues))
		deallocate(AB, eigenvalues)
	end subroutine sband_extremeEv
	
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

	subroutine qw_cheby(psi,psiT,dt,H,Emin,Emax,N)
		integer, intent(in)	:: N
		complex(8), intent(in)	:: psi(N)
		real(8), intent(in)	:: dt, H(N,N), Emax, Emin
		complex(8), intent(out)	:: psiT(N)

		! local variables
		integer			:: m, terms, i, j
		real(8)			:: alpha
		complex(8),allocatable	:: phi0(:), phi1(:), phi2(:), U(:)

		allocate(phi0(N),phi1(N),phi2(N),U(N))

		! set Chebyshev variables
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

			write(26,"(E24.15E3)")abs(2.d0*dbesjn(m,alpha))

			phi0 = phi1
			phi1 = phi2
		end do
		close(26)
		
		psiT = exp(-ii*(Emax+Emin)*dt/2.d0)*U
		deallocate(phi0,phi1,phi2,U)
	end subroutine qw_cheby

	! finds the matrix exponential exp(Ht) using Chebyshev
	function matrixExp_cheby(H,t)
		real(8), intent(in)	::	H(:,:)
		real(8), intent(in)	::	t

		! local variables
		integer			::	N, m, terms, i, j
		real(8)			::	alpha, Emax, Emin
		real(8), allocatable	::	ident(:,:)
		complex(8), allocatable	::	phi0(:,:), phi1(:,:), phi2(:,:), &
						 & U(:,:), matrixExp_cheby(:,:)

		N = size(H,1)
		allocate(phi0(N,N),phi1(N,N),phi2(N,N),U(N,N),matrixExp_cheby(N,N),ident(N,N))

		! set Chebyshev variables
		call extremeEv(H,Emin,Emax)
		alpha = (Emax-Emin)*t/2.d0

		call identity(ident,N)

		phi0 = ident
		phi1 = (2.d0*H-(Emax+Emin)*ident)/(Emax-Emin)
		U = dbesjn(0,alpha)*phi0 + 2.d0*dbesjn(1,alpha)*phi1

		terms = 0
		do while (abs(2.d0*dbesjn(terms,alpha)) > 1.d-45)
			terms = terms + 1
		end do

		do m = 2, terms
			phi2 = 2.d0*(2.d0*matmul(H,phi1)-(Emax+Emin)*phi1)/(Emax-Emin) + phi0
			U = U + 2.d0*dbesjn(m,alpha)*phi2

			phi0 = phi1
			phi1 = phi2
		end do

		matrixExp_cheby = exp((Emax+Emin)*t/2.d0)*U
		deallocate(phi0,phi1,phi2,U,ident)
	end function matrixExp_cheby
	
	subroutine pyexpm_cheby(H,t,expHt,N)
		integer, intent(in)	:: N
		real(8), intent(in)	:: H(N,N), t
		complex(8), intent(out)	:: expHt(N,N)
		
		expHt = matrixExp_cheby(H,t)		
	end subroutine pyexpm_cheby
	
	subroutine qw_Burkadt(psi0,psi,dt,H,N)
		integer, intent(in)	:: N
		complex(8), intent(in)	:: psi0(N)
		real(8), intent(in)	:: dt, H(N,N)
		complex(8), intent(out)	:: psi(N)
		
		complex(8),allocatable	:: U(:,:)
		
		allocate(U(N,N))
		
		call c8mat_expm1(N,-ii*H*dt,U)
		psi = matmul(U,psi0)
		
		deallocate(U)
	end subroutine qw_Burkadt

	subroutine progressbar(i,NN)
		integer, intent(in)	::	i, NN
		integer			::	k
		character(len=39)	::	bar="  \r?% |                              |"

		! updates the fraction of calculation done
		write(unit=bar(3:5),fmt="(i3)") (100*i)/NN
		do k = 2, (i*30)/NN+1
			bar(7+k:7+k)="="
		enddo

		open (unit=6, carriagecontrol='fortran')
		write(6,'(3a)',advance='no')' ',CHAR(13),bar
		flush 6

		! once the progress bar is full, clear it for further use
		if ((100*i)/NN == 100) bar="  \r?% |                              |"
	end subroutine progressbar

end module ctqw
