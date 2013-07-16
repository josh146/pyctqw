module fileOps
	implicit none
	contains

	subroutine writestate(psi,step,res)
		integer, intent(in)	::	step
		complex(8), intent(in)	::	psi(:)
		character,intent(in)	::	res

		! local variables
		character(len=50)	:: level1, filename
		integer			:: m, N, minN, maxN

		N = size(psi)
		minN = floor(1.d0-N/2.d0);	maxN = N/2

		! adjust the filename to include 't'
		write(level1,*)step
		filename = "output/output_"//trim(adjustl(res))//"_t"//trim(adjustl(level1))// ".txt"

		! write the wavefunction to a file
		open(12,file=filename)
		write(12,"(1a10,3a22)")"node","Re(psi)","Im(psi)","psi^2"
		do m = 1, N
			write(12,15)minN+m-1, psi(m), abs(psi(m))**2.d0
		end do
		15 format(i5,3E24.15E3)
		close(12)

!		write(*,*)"File "//trim(filename)//" created."
	end subroutine writestate

	subroutine writemarginal(psi,step,res)
		integer, intent(in)	::	step
		complex(8), intent(in)	::	psi(:)
		character,intent(in)	::	res

		! local variables
		character(len=50)	:: level1, filename
		integer			:: m, N, minN, maxN

		N = size(psi)
		minN = floor(1.d0-N/2.d0);	maxN = N/2

		! adjust the filename to include 't'
		write(level1,*)step
		filename = "output/output_"//trim(adjustl(res))//"_t"//trim(adjustl(level1))// ".txt"

		! write the wavefunction to a file
		open(12,file=filename)
		write(12,"(1a10,1a22)")"node","psi^2"
		do m = 1, N
			write(12,15)minN+m-1, psi(m)
		end do
		15 format(i5,3E24.15E3)
		close(12)

!		write(*,*)"File "//trim(filename)//" created."
	end subroutine writemarginal

end module fileOps
