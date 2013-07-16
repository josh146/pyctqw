module FFT
	implicit none
	real(8), parameter		::	PI = 4.d0*atan(1.d0)
	complex(8), parameter	::	ii = (0.d0,1.d0)
	include "fftw3.f"
	integer(8)					::	plan

	contains

	function CFFT(f)
		complex(8), intent(in)	::	f(:)
		complex(8), allocatable	::	temp(:), temp2(:), CFFT(:)

		allocate(temp(size(f)),temp2(size(f)),CFFT(size(f)))

		! perform the DFT
		temp = cshift(f,size(f)/2)
		call dfftw_plan_dft_1d(plan,size(f),temp,temp2,FFTW_FORWARD,FFTW_ESTIMATE)
		call dfftw_execute_dft(plan, temp, temp2)
		call dfftw_destroy_plan(plan)

		CFFT = cshift(temp2,size(f)/2)
		deallocate(temp,temp2)
	end function CFFT

	function ICFFT(f)
		complex(8), intent(in)	::	f(:)
		complex(8), allocatable	::	temp(:), temp2(:), ICFFT(:)

		allocate(temp(size(f)),temp2(size(f)),ICFFT(size(f)))

		! perform the DFT
		temp = cshift(f,size(f)/2)
		call dfftw_plan_dft_1d(plan,size(f),temp,temp2,FFTW_BACKWARD,FFTW_ESTIMATE)
		call dfftw_execute_dft(plan, temp, temp2)
		call dfftw_destroy_plan(plan)

		ICFFT = cshift(temp2,size(f)/2)/size(f)
		deallocate(temp,temp2)
	end function ICFFT

	! DFT differentiation
	function diff(f,dx,order)
		complex(8), intent(in)	::	f(:)
		real(8), intent(in)		::	dx
		integer, intent(in)		::	order
		! local variables
		integer					::	i, j, N
		complex(8), allocatable	::	freq(:), temp(:), diff(:)

		N = size(f)
		allocate(freq(N),temp(N),diff(N))

		! perform the centered DFT
		temp = CFFT(f)

		! multiply the DFT by 2pi*I*k*dk/N to find
		! the derivative in Fourier space. Note that
		! dk=1/(Ndx), and that the k values repeat
		freq = ii*[(-pi/dx+2.d0*pi*i/(N*dx), i=0, N-1)]
		temp = (freq**order) * temp

		! find the inverse fourier transform
		diff = ICFFT(temp)
	end function diff
end module FFT
