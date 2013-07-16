#!/usr/bin/make -f

ifeq ($(FC),intel)
	compiler = ifort
	libName = "_intel"
else
	compiler = gfortran
	libName = "_gcc"
endif

switch = -O2

#Makefile
all: libmatExp.a libc8lib.a libr8lib.a cleanTemp

libmatExp.a: temp_mat/c8mat_expm1.f90
	cd temp_mat; for FILE in $(shell cd temp_mat && ls *.f90); \
		do $(compiler) -c -g $(switch) $$FILE; \
	done
	ar qc libmatExp.a temp_mat/*.o
	mkdir -p ../../lib
	mv libmatExp.a ../../lib/libmatExp$(libName).a
	
temp_mat/c8mat_expm1.f90: matrix_exponential.f90 f90split
	mkdir -p temp_mat; cd temp_mat; \
	../f90split ../matrix_exponential.f90
	
libr8lib.a: temp_r8/gamma_values.f90
	cd temp_r8; for FILE in $(shell cd temp_r8 && ls *.f90); \
		do $(compiler) -c -g $(switch) $$FILE; \
	done
	ar qc libr8lib.a temp_r8/*.o
	mv libr8lib.a ../../lib/libr8$(libName).a
	
temp_r8/gamma_values.f90: r8lib.f90 f90split
	mkdir -p temp_r8; cd temp_r8; \
	../f90split ../r8lib.f90
	
libc8lib.a: temp_c8/c8_abs.f90
	cd temp_c8; for FILE in $(shell cd temp_c8 && ls *.f90); \
		do $(compiler) -c -g $(switch) $$FILE; \
	done
	ar qc libc8lib.a temp_c8/*.o
	mv libc8lib.a ../../lib/libc8$(libName).a

temp_c8/c8_abs.f90: c8lib.f90 f90split
	mkdir -p temp_c8; cd temp_c8; \
	../f90split ../c8lib.f90
	
f90split: f90split.f90
	$(compiler) $(switch) f90split.f90 -o f90split
	
#Cleaning files
cleanTemp:
	rm -rf temp*
	
clean:
	rm f90split
	rm -rf temp*



