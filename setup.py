#!/usr/bin/env python2.7
import sys, os
from os.path import join, isdir
	
def MPIpackage(config):
	try:
		from numpy.distutils.fcompiler import FCompiler
		def runtime_library_dir_option(self, dir):
			return self.c_compiler.runtime_library_dir_option(dir)
		FCompiler.runtime_library_dir_option = \
			runtime_library_dir_option
	except Exception:
		pass
	
	INCLUDE_DIRS = []
	LIBRARY_DIRS = []
	LIBRARIES	= []

	# PETSc
	try:
		PETSC_DIR  = os.environ['PETSC_DIR']
	except:
		print "ERROR: PETSC_DIR environment variable not set, pyCTQW.MPI will not be installed"
		return
		
	PETSC_ARCH = os.environ.get('PETSC_ARCH', '')
	
	if PETSC_ARCH and isdir(join(PETSC_DIR, PETSC_ARCH)):
		INCLUDE_DIRS += [join(PETSC_DIR, PETSC_ARCH, 'include'),
						 join(PETSC_DIR, 'include')]
		LIBRARY_DIRS += [join(PETSC_DIR, PETSC_ARCH, 'lib')]
	else:
		if PETSC_ARCH: pass
		INCLUDE_DIRS += [join(PETSC_DIR, 'include')]
		LIBRARY_DIRS += [join(PETSC_DIR, 'lib')]
	LIBRARIES += ['petsc']

	# PETSc for Python
	try:
		import petsc4py
	except:
		print "ERROR: petsc4py not installed, pyCTQW.MPI will not be installed"
		return
		
	INCLUDE_DIRS += [petsc4py.get_include()]
				  
	# SLEPc
	
	try:
		SLEPC_DIR  = os.environ['SLEPC_DIR']
	except:
		print "ERROR: SLEPC_DIR environment variable not set, pyCTQW.MPI will not be installed"
		return
	
	if PETSC_ARCH and isdir(join(PETSC_DIR, PETSC_ARCH)):
		INCLUDE_DIRS += [join(SLEPC_DIR, PETSC_ARCH, 'include'),
						 join(SLEPC_DIR, 'include')]
		LIBRARY_DIRS += [join(SLEPC_DIR, PETSC_ARCH, 'lib')]
	else:
		INCLUDE_DIRS += [join(SLEPC_DIR, 'include')]
		LIBRARY_DIRS += [join(SLEPC_DIR, 'lib')]
	LIBRARIES += ['slepc']

#	# SLEPC for Python
#	import slepc4py
#	INCLUDE_DIRS += [slepc4py.get_include()]
	
	os.environ["CC"] = "mpicc"
	os.environ["F90"] = "mpif90"
	os.environ["LDSHARED"] = "mpif90"
	config.add_extension('libpyctqw_MPI',
				 sources = ['src/libctqw-MPI.F90','src/libctqw-MPI.pyf'],
				 f2py_options=['--quiet'],
				 extra_f90_compile_args=['-Wno-unused-variable'],
				 #define_macros=[('F2PY_REPORT_ON_ARRAY_COPY',1)],
				 include_dirs=INCLUDE_DIRS + [os.curdir],
				 libraries=LIBRARIES,
				 library_dirs=LIBRARY_DIRS,
				 runtime_library_dirs=LIBRARY_DIRS)
	return config

def configuration(parent_package='', top_path=''):
	from numpy.distutils.misc_util import Configuration

	config = Configuration(None, parent_package, top_path)
	
	try:
		MPIpackage(config)
	except:
		print "Can't find PETSc/SLEPc/petsc4py!\
			\npyCTQW.MPI will not be installed"
	
	config.add_data_files(  ('pyCTQW/examples', 'examples/*'),
				('pyCTQW/graphs', 'graphs/*') )
	
	return config


def run_setup():
	from numpy.distutils.core import setup

	setup(
		name='pyCTQW',
		version='0.1.0',
		author='Josh Izaac',
		author_email='josh.izaac@uwa.edu.au',
		packages=['pyCTQW', 'pyCTQW.MPI'],
		url='http://pypi.python.org/pypi/pyCTQW/',
		#license='LICENSE.txt',
		description='An MPI enabled CTQW simulator',
		#long_description=open('README.txt').read(),
		provides=["pyctqw"],
#		scripts=["examples/exampleMPI.py","examples/example.py"],
		**(configuration(top_path='').todict())
	)

if __name__ == '__main__':
	run_setup()

