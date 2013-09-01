#!/usr/bin/env python2.7
import sys, os

info={
	'name'			: 'pyCTQW',
	'version'		: '0.1.0',
	'author'		: 'Josh Izaac',
	'author_email' 	: 'josh.izaac@uwa.edu.au',
	'url' 			: 'http://pypi.python.org/pypi/pyCTQW/',
	'license' 		: 'LICENSE',
	'description' 	: 'An MPI enabled CTQW simulator',
	'long_description' : open('README').read(),
	'provides' 		: ["pyCTQW"]
  }

classifiers=[
	'Development Status :: 2 - Pre-Alpha',
	'Environment :: Console',
	'Intended Audience :: Science/Research',
	'License :: OSI Approved :: BSD License',
	'Natural Language :: English',
	'Operating System :: MacOS :: MacOS X',
	'Operating System :: Microsoft :: Windows',
	'Operating System :: POSIX',
	'Programming Language :: Python',
	'Programming Language :: Fortran',
	'Topic :: Scientific/Engineering :: Physics'
  ],

class noMPI(Exception):
	pass

class dirError(Exception):
	pass

def warn(msg=''):
	return sys.stderr.write(msg+'\n')

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
		PETSC_ARCH = os.environ.get('PETSC_ARCH', '')
		PETSC_DIR  = os.environ['PETSC_DIR']
		if PETSC_DIR == '':
			raise KeyError
		elif not os.path.isfile(PETSC_DIR+PETSC_ARCH+'/lib/libpetsc.so'):
			raise dirError
	except KeyError:
		raise noMPI("WARNING: PETSC_DIR environment variable not set")
	except dirError:
		raise noMPI("WARNING: PETSC_DIR does not point towards a valid directory")
	
	if PETSC_ARCH and os.path.isdir(os.path.join(PETSC_DIR, PETSC_ARCH)):
		INCLUDE_DIRS += [os.path.join(PETSC_DIR, PETSC_ARCH, 'include'),
						 os.path.join(PETSC_DIR, 'include')]
		LIBRARY_DIRS += [os.path.join(PETSC_DIR, PETSC_ARCH, 'lib')]
	else:
		if PETSC_ARCH: pass
		INCLUDE_DIRS += [os.path.join(PETSC_DIR, 'include')]
		LIBRARY_DIRS += [os.path.join(PETSC_DIR, 'lib')]
	LIBRARIES += ['petsc']

	# PETSc for Python
	try:
		import petsc4py
	except ImportError:
		raise noMPI("WARNING: petsc4py not installed")
		
	INCLUDE_DIRS += [petsc4py.get_include()]
				  
	# SLEPc
	
	try:
		SLEPC_DIR  = os.environ['SLEPC_DIR']
		if SLEPC_DIR == '':
			raise KeyError
		elif not os.path.isfile(SLEPC_DIR+PETSC_ARCH+'/lib/libslepc.so'):
			raise dirError
	except KeyError:
		raise noMPI("WARNING: SLEPC_DIR environment variable not set")
	except dirError:
		raise noMPI("WARNING: SLEPC_DIR does not point towards a valid directory")
	
	if PETSC_ARCH and os.path.isdir(os.path.join(PETSC_DIR, PETSC_ARCH)):
		INCLUDE_DIRS += [os.path.join(SLEPC_DIR, PETSC_ARCH, 'include'),
						 os.path.join(SLEPC_DIR, 'include')]
		LIBRARY_DIRS += [os.path.join(SLEPC_DIR, PETSC_ARCH, 'lib')]
	else:
		INCLUDE_DIRS += [os.path.join(SLEPC_DIR, 'include')]
		LIBRARY_DIRS += [os.path.join(SLEPC_DIR, 'lib')]
	LIBRARIES += ['slepc']
	
	if os.getenv('LDSHARED', '') == '':
		os.environ["CC"] = "mpicc"
		os.environ["F90"] = "mpif90"
		os.environ["LDSHARED"] = "mpif90"

	config.add_installed_library('refsor', ['src/sort.f90'], '.')
		
	config.add_extension('libpyctqw_MPI',
				 sources = ['src/sort.f90','src/libctqw-MPI.F90','src/libctqw-MPI.pyf'],
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
	
	if os.getenv('use_mpi','True') == 'True':
		try:
			MPIpackage(config)
			mpi = True
		except noMPI, err:
			warn(str(err))
			warn("WARNING: error when looking for PETSc/SLEPc/petsc4py!\
				 \nIf they are installed, double check $PETSC_DIR, $SLEPC_DIR and $PETSC_ARCH environment varibles.\
				 \n\npyCTQW.MPI will not be installed\n")
			mpi = False
	else:
		mpi = False
	
	# config.add_data_files(  ('pyCTQW/examples', 'examples/*.py'),
	# 			('pyCTQW/graphs/cayley', 'graphs/cayley/*'),
	# 			('pyCTQW/graphs/strong-regular-25-12-5-6', 'graphs/strong-regular-25-12-5-6/*') )
	
	return config.todict(), mpi


def run_setup():
	try:
		from numpy.distutils.core import setup
	except ImportError:
		raise DistutilsError("requires NumPy>=1.6")

	config, mpi = configuration(top_path='')
	if mpi:
		packages=['pyCTQW', 'pyCTQW.MPI']
	else:
		packages=['pyCTQW']

	config = dict(info,**config)

	setup(
		packages=packages,
		classifiers=classifiers,
		**(config)
	)

if __name__ == '__main__':
	run_setup()

