#!/usr/bin/env python2.7
import sys, os

info={
    'name'          : 'pyCTQW',
    'version'       : '1.0.0',
    'author'        : 'Josh Izaac',
    'author_email'  : 'josh@iza.ac',
    'url'           : 'http://pyctqw.readthedocs.org',
    'license'       : 'GPLv3',
    'description'   : 'An MPI enabled CTQW simulator',
    'long_description' : open('README.rst').read(),
    'provides'      : ["pyCTQW"],
    'install_requires' : ['numpy','scipy','matplotlib','networkx==1.7','mpi4py','petsc4py==3.4'],
    'dependency_links' : ['http://networkx.lanl.gov/download/networkx/networkx-1.7.tar.gz#egg=networkx-1.7',
        'https://bitbucket.org/petsc/petsc4py/downloads/petsc4py-3.4.tar.gz#egg=petsc4py-3.4']
  }

classifiers=[
    "Development Status :: 5 - Production/Stable",
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    "Natural Language :: English",
    "Operating System :: POSIX",
    "Operating System :: MacOS :: MacOS X",
    "Operating System :: POSIX :: Linux",
    "Programming Language :: Python",
    "Programming Language :: Fortran",
    "Topic :: Scientific/Engineering :: Physics"
  ]

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
    LIBRARIES   = []
    
    # Get PETSc environment variables
    try:
        PETSC_ARCH = os.environ.get('PETSC_ARCH', '')
        PETSC_DIR  = os.environ['PETSC_DIR']
        if PETSC_DIR == '':
            raise KeyError
        elif not (os.path.isfile(PETSC_DIR+'/'+PETSC_ARCH+'/lib/libpetsc.so') or os.path.isfile(PETSC_DIR+'/'+PETSC_ARCH+'/lib/libpetsc.a')):
            raise dirError
    except KeyError:
        raise Exception("ERROR: PETSC_DIR environment variable not set")
    except dirError:
        raise Exception("ERROR: PETSC_DIR does not point towards a valid directory")
    
    # get PETSc include and library directories
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
        raise Exception("""ERROR: petsc4py 3.4 not installed. 
This can be installed by running 
pip install https://bitbucket.org/petsc/petsc4py/downloads/petsc4py-3.4.tar.gz
Please see http://pythonhosted.org/petsc4py/ for more info.""")
        
    INCLUDE_DIRS += [petsc4py.get_include()]
                  
    # Get SLEPc environment variables
    try:
        SLEPC_DIR  = os.environ['SLEPC_DIR']
        if SLEPC_DIR == '':
            raise KeyError
        elif not (os.path.isfile(SLEPC_DIR+'/'+PETSC_ARCH+'/lib/libslepc.so') or os.path.isfile(SLEPC_DIR+'/'+PETSC_ARCH+'/lib/libslepc.a')):
            raise dirError
    except KeyError:
        raise Exception("ERROR: SLEPC_DIR environment variable not set")
    except dirError:
        raise Exception("ERROR: SLEPC_DIR does not point towards a valid directory")
    
    # get SLEPc include and library directories 
    if PETSC_ARCH and os.path.isdir(os.path.join(PETSC_DIR, PETSC_ARCH)):
        INCLUDE_DIRS += [os.path.join(SLEPC_DIR, PETSC_ARCH, 'include'),
                         os.path.join(SLEPC_DIR, 'include')]
        LIBRARY_DIRS += [os.path.join(SLEPC_DIR, PETSC_ARCH, 'lib')]
    else:
        INCLUDE_DIRS += [os.path.join(SLEPC_DIR, 'include')]
        LIBRARY_DIRS += [os.path.join(SLEPC_DIR, 'lib')]
    LIBRARIES += ['slepc']
    
    # set the compiler to mpi
    if os.getenv('LDSHARED', '') == '':
        os.environ["CC"] = "mpicc"
        os.environ["F90"] = "mpif90"
        os.environ["LDSHARED"] = "mpif90"
    
    # create the extension  
    config.add_extension('libpyctqw_MPI',
                 sources = ['src/ctqwMPI.F90','src/ctqwMPI.pyf'],
                 f2py_options=['--quiet'],
                 extra_f90_compile_args=['-Wno-unused-variable','-Wno-conversion','-Wno-unused-dummy-argument'],
                 #define_macros=[('F2PY_REPORT_ON_ARRAY_COPY',1)],
                 include_dirs=INCLUDE_DIRS + [os.curdir],
                 libraries=LIBRARIES,
                 library_dirs=LIBRARY_DIRS,
                 runtime_library_dirs=LIBRARY_DIRS)
    return config

def configuration(parent_package='', top_path=''):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    MPIpackage(config)
    mpi = True
    
    # if os.getenv('use_mpi','True') == 'True':
    #   try:
    #       MPIpackage(config)
    #       mpi = True
    #   except noMPI, err:
    #       warn(str(err))
    #       warn("WARNING: error when looking for PETSc/SLEPc/petsc4py!\
    #            \nIf they are installed, double check $PETSC_DIR, $SLEPC_DIR and $PETSC_ARCH environment varibles.\
    #            \n\npyCTQW.MPI will not be installed\n")
    #       mpi = False
    # else:
    #   mpi = False
    
    # config.add_data_files(  ('pyCTQW/examples', 'examples/*.py'),
    #           ('pyCTQW/graphs/cayley', 'graphs/cayley/*'),
    #           ('pyCTQW/graphs/strong-regular-25-12-5-6', 'graphs/strong-regular-25-12-5-6/*') )
    
    return config.todict(), mpi


def run_setup():
    try:
        import setuptools
        from numpy.distutils.core import setup
    except ImportError:
        raise Exception("Requires NumPy>=1.6")

    config, mpi = configuration(top_path='')
    if mpi:
        packages=['pyCTQW', 'pyCTQW.MPI']
    else:
        packages=['pyCTQW']

    config = dict(info,**config)

    setup(
        packages=packages,
        classifiers=classifiers,
        # ext_package='pyCTQW.MPI',
        **(config)
    )

if __name__ == '__main__':
    run_setup()

