Fortran library and Python module to calculate continuous-time quantum walks

This is intended to provide a framework to quickly and easily work with quantum
walkers, take advantage of high performance computing, and allow easy visualisation.

Features:
============
* Fortran and Python bindings (in the form of a library and module respectively)
* Supports MPI through the use of the PETSc and SLEPc high-performance sparse
  matrix libraries (CUDA support planned)
* Has in-built support for infinite-line hamiltonians
* Import and export matrices/states in binary or text 
* Can import custom adjancency matrices
* Supports both one and two walkers (non-interacting)
* Import and export matrices/states in binary, text or matlab format
* Python module supports plotting and visualisation using matplotlib and networkx
* Entanglement calculations
* Ability to place diagonal defects/barriers on graph vertices
* MPI graph isomorphism methods, for both 2 and 3 interacting particles

TODO:
============
* CUDA support
* Calculation of transmission
* Add fallback modes when PETSc/SLEPc are not present

Dependencies
============

In addition to an MPI implementation (e.g. `MPICH <http://www.mpich.org/>`_ or `Open MPI <http://www.open-mpi.org/>`_), the python :py:mod:`pyCTQW.MPI` module depends on the following components:

	- `Python <http://www.python.org/>`_ >= 2.7
	- `NumPy <http://www.numpy.org/>`_ >= 1.6.0
	- `PETSc <http://www.mcs.anl.gov/petsc/>`_ >= 3.4.2	
	- `SLEPc <http://www.grycap.upv.es/slepc/>`_ >= 3.4.1	
	- `petsc4py 3.4 or petsc4py-dev <https://pypi.python.org/pypi/petsc4py/3.4>`_
	- `mpi4py <http://mpi4py.scipy.org/>`_		(recommended, used for some plotting)
	- `matplotlib <http://matplotlib.org/>`_	(recommended, for node plotting and graph visualisation)
	- `SciPy <http://www.scipy.org/>`_			(recommended, for some I/O operations)
	- `NetworkX <http://networkx.github.io/>`_		(recommended, graph visualisation)


Installation
============

1) Ensure all dependencies required above are installed

2) Extract the :mod:`pyCTQW` tarball, and ``cd`` into the extracted directory:

	.. code-block:: bash
		
		$ tar -zxvf pyCTQW-X.Y.tar.gz
		$ cd pyCTQW-X.Y

3) Ensure that your PETSc and SLEPc environment variables are correctly set; for example,

	.. code-block:: bash

		$ export PETSC_DIR=/path/to/petsc
		$ export PETSC_ARCH=linux-gnu
		$ export SLEPC_DIR=/path/to/slepc

	If you are unsure what your PETSc or SLEPc variables should be, please refer to their documentation.

	.. important::
		If you plan to install :py:mod:`pyCTQW.MPI` using ``root`` to a **system** directory, the PETSc and SLEPc environment variables must be available to the root user.

4) Compile the Python module :py:mod:`pyCTQW.MPI` by running

	.. code-block:: bash
		
		$ python setup.py build

5) System-wide install:

	.. code-block:: bash
		
		$ sudo -E python setup.py install

	where the command ``-E`` ensures that the environment variables set in step 3 are passed to the root.

	.. note::
		If you do not have root access, or the above command does not appear to work, you can install the package locally by running

			.. code-block:: bash
				
				$ python setup.py install --user

	Now, have a go running some of the :doc:`examples`!

*Optional*: compiling :mod:`libctqwMPI`
===========================================================

In addition to an MPI implementation (e.g. `MPICH <http://www.mpich.org/>`_ or `Open MPI <http://www.open-mpi.org/>`_), the Fortran library :mod:`libctqwMPI` depends on the following components:
    - `PETSc <http://www.mcs.anl.gov/petsc/>`_ >= 3.4.2   
    - `SLEPc <http://www.grycap.upv.es/slepc/>`_ >= 3.4.1

Once these dependencies are installed, simply open a terminal in the root directory of :file:`pyCTQW-X.Y` and run
     
    .. code-block:: bash    
        
        $ make fortran [options]

where available options include

.. list-table::
    :widths: 3 3 30
    :header-rows: 1

    * - Option
      - Values
      - Description

    * - :attr:`shared_lib`
      - 0 (default), 1
      - whether to build :f:mod:`libctqwMPI` as a shared library (:attr:`shared_lib=1`, producing :file:`libctqwMPI.so`) or a static library (:attr:`shared_lib=0` (default), producing :file:`libctqwMPI.a`).

        If built as a shared library, compiled programs will be smaller, but :file:`libctqwMPI.so` will need to be added to a directory used by ``ld`` (either by setting the environment variable ``LD_LIBRARY_PATH`` or by placing :file:`libctqwMPI.so` in :file:`/usr/local/lib` etc).
        
The fortran library (:file:`libctqwMPI.so` or :file:`libctqwMPI.a`) can be found in the :file:`pyCTQW-X.Y/lib` directory, with required module files found in the :file:`pyCTQW-X.Y/include` directory.


**Optional:** build documentation 
=======================================

If `Sphinx <http://sphinx-doc.org/>`_ is installed, the documentation can be compiled by running

	.. code-block:: bash
		
		$ make docs-html

Documentation
===============

.. seealso::
	For more information on how to use this package, please see the `online documentation <http://pyctqw.readthedocs.org>`_