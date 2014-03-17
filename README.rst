Fortran library and Python module to calculate continuous-time quantum walks

This is intended to provide a framework to quickly and easily work with quantum
walkers, take advantage of high performance computing, and allow easy visualisation.

For more information on how to use this package, please see the `online documentation <http://pyctqw.readthedocs.org>`_

Features
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

TODO
============
* CUDA support
* Calculation of transmission
* Add fallback modes when PETSc/SLEPc are not present

Dependencies
============

In addition to an MPI implementation (e.g. `MPICH <http://www.mpich.org/>`_ or `Open MPI <http://www.open-mpi.org/>`_), the python ``pyCTQW.MPI`` module depends on the following components:

- `Python <http://www.python.org/>`_ >= 2.7
- `NumPy <http://www.numpy.org/>`_ >= 1.6.0
- `PETSc <http://www.mcs.anl.gov/petsc/>`_ >= 3.4.2	
- `SLEPc <http://www.grycap.upv.es/slepc/>`_ >= 3.4.1	
- `petsc4py 3.4 or petsc4py-dev <https://bitbucket.org/petsc/petsc4py>`_
- `mpi4py <http://mpi4py.scipy.org/>`_		(recommended, used for some plotting)
- `matplotlib <http://matplotlib.org/>`_	(recommended, for node plotting and graph visualisation)
- `SciPy <http://www.scipy.org/>`_			(recommended, for some I/O operations)
- `NetworkX <http://networkx.github.io/>`_		(recommended, graph visualisation)


Installation using `pip`
===========================

After ensuring NumPy and petsc4py are installed (and all PETSc, SLEPc and MPI environment variables are properly set), pyCTQW can be installed using `pip`:

.. code-block:: bash
	
	$ pip install pyCTQW


Installation from source code
==============================

Alternatively, the source code can be downloaded and compiled manually:

1) Ensure all dependencies required above are installed

2) Extract the ``pyCTQW`` folder, and ``cd`` into the extracted directory:

	.. code-block:: bash
		
		$ tar xvzf pyctqw-1.0.0.tar.gz
		$ cd pyctqw-1.0.0

3) Ensure that your PETSc and SLEPc environment variables are correctly set; for example,

	.. code-block:: bash

		$ export PETSC_DIR=/path/to/petsc
		$ export PETSC_ARCH=linux-gnu
		$ export SLEPC_DIR=/path/to/slepc

	If you are unsure what your PETSc or SLEPc variables should be, please refer to their documentation.

	.. important::
		If you plan to install ``pyCTQW.MPI`` using ``root`` to a **system** directory, the PETSc and SLEPc environment variables must be available to the root user.

4) Compile the Python module ``pyCTQW.MPI`` by running

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

	Now, have a go running some of the examples!

*Optional*: compiling ``libctqwMPI``
===========================================================

In addition to an MPI implementation (e.g. `MPICH <http://www.mpich.org/>`_ or `Open MPI <http://www.open-mpi.org/>`_), the Fortran library ``libctqwMPI`` depends on the following components:
- `PETSc <http://www.mcs.anl.gov/petsc/>`_ >= 3.4.2   
- `SLEPc <http://www.grycap.upv.es/slepc/>`_ >= 3.4.1

Once these dependencies are installed, simply open a terminal in the root directory of ``pyCTQW-X.Y`` and run
     
.. code-block:: bash    
    
    $ make fortran [options]

where available options include

.. list-table::
    :widths: 3 3 30
    :header-rows: 1

    * - Option
      - Values
      - Description

    * - ``shared_lib``
      - 0 (default), 1
      - whether to build ``libctqwMPI`` as a shared library (``shared_lib=1``, producing ``libctqwMPI.so``) or a static library (``shared_lib=0`` (default), producing ``libctqwMPI.a``).

        If built as a shared library, compiled programs will be smaller, but ``libctqwMPI.so`` will need to be added to a directory used by ``ld`` (either by setting the environment variable ``LD_LIBRARY_PATH`` or by placing ``libctqwMPI.so`` in ``/usr/local/lib`` etc).
        
The fortran library (``libctqwMPI.so`` or ``libctqwMPI.a``) can be found in the ``pyCTQW-X.Y/lib`` directory, with required module files found in the ``pyCTQW-X.Y/include`` directory.


*Optional:* build documentation 
=======================================

If `Sphinx <http://sphinx-doc.org/>`_ is installed, the documentation can be compiled by running

.. code-block:: bash
	
	$ make docs-html

Documentation
===============

For more information on how to use this package, please see the `online documentation <http://pyctqw.readthedocs.org>`_

Acknowledgements
===========================

The graph isomorphism subroutine ``GraphISCert`` uses the external subroutine ``d_refsor``, a highly optimised Fortran sorting implementation written by Michel Olagnon and part of the `ORDERPACK 2.0 <http://www.fortran-2000.com/rank/>`_ suite of ranking and sorting algorithms for Fortran 90.