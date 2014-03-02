==================================
Installation
==================================

.. note::
	At the moment, only Unix based systems are supported.

Dependencies
============

In addition to an MPI implementation (e.g. `MPICH <http://www.mpich.org/>`_ or `Open MPI <http://www.open-mpi.org/>`_), the :py:mod:`pyCTQW.MPI` module depends on the following components:

	- `Python <http://www.python.org/>`_ >= 2.7
	- `NumPy <http://www.numpy.org/>`_ >= 1.6.0
	- `PETSc <http://www.mcs.anl.gov/petsc/>`_ >= 3.4.2	
	- `SLEPc <http://www.grycap.upv.es/slepc/>`_ >= 3.4.1	
	- `petsc4py 3.4 or petsc4py-dev <https://pypi.python.org/pypi/petsc4py/3.4>`_
	- `mpi4py <http://mpi4py.scipy.org/>`_		(recommended, used for some plotting)
	- `matplotlib <http://matplotlib.org/>`_	(recommended, for node plotting and graph visualisation)
	- `SciPy <http://www.scipy.org/>`_			(recommended, for some I/O operations)
	- `NetworkX <http://networkx.github.io/>`_		(recommended, graph visualisation)


Compiling and installing :py:mod:`pyCTQW.MPI`
=============================================

1) Ensure all dependencies required above are installed
   
2) Download the latest version of :mod:`pyctqw`:
   	
   	.. code-block:: bash

   		$ wget https://github.com/josh146/pyctqw/archive/v1.0.zip -O pyctqw.zip

3) Extract the :mod:`pyCTQW` zip folder, and ``cd`` into the extracted directory:

	.. code-block:: bash
		
		$ unzip v1.0.zip
		$ cd pyctqw-v1.0

4) Ensure that your PETSc and SLEPc environment variables are correctly set; for example,

	.. code-block:: bash

		$ export PETSC_DIR=/path/to/petsc
		$ export PETSC_ARCH=linux-gnu
		$ export SLEPC_DIR=/path/to/slepc

	If you are unsure what your PETSc or SLEPc variables should be, please refer to their documentation.

	.. important::
		If you plan to install :py:mod:`pyCTQW.MPI` using ``root`` to a **system** directory, the PETSc and SLEPc environment variables must be available to the root user.

5) Compile the Python module :py:mod:`pyCTQW.MPI` by running

	.. code-block:: bash
		
		$ python setup.py build

6) System-wide install:

	.. code-block:: bash
		
		$ sudo -E python setup.py install

	where the command ``-E`` ensures that the environment variables set in step 3 are passed to the root.

	.. note::
		If you do not have root access, or the above command does not appear to work, you can install the package locally by running

			.. code-block:: bash
				
				$ python setup.py install --user

	Now, have a go running some of the :doc:`examples`!


**Optional:** build documentation 
=======================================

If `Sphinx <http://sphinx-doc.org/>`_ is installed, the documentation can be compiled by running

	.. code-block:: bash
		
		$ make docs-html

or

	.. code-block:: bash
		
		$ make docs-pdf

Note that in order to compile the PDF documentation, texlive must be
installed.

Known Issues
==============

* Non-mpi fallback modes not present yet