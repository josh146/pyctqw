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
- `mpi4py <http://mpi4py.scipy.org/>`_      (recommended, used for some plotting)
- `matplotlib <http://matplotlib.org/>`_    (recommended, for node plotting and graph visualisation)
- `SciPy <http://www.scipy.org/>`_          (recommended, for some I/O operations)
- `NetworkX <http://networkx.github.io/>`_      (recommended, graph visualisation)

Installation using `pip`
===========================

After ensuring NumPy and petsc4py are installed (and all PETSc, SLEPc and MPI environment variables are properly set), :mod:`pyctqw` can be installed using `pip`:

.. code-block:: bash
    
    $ pip install pyCTQW

.. note::
    The current development version `pyCTQW-dev <http://github.com/josh146/pyctqw/archive/master.tar.gz#egg=pyctqw-dev>`_ can also be installed using `pip`:

    .. code-block:: bash
    
        $ pip install --allow-unverified pyCTQW pyCTQW==dev

Installing :py:mod:`pyCTQW.MPI` from source code
=================================================

Alternatively, the source code can be downloaded and compiled manually:
   
1) :doc:`Download <downloads>` the latest version of :mod:`pyctqw`, extract the :mod:`pyCTQW` archive, and ``cd`` into the extracted directory:
    
    .. code-block:: bash

        $ wget https://github.com/josh146/pyctqw/archive/1.0.0.tar.gz -O pyctqw-1.0.0.tar.gz
        $ tar xzvf pyctqw-1.0.0.tar.gz
        $ cd pyctqw-1.0.0

2) Ensure that your PETSc and SLEPc environment variables are correctly set; for example,

    .. code-block:: bash

        $ export PETSC_DIR=/path/to/petsc
        $ export PETSC_ARCH=linux-gnu
        $ export SLEPC_DIR=/path/to/slepc

    If you are unsure what your PETSc or SLEPc variables should be, please refer to their documentation.

    .. important::
        If you plan to install :py:mod:`pyCTQW.MPI` using ``root`` to a **system** directory, the PETSc and SLEPc environment variables must be available to the root user.

3) Compile the Python module :py:mod:`pyCTQW.MPI` by running

    .. code-block:: bash
        
        $ python setup.py build

4) System-wide install:

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
    
    $ make docs

Known Issues
==============

* Non-mpi fallback modes not present yet
