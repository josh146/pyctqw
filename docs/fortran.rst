==================================
Fortran Library
==================================

.. caution::
	It is significantly more difficult to compile and link programs against the Fortran library directly, due to the myriad of pre-processing and the number of required headers/environment variables. As such, it is recommended users interface with this library via the much easier to use Python module :py:mod:`pyCTQW.MPI`.

	However, if you have experience with using PETSc and SLEPc with Fortran, and wish to proceed. A makefile template is provided that can be modified to link your program against :f:mod:`libctqwMPI`.

Dependencies
============

In addition to an MPI implementation (e.g. `MPICH <http://www.mpich.org/>`_ or `Open MPI <http://www.open-mpi.org/>`_), the Fortran library :f:mod:`libctqwMPI` depends on the following components:

	- `PETSc <http://www.mcs.anl.gov/petsc/>`_ >= 3.4.2	
	- `SLEPc <http://www.grycap.upv.es/slepc/>`_ >= 3.4.1


Compiling :f:mod:`libctqwMPI`
===========================================================

If you wish to write programs linking directly to the Fortran library, follow steps 1-3 above, before simply open a terminal in the root directory of :file:`pyCTQW-X.Y` and run
     
    .. code-block:: bash    
        
        $ make fortran
        
The fortran libraries :file:`libctqwMPI.so` and :file:`librefsor.a` can be found in the :file:`pyCTQW-X.Y/lib` directory, with required interfaces found in the :file:`pyCTQW-X.Y/include` directory.


Using :f:mod:`libctqwMPI`
===========================

Makefile
---------

Functions Documentation
========================
.. toctree::
	
	libctqw-MPI


Known Issues
==============

* Non-mpi fallback modes not present yet