==================================
Fortran Library
==================================

.. list-table::
    :widths: 10 30


    * - .. admonition:: Function documentation

            .. toctree::
    
                libctqw-MPI
      - .. caution::                                                                                                                                                   
            It is slightly more difficult to compile and link programs against the Fortran library directly, due to the myriad of pre-processing and the number of \ 
            required headers/environment variables. As such, it is recommended users interface with this library via the much easier to use Python module :py:mod:`pyCTQW.MPI`.
            
            However, if you have experience using PETSc and SLEPc with Fortran, and wish to proceed, a makefile template is provided that \ 
            can be modified to link your program against :f:mod:`libctqwMPI`.


Compiling :f:mod:`libctqwMPI`
===========================================================

In addition to an MPI implementation (e.g. `MPICH <http://www.mpich.org/>`_ or `Open MPI <http://www.open-mpi.org/>`_), the Fortran library :f:mod:`libctqwMPI` depends on the following components:
    - `PETSc <http://www.mcs.anl.gov/petsc/>`_ :math:`\geq` 3.4.2   
    - `SLEPc <http://www.grycap.upv.es/slepc/>`_ :math:`\geq` 3.4.1

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

|


Using :f:mod:`libctqwMPI`
=============================

To call functions and subroutines from :f:mod:`libctqwMPI`, your program should have the following structure:

.. code-block :: fortran
    :linenos:
    
    program main
        ! interfaces for libctqwMPI
        use ctqwMPI

    ! PETSc headers required by the preprocessor.
    ! Note the lack of indentation.
    #include <finclude/petsc.h>

        PetscErrorCode :: ierr
        PetscMPIInt    :: rank
        !--------------------
        ! your variable
        ! declarations 
        ! -------------------

        ! initialize SLEPc and PETSc
        call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
        call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)

        !--------------------
        ! your program code
        ! -------------------

        ! finalise PETSc
        call PetscFinalize(ierr)

    end program main

CTQW subroutines can now be called directly; for example,

.. code-block :: fortran

    call qw_cheby(psi0,psi,t,H,0.*PETSc_i,Emax,rank,n)

.. seealso::

    * For more details on the available subroutines, see :f:mod:`libctqwMPI`.
    * For an example of a Fortran program calling :f:mod:`libctqwMPI`, see :ref:`exampleMPI.F90 <example-MPI>`.

|

*Command Line Options*
-----------------------

Similarly to the :ref:`Python scripts <cmd-line-options>`, PETSc allows command line arguments to be easily added to your Fortran programs, through the use of the following subroutines:

.. code-block:: fortran

    call PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-n",n,flag,ierr)
    CHKERRQ(ierr)
    if (flag .eqv. .false.) n = 100

In this case, the command line option :attr:`-n` is added, and the command stored in variable :f:var:`n` (of type :f:type:`PETSC_INT`). Furthermore, if the flag is not provided, the program assigns :f:var:`n` a default value of ``100``.

Furthermore, most PETSc and SLEPc subroutines accept command line options which modify their settings; for instance, when using the SLEPc EPS eigensolver, the EPS type can be changed dynamically when run through the following options:

.. code-block:: bash

    $ mpirun -np 2 <program> -eps_type='lapack'

.. seealso::
    For more details, please refer to the `PETSc <http://www.mcs.anl.gov/petsc/>`_  and `SLEPc <http://www.grycap.upv.es/slepc/>`_ documentation.

|

*Code Profiling*
-----------------

PETSc also allows for easy code profiling using :f:type:`PetscLogStage` variables; however, unlike the built in log stages found in :ref:`pyCTQW.MPI <log-stage>`, these must be manually created for fortran programs.

For example, to create enable profiling of they Chebyshev method used for CTQW propagation, we would declare a log stage variable :f:var:`stage`,

.. code-block:: fortran

    PetscLogStage  :: stage
    PetscLogDouble :: tc0, tc1
    PetscErrorCode :: ierr

and then wrap the Chebyshev subroutine like so:

.. code-block:: fortran

    ! measure the initial wall time
    call PetscTime(tc0,ierr)

    ! register the log stage, named 'Chebyshev'
    call PetscLogStageRegister('Chebyshev',stage,ierr)
    call PetscLogStagePush(stage,ierr)

    ! CTQW propagation
    call qw_cheby(psi0,psi,t,H,0.*PETSc_i,Emax,rank,n)

    ! view the final statespace
    call VecView(psi,PETSC_VIEWER_STDOUT_WORLD,ierr)

    ! wait here until all MPI nodes arrive
    call PetscBarrier(psi,ierr)
    ! end the log stage
    call PetscLogStagePop(ierr)

    ! measure the end wall time
    call PetscTime(tc1,ierr)

Then, when calling the program with the :f:var:`-log_summary` command line option,

.. code-block:: bash

    $ mpirun -np 2 <program> -log_summary

the produced code profile will include a summary of performance within the 'Chebyshev' log stage.

.. seealso::
    For more details, please refer to the `PETSc <http://www.mcs.anl.gov/petsc/>`_  and `SLEPc <http://www.grycap.upv.es/slepc/>`_ documentation.

|

Linking :f:mod:`libctqwMPI`
================================

To compile your program and link it against :f:mod:`libctqwMPI`, you will need to include the :file:`ctqw_common` file (present in the root directory of the :file:`pyCTQW-X.Y` folder) in your makefile - this helps set up the required library/include variables.

For example, consider the makefile below (used to compile :ref:`exampleMPI.F90 <example-MPI>`):

*Makefile*
-----------
:download:`[Download makefile] <../examples/fortran/makefile>`

.. literalinclude:: ../examples/fortran/makefile
    :linenos:
    :language: makefile

A brief summary of the variables:

.. list-table::
    :widths: 3 30

    * - ``CTQW_DIR``
      - the path to the :file:`pyCTQW-X.Y` directory (used to correctly set up ``CTQW_LIBS`` and ``CTQW_INCLUDE_DIRS``).
    * - ``CTQW_LIBS``
      - the libraries required to link against :f:mod:`libctqwMPI`; these include :file:`libctqwMPI.so` (by default, this is compiled to ``CTQW_DIR/lib``) as well as the PETSc and SLEPc libraries.

        These are required when linking your object files.
    * - ``CTQW_INCLUDE_DIRS``
      - the location of the header/module files required - by default, these are compiled to ``CTQW_DIR/include``. Note that this also contains required PETSc/SLEPc include directories.
        
        These are required when compiling your code to an object file.

|

Running your program
===========================

To run your program, simply run in a terminal

.. code-block:: bash
    
    $ mpirun -np X <program> [options]

where ``X`` is the number of MPI nodes you would like to run it on. 

.. important::
    If you compiled your program using the shared library :file:`libctqwMPI.so`, make sure that the shared library is available at runtime, either by:

        * setting the environment variable ``LD_LIBRARY_PATH`` to include the directory containing :file:`libctqwMPI.so`,
        * or by placing :file:`libctqwMPI.so` in a directory used by ``ld`` (e.g. :file:`/usr/local/lib` etc).

Acknowledgements
===========================

The graph isomorphism subroutine :f:func:`graphiscert` uses the external subroutine :f:func:`d_refsor`, a highly optimised Fortran sorting implementation written by Michel Olagnon and part of the `ORDERPACK 2.0 <http://www.fortran-2000.com/rank/>`_ suite of ranking and sorting algorithms for Fortran 90.