==================================
Quickstart Guide
==================================

Creating Python Scripts
=======================

Once :mod:`pyCTQW.MPI` is installed, calculations involving CTQWs can be called in several ways:

	* as an executable Python script
	* as an interactive session, using `IPython <http://ipython.org/>`_ or `IPython Notebook <http://ipython.org/notebook.html>`_

The syntax in both cases are identical, however this guide will be deal with executable Python scripts.

|

Initialization
---------------

The first step is to initialize the PETSc environment, and import the :mod:`pyCTQW.MPI` module:

::

	import sys, petsc4py
	petsc4py.init(sys.argv)
	from petsc4py import PETSc
	import pyCTQW.MPI as qw

|

.. _cmd-line-options:

Command Line Options
---------------------

Next, PETSc can be used to create command line options for the script; for example, the following code creates two command line options, :f:var:`t` and :f:var:`N`, with default values of ``100`` and ``20`` respectively:

::

	OptDB = PETSc.Options()
	N = OptDB.getInt('N', 100)
	t = OptDB.getReal('t', 20)

Furthermore, most PETSc and SLEPc subroutines accept command line options which modify their settings; for instance, when using the SLEPc EPS eigensolver, the EPS type can be changed dynamically when run through the following options:

.. code-block:: bash

    $ mpirun -np 2 <program> -eps_type='lapack'

.. seealso::
    For more details, please refer to the `PETSc <http://www.mcs.anl.gov/petsc/>`_  and `SLEPc <http://www.grycap.upv.es/slepc/>`_ documentation.

|

Rank and I/O
---------------------

When running on multiple nodes, we sometimes want only one node to perform a calculation or operation -- for instance, for I/O operations where all nodes already have the same information. Using PETSc, the ``rank`` (MPI process number) can be determined for each process, and conditional statements used to control which node performs the I/O operation:

::
	
	rank =  PETSc.Comm.Get_rank(PETSc.COMM_WORLD)

	if rank == 0:
		print '1P Line\n'

.. caution::
	* Be careful though: all of the methods and functions available in :mod:`pyCTQW.MPI` are designed to work globally on *all* processes, and should not be created or called on a subset of all available nodes. Doing so may cause your program to hang.
	* Most objects in :mod:`pyCTQW.MPI` contain I/O methods, for instance :func:`pyCTQW.MPI.Graph.exportState`; these are global over all nodes (as mentioned above) and should be used over custom methods when possible.
	* I/O operations involving PETSc objects (for instance :func:`PETSc.Vec.view()`) are also designed to be global/collective over all nodes; see `petsc4py <https://pypi.python.org/pypi/petsc4py/3.4>`_ for relevant documentation.

|

Creating CTQW objects
---------------------

All available objects and methods are detailed over at the API documentation pages :mod:`pyCTQW.MPI`.

.. seealso::
	You can also refer to the :ref:`examples <examples>` to see how CTQW objects can be manipulated

|

Running your script
====================

Once your script is complete, save it with a ``.py`` extension, and make it executable by running

.. code-block:: bash
	
	$ chmod +x <script>.py

where ``<script>`` is the file path of your script. Then, to run your program, simply run in a terminal

.. code-block:: bash
    
    $ mpirun -np X <script>.py [options]

where ``X`` is the number of MPI processes to run. 

|

.. _log-stage:

Code profiling
===============

PETSc also allows for easy code profiling by supplying the command line option :f:var:`-log_summary` when executing your script, and this is built in to :mod:`pyCTQW.MPI`; for instance, the methods available in :func:`pyCTQW.MPI.Line2P` automatically create log stages for creating the Hamiltonian, initial state, finding the eigenvalues, propagation etc.

If you wish to create custom log stages, this can also be done:

::
	
	stage1 = _PETSc.Log.Stage('First Stage')
	stage1.push()
	# place stage 1 functions/operations here
	stage1.pop()

.. seealso::
    For more details, please refer to the `petsc4py <https://pypi.python.org/pypi/petsc4py/3.4>`_ documentation.

