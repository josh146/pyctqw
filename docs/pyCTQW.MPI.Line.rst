=========================
:class:`pyCTQW.MPI.Line`
=========================

.. currentmodule:: pyCTQW.MPI

.. autoclass:: pyCTQW.MPI.Line
	:show-inheritance:

Method Summary
---------------

.. autosummary::
	:toctree: stub
	:template: method.rst

	Line.createH
	Line.createInitState
	Line.destroy
	Line.exportState
	Line.importInitState
	Line.plot
	Line.plotNodes
	Line.propagate
	Line.psiToInit
	Line.watch


Attribute Summary
------------------

.. list-table::
	:widths: 15 10 30
	:header-rows: 0

	* - *Attribute*
	  - *Type*
	  - *Description*
	* - :attr:`Line.psi0`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N` element PETSc vector containing the initial state
	* - :attr:`Line.psi`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N` element PETSc vector containing the final state after the last propagation.
	* - :attr:`Line.H`
	  - :func:`pyCTQW.MPI.ctqw.Hamiltonian`
	  - Hamiltonian matrix
	* - :attr:`Line.EigSolver`
	  - :func:`pyCTQW.MPI.ctqw.EigSolver`
	  - The Hamiltonian eigensolver
	* - :attr:`Line.handle`
	  - :func:`pyCTQW.MPI.ctqw.nodeHandle`
	  - A node handle, created if nodes are being watched for probability.