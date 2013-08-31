==========================
:class:`pyCTQW.MPI.Line3P`
==========================

.. currentmodule:: pyCTQW.MPI

.. autoclass:: pyCTQW.MPI.Line3P
	:show-inheritance:

Method Summary
---------------

.. autosummary::
	:toctree: stub
	:template: method.rst

	Line3P.createH
	Line3P.createInitState
	Line3P.destroy
	Line3P.exportState
	Line3P.importInitState
	Line3P.plot
	Line3P.plotEntanglement
	Line3P.plotNode
	Line3P.plotNodes
	Line3P.propagate
	Line3P.psiToInit
	Line3P.watch

Attribute Summary
------------------

.. list-table::
	:widths: 15 10 30
	:header-rows: 0

	* - *Attribute*
	  - *Type*
	  - *Description*
	* - :attr:`Line3P.psi0`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N^3` element PETSc vector containing the initial state
	* - :attr:`Line3P.psi`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N^3` element PETSc vector containing the final state after the last propagation.
	* - :attr:`Line3P.psiX`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N` element PETSc vector containing the marginal propability of particle 1
	* - :attr:`Line3P.psiY`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N` element PETSc vector containing the marginal propability of particle 2
	* - :attr:`Line3P.psiZ`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N` element PETSc vector containing the marginal propability of particle 3
	* - :attr:`Line3P.H`
	  - :func:`pyCTQW.MPI.ctqw.Hamiltonian`
	  - Hamiltonian matrix
	* - :attr:`Line3P.EigSolver`
	  - :func:`pyCTQW.MPI.ctqw.EigSolver`
	  - The Hamiltonian eigensolver
	* - :attr:`Line3P.handle`
	  - :func:`pyCTQW.MPI.ctqw.nodeHandle`
	  - A node handle, created if nodes are being watched for probability.
	* - :attr:`Line3P.entanglementHandle`
	  - :func:`pyCTQW.MPI.ctqw.entanglementHandle`
	  - A node handle, created if the entanglement is being watched.
	
|