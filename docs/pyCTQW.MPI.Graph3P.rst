===========================
:class:`pyCTQW.MPI.Graph3P`
===========================

.. currentmodule:: pyCTQW.MPI

.. autoclass:: pyCTQW.MPI.Graph3P
	:show-inheritance:

Method Summary
---------------

.. autosummary::
	:toctree: stub
	:template: method.rst

	Graph3P.createH
	Graph3P.createInitState
	Graph3P.destroy
	Graph3P.exportState
	Graph3P.importInitState
	Graph3P.plot
	Graph3P.plotEntanglement
	Graph3P.plotNode
	Graph3P.plotNodes
	Graph3P.propagate
	Graph3P.psiToInit
	Graph3P.watch


Attribute Summary
------------------

.. list-table::
	:widths: 15 10 30
	:header-rows: 0

	* - *Attribute*
	  - *Type*
	  - *Description*
	* - :attr:`Graph3P.psi0`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N^3` element PETSc vector containing the initial state
	* - :attr:`Graph3P.psi`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N^3` element PETSc vector containing the final state after the last propagation.
	* - :attr:`Graph3P.psiX`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N` element PETSc vector containing the marginal propability of particle 1
	* - :attr:`Graph3P.psiY`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N` element PETSc vector containing the marginal propability of particle 2
	* - :attr:`Graph3P.psiZ`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N` element PETSc vector containing the marginal propability of particle 3
	* - :attr:`Graph3P.H`
	  - :func:`pyCTQW.MPI.ctqw.Hamiltonian`
	  - Hamiltonian matrix
	* - :attr:`Graph3P.EigSolver`
	  - :func:`pyCTQW.MPI.ctqw.EigSolver`
	  - The Hamiltonian eigensolver
	* - :attr:`Graph3P.handle`
	  - :func:`pyCTQW.MPI.ctqw.nodeHandle`
	  - A node handle, created if nodes are being watched for probability.
	* - :attr:`Graph3P.entanglementHandle`
	  - :func:`pyCTQW.MPI.ctqw.entanglementHandle`
	  - A node handle, created if the entanglement is being watched.