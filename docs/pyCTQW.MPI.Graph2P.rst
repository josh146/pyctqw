===========================
:class:`pyCTQW.MPI.Graph2P`
===========================

.. currentmodule:: pyCTQW.MPI

.. autoclass:: pyCTQW.MPI.Graph2P
	:show-inheritance:

Method Summary
---------------

.. autosummary::
	:toctree: stub
	:template: method.rst

	Graph2P.clearLiveGraph
	Graph2P.createH
	Graph2P.createInitState
	Graph2P.destroy
	Graph2P.exportPartialTrace
	Graph2P.exportState
	Graph2P.importInitState
	Graph2P.plot
	Graph2P.plotEntanglement
	Graph2P.plotGraph
	Graph2P.plotLiveGraph
	Graph2P.plotNode
	Graph2P.plotNodes
	Graph2P.propagate
	Graph2P.psiToInit
	Graph2P.watch


Attribute Summary
------------------

.. list-table::
	:widths: 15 10 30
	:header-rows: 0

	* - *Attribute*
	  - *Type*
	  - *Description*
	* - :attr:`Graph2P.psi0`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N^2` element PETSc vector containing the initial state
	* - :attr:`Graph2P.psi`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N^2` element PETSc vector containing the final state after the last propagation.
	* - :attr:`Graph2P.psiX`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N` element PETSc vector containing the marginal propability of particle 1
	* - :attr:`Graph2P.psiY`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N` element PETSc vector containing the marginal propability of particle 2
	* - :attr:`Graph2P.H`
	  - :func:`pyCTQW.MPI.ctqw.Hamiltonian`
	  - Hamiltonian matrix
	* - :attr:`Graph2P.EigSolver`
	  - :func:`pyCTQW.MPI.ctqw.EigSolver`
	  - The Hamiltonian eigensolver
	* - :attr:`Graph2P.handle`
	  - :func:`pyCTQW.MPI.ctqw.nodeHandle`
	  - A node handle, created if nodes are being watched for probability.
	* - :attr:`Graph2P.entanglementHandle`
	  - :func:`pyCTQW.MPI.ctqw.entanglementHandle`
	  - A node handle, created if the entanglement is being watched.
	
|