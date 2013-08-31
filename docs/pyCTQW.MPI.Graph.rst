=========================
:class:`pyCTQW.MPI.Graph`
=========================

.. currentmodule:: pyCTQW.MPI

.. autoclass:: pyCTQW.MPI.Graph
	:show-inheritance:


Method Summary
---------------

.. autosummary::
	:toctree: stub
	:template: method.rst

	Graph.clearLiveGraph
	Graph.createH
	Graph.createInitState
	Graph.destroy
	Graph.exportState
	Graph.importInitState
	Graph.plot
	Graph.plotGraph
	Graph.plotLiveGraph
	Graph.plotNodes
	Graph.propagate
	Graph.psiToInit
	Graph.watch


Attribute Summary
------------------

.. list-table::
	:widths: 15 10 30
	:header-rows: 0

	* - *Attribute*
	  - *Type*
	  - *Description*
	* - :attr:`Graph.psi0`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N` element PETSc vector containing the initial state
	* - :attr:`Graph.psi`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N` element PETSc vector containing the final state after the last propagation.
	* - :attr:`Graph.H`
	  - :func:`pyCTQW.MPI.ctqw.Hamiltonian`
	  - Hamiltonian matrix
	* - :attr:`Graph.EigSolver`
	  - :func:`pyCTQW.MPI.ctqw.EigSolver`
	  - The Hamiltonian eigensolver
	* - :attr:`Graph.handle`
	  - :func:`pyCTQW.MPI.ctqw.nodeHandle`
	  - A node handle, created if nodes are being watched for probability.
