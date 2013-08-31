==========================
:class:`pyCTQW.MPI.Line2P`
==========================

.. currentmodule:: pyCTQW.MPI

.. autoclass:: pyCTQW.MPI.Line2P
	:show-inheritance:

Method Summary
---------------

.. autosummary::
	:toctree: stub
	:template: method.rst

	Line2P.createH
	Line2P.createInitState
	Line2P.destroy
	Line2P.exportPartialTrace
	Line2P.exportState
	Line2P.importInitState
	Line2P.plot
	Line2P.plotEntanglement
	Line2P.plotNode
	Line2P.plotNodes
	Line2P.propagate
	Line2P.psiToInit
	Line2P.watch


Attribute Summary
------------------

.. list-table::
	:widths: 15 10 30
	:header-rows: 0

	* - *Attribute*
	  - *Type*
	  - *Description*
	* - :attr:`Line2P.psi0`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N^2` element PETSc vector containing the initial state
	* - :attr:`Line2P.psi`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N^2` element PETSc vector containing the final state after the last propagation.
	* - :attr:`Line2P.psiX`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N` element PETSc vector containing the marginal propability of particle 1
	* - :attr:`Line2P.psiY`
	  - :func:`petsc4py.PETSc.Vec`
	  - :math:`N` element PETSc vector containing the marginal propability of particle 2
	* - :attr:`Line2P.H`
	  - :func:`pyCTQW.MPI.ctqw.Hamiltonian`
	  - Hamiltonian matrix
	* - :attr:`Line2P.EigSolver`
	  - :func:`pyCTQW.MPI.ctqw.EigSolver`
	  - The Hamiltonian eigensolver
	* - :attr:`Line2P.handle`
	  - :func:`pyCTQW.MPI.ctqw.nodeHandle`
	  - A node handle, created if nodes are being watched for probability.
	* - :attr:`Line2P.entanglementHandle`
	  - :func:`pyCTQW.MPI.ctqw.entanglementHandle`
	  - A node handle, created if the entanglement is being watched.
	
|