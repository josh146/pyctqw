
.. f:module:: libctqwMPI
	:synopsis: Fortran interface to CTQW and Graph Isomorphism routines

.. f:currentmodule:: ctqwmpi

===================================
:f:mod:`libctqwMPI` API
===================================

(source: :file:`libctqw-MPI.F90`)

.. f:automodule: ctqwmpi

I/O subroutines
----------------

.. f:autofunction:: exportVec

.. f:autofunction:: importVec

.. f:autofunction:: exportMat

.. f:autofunction:: importMat

Matrix Array subroutines
-------------------------

.. f:autofunction:: identity

.. f:autofunction:: kron


Hamiltonian subroutines
-------------------------

.. f:autofunction:: importAdjToH

.. f:autofunction:: adjToH

.. f:autofunction:: hamiltonian_1p_line

.. f:autofunction:: hamiltonian_3p_line

 
Statespace subroutines
-------------------------------------

.. f:autofunction:: coord

.. f:autofunction:: coord3P

.. f:autofunction:: p1prob

.. f:autofunction:: marginal

.. f:autofunction:: marginal3

.. f:autofunction:: p1_init

.. f:autofunction:: p2_init

.. f:autofunction:: p3_init

 
MatrixExp and Eigenvalues
-------------------------------------


.. f:autofunction:: expm

.. f:autofunction:: min_max_eigs

.. f:autofunction:: qw_cheby

 
Entanglement subroutines
-------------------------------------


.. f:autofunction:: partial_trace_array

.. f:autofunction:: partial_trace_mat

.. f:autofunction:: entanglement


GraphIso subroutines
-------------------------------------



.. f:autofunction:: number_of_edges

.. f:autofunction:: getEdgeState

.. f:autofunction:: getAllEdgeStates

.. f:autofunction:: getAllEdgeStates3P

.. f:autofunction:: GraphISCert

.. f:autofunction:: d_refsor


