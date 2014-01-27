
.. f:module:: libctqwMPI
	:synopsis: Fortran interface to CTQW and Graph Isomorphism routines

.. f:currentmodule:: ctqwmpi

===================================
:f:mod:`libctqwMPI` API
===================================

(source: :file:`ctqwMPI.F90`)

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

.. f:autofunction:: hamiltonian_p1_line

.. f:autofunction:: hamiltonian_p2_line

.. f:autofunction:: hamiltonian_p3_line

 
Statespace subroutines
-------------------------------------

.. f:autofunction:: coord

.. f:autofunction:: coord3P

.. f:autofunction:: marginal1

.. f:autofunction:: marginal2

.. f:autofunction:: marginal3

.. f:autofunction:: p1_init

.. f:autofunction:: p2_init

.. f:autofunction:: p3_init

 
MatrixExp and Eigenvalues
-------------------------------------


.. f:autofunction:: qw_krylov

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

.. note::
	:f:func:`d_refsor`, a highly optimised Fortran sorting implementation written by Michel Olagnon and part of the `ORDERPACK 2.0 <http://www.fortran-2000.com/rank/>`_ suite of ranking and sorting algorithms for Fortran 90.

