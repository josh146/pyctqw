#!/usr/bin/env python2.7
# initialize PETSc
import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np

# import pyCTQW as qw
import pyCTQW.MPI as qw

# get the MPI rank
rank =  PETSc.Comm.Get_rank(PETSc.COMM_WORLD)

# create a graph isomorphism object
gi = qw.GraphISO()

# create a comparison table of two 3-cayley permutations
# present in the '../graphs/cayley/3-cayley-permutations'
# folder. As they are permutations, they are isomorphic,
# so the result should be a 2x2 matrix composed of 1.
comparisonTable = gi.AllIsomorphicQ('../graphs/cayley/3-cayley-permutations',info=False)

if rank==0:
	print '1) Testing isomorphism of all pairings:'
	print comparisonTable


# create a comparison table of a 3-cayley graph, and a
# a 3-cayley graph with an additional edge added.
# These are *not* isomorphic,
# so the result should be a 2x2 matrix identity matrix.
comparisonTable = gi.AllIsomorphicQ('../graphs/cayley/3-cayley-variant',info=False)

if rank==0:
	print '\n2) Testing isomorphism of all pairings:'
	print comparisonTable