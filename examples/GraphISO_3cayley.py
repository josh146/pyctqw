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
# so the result should be a 2x2 identity matrix.
comparisonTable = gi.AllIsomorphicQ('../graphs/cayley/3-cayley-permutations')

if rank==0: print comparisonTable