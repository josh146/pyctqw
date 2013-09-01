#!/usr/bin/env python2.7
# initialize PETSc
import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np

# import pyCTQW as qw
import pyCTQW.MPI as qw

# enable command line arguments -t and -N
OptDB = PETSc.Options()
N = OptDB.getInt('N', 20)
t = OptDB.getReal('t', 2)

# get the MPI rank
rank =  PETSc.Comm.Get_rank(PETSc.COMM_WORLD)

if rank == 0:
	print '3P Line\n'

# initialise an N (default 20) node graph CTQW
walk = qw.Line3P(N)

# Create a Hamiltonian with 2P very strong interaction.
walk.createH(interaction=10.)

# create the initial state |0,4,4>
init_state = [[0,4,4,1.0]]
walk.createInitState(init_state)

# set the eigensolver properties.
walk.EigSolver.setEigSolver(tol=1.e-2)
# underestimate the minimum eigenvalue
walk.EigSolver.setEigSolver(emin_estimate=0)

# Propagate the CTQW using the Chebyshev method for t=2s
walk.propagate(t,method='chebyshev')

# plot the marginal probabilities
# after propagation over all nodes
walk.plot('out/3p_line_plot.png')

# destroy the quantum walk
walk.destroy()