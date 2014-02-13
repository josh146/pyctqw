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
N = OptDB.getInt('N', 100)
t = OptDB.getReal('t', 20)

# get the MPI rank
rank =  PETSc.Comm.Get_rank(PETSc.COMM_WORLD)

if rank == 0:
	print '1P Line\n'

# initialise an N (default 100) node graph CTQW
walk = qw.Line(N)

# Create a Hamiltonian with defect and amplitude as below.
d = [3,4]
amp = [2.0,1.5]
walk.createH(d,amp)

# create the initial state (1/sqrt(2)) (|0>+|1>)
init_state = [[0.,1.0/np.sqrt(2.0)], [1.,1.0/np.sqrt(2.0)]]
walk.createInitState(init_state)

# set the eigensolver properties.
walk.EigSolver.setEigSolver(tol=1.e-3)

# create a handle to watch the probability at nodes -5,0,1:
walk.watch([0,1,-5])

# Propagate the CTQW using the Chebyshev method
# for t=100s in timesteps of dt=0.01
# Note that psiToInit() is being used rather than global timesteps.
for i in range(int(t/0.01)):
	walk.propagate(0.01,method='chebyshev')
	walk.psiToInit()

# plot the marginal probabilities
# after propagation over all nodes
walk.plot('out/1p_line_plot.png')

# plot the probability over time for the watched nodes
walk.plotNodes('out/1p_line_nodes.png')

# export final state
walk.exportState("out/1p_final_state.txt", "txt")

# destroy the quantum walk
walk.destroy()