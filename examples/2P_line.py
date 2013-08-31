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
	print '2P Line\n'

# initialise an N (default 100) node graph CTQW
walk = qw.Line2P(N)

# Create a Hamiltonian with 2P interaction.
walk.createH(interaction=1.)

# create the initial state (1/sqrt(2)) (|0,0>-|1,1>)
init_state = [[0,0,1.0/np.sqrt(2.0)], [1,1,-1.0/np.sqrt(2.0)]]
walk.createInitState(init_state)

# set the eigensolver properties.
walk.EigSolver.setEigSolver(tol=1.e-2)
# underestimate the minimum eigenvalue
walk.EigSolver.setEigSolver(emin_estimate=0)

# create a handle to watch the probability at nodes -5,0,1:
walk.watch([0,1,-5])
# create a handler to watch the entanglement
walk.watch(None,watchtype='entanglement',verbose=False)

# Propagate the CTQW using the Chebyshev method
# for t=100s in timesteps of dt=0.1
for i in np.arange(0.1,t+0.1,0.1):
	walk.propagate(i,method='chebyshev')

# plot the marginal probabilities
# after propagation over all nodes
walk.plot('out/2p_line_plot.png')

# plot the probability over time for the watched nodes
walk.plotNodes('out/2p_line_nodes.png')

# plot the entanglement over time
walk.plotEntanglement('out/2p_line_ent.png')

# export the final state
walk.exportState('out/2p_line_state.bin','bin')

# destroy the quantum walk
walk.destroy()