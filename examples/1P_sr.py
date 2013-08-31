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

if rank == 0:
	print '1P Strong Regular\n'

# initialise a 25 node graph CTQW, 
# and create a Hamiltonian from a file.
walk = qw.Graph(25)
walk.createH('../graphs/strong-regular-25-12-5-6/1.txt','txt',layout='circle')

# create the initial state |0,1>
init_state = [[0,1]]
walk.createInitState(init_state)

# Propagate the CTQW using Krylov subspace methods for t=100s
walk.propagate(100,method='krylov')

# plot the marginal probabilities
# after propagation over all nodes
walk.plot('out/1p_sr_plot.png')

# plot a 3D graph representation with bars
# representing the marginal probabilities
walk.plotGraph(output='out/1p_sr_graph.png')

# destroy the quantum walk
walk.destroy()