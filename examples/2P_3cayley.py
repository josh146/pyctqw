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
N = OptDB.getInt('N', 10)
t = OptDB.getReal('t', 5)

# get the MPI rank
rank =  PETSc.Comm.Get_rank(PETSc.COMM_WORLD)

if rank == 0:
	print '2P 3-cayley Tree CTQW\n'

# initialise a 10 node graph CTQW
walk = qw.Graph2P(N)

# Create a 2 particle interacting Hamiltonian from a file.
# A defect of amplitude 2 and 0.5 has
# been added to nodes 3 and 4 respectively.
d = [3,4]
amp = [2.0,1.5]
walk.createH('../graphs/cayley/3-cayley.txt','txt',d=d,amp=amp,layout='spring',interaction=0.5)

# create the initial state (1/sqrt(2)) (|0,1>+i|1,1>)
init_state = [[0,1,1.0/np.sqrt(2.0)], [1,1,1.0j/np.sqrt(2.0)]]
walk.createInitState(init_state)

# set the eigensolver properties.
#Note that Emin has been set exactly.
walk.EigSolver.setEigSolver(tol=1.e-2,verbose=False,emin_estimate=0.)

# create a handle to watch the probability at nodes 0-4,9:
walk.watch([0,1,2,3,4,9])
# create a handle to watch the entanglement
walk.watch(None,watchtype='entanglement',verbose=False,esolver='lapack',max_it=100)

# Propagate the CTQW using the Chebyshev method
# for t=5s in timesteps of dt=0.01
for dt in np.arange(0.01,t+0.01,0.01):
	walk.propagate(dt,method='chebyshev')

# plot the p1 and p2 marginal probabilities
# after propagation over all nodes
walk.plot('out/2p_3cayley_plot.png')

# plot a 3D graph representation with bars
# representing the p1 and p2 marginal probabilities
walk.plotGraph(output='out/2p_3cayley_graph.png')

# plot the p1 and p2 probability over time for watched node1
walk.plotNode('out/2p_3cayley_node1.png',1)
# plot the particle 1 probability over all watched nodes
walk.plotNodes('out/2p_3cayley_nodes_particle1.png',p=1)
# plot the particle 2 probability over all watched nodes
walk.plotNodes('out/2p_3cayley_nodes_particle2.png',p=2)

# plot the entanglement vs. time
walk.plotEntanglement('out/2p_3cayley_ent.png')

# export the partial trace
walk.exportPartialTrace('out/3-cayley-2p-rhoX.txt','txt',p=1)

# destroy the quantum walk
walk.destroy()