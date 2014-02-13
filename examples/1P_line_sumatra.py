#!/usr/bin/env python2.7
# initialize PETSc
import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np

from sumatra.parameters import build_parameters
from sumatra.decorators import capture

# import pyCTQW as qw
import pyCTQW.MPI as qw

__file__ = "pyctqw/trunk/examples/1P_line_sumatra.py"

@capture
def main(parameters):
	# get the MPI rank
	rank = PETSc.Comm.Get_rank(PETSc.COMM_WORLD)

	N = parameters["vertices"]
	t = parameters["propagation_time"]

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
	walk.plot('out/{}plot.png'.format(parameters["sumatra_label"]))

	# plot the probability over time for the watched nodes
	walk.plotNodes('out/{}nodes.png'.format(parameters["sumatra_label"]))

	# export final state
	walk.exportState("out/{}state.txt".format(parameters["sumatra_label"]), "txt")

	# destroy the quantum walk
	walk.destroy()

parameter_file = sys.argv[1]
parameters = build_parameters(parameter_file)
main(parameters)