#!/usr/bin/env python2.7
# initialize PETSc
import sys
import os
import time
import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np

import sumatra
from sumatra.launch import DistributedLaunchMode
from sumatra.launch import SerialLaunchMode
from sumatra.parameters import build_parameters

# import pyCTQW as qw
import pyCTQW.MPI as qw

class capture(object):

    def __init__(self,np=1,reason=None,label=None):
        self.np = np
        self.reason = reason
        self.label = label + '_' + time.strftime("%d%m%y-%H%M%S")

        if self.np > 1:
            self.launchmode = DistributedLaunchMode(self.np)
            self.executable = sumatra.programs.Executable(path="mpirun")
        else:
            self.launchmode = SerialLaunchMode()
            self.executable = sumatra.programs.PythonExecutable(path=sys.executable)

    def __call__(self,main):

        def wrapped_main(parameters, *args, **kwargs):
            import sumatra.projects

            project = sumatra.projects.load_project()
            main_file = os.path.relpath(sys.modules['__main__'].__file__)

            record = project.new_record(parameters=parameters,
                                        main_file=main_file,
                                        executable=self.executable,
                                        launch_mode=self.launchmode,
                                        label=self.label,
                                        reason=self.reason)

            parameters.update({"sumatra_label": record.label})
            start_time = time.time()
            main(parameters, *args, **kwargs)
            record.duration = time.time() - start_time
            record.output_data = record.datastore.find_new_data(record.timestamp)
            project.add_record(record)
            project.save()

        return wrapped_main

@capture(np=4,label='1P_line')
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
    walk.plot('out/{}-plot.png'.format(parameters["sumatra_label"]))

    # plot the probability over time for the watched nodes
    walk.plotNodes('out/{}-nodes.png'.format(parameters["sumatra_label"]))

    # export final state
    walk.exportState("out/{}-state.txt".format(parameters["sumatra_label"]), "txt")

    # destroy the quantum walk
    walk.destroy()

parameter_file = sys.argv[1]
parameters = build_parameters(parameter_file)
main(parameters)