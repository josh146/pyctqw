#!/usr/bin/python
import sys, os, errno, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np

import pyctqw_MPI

OptDB = PETSc.Options()
N = OptDB.getInt('N', 100)
t = OptDB.getInt('t', 20.0)
d = [3,4]
amp = [2.0,1.5]
init_state = [[0.,1.,1.0/np.sqrt(2.0)], [1.,1.,1.0/np.sqrt(2.0)]]

walk = pyctqw_MPI.Line2P(N)
walk.createH(d,amp)
walk.createInitState(init_state)
walk.propagate(t,method='chebyshev')
walk.psiX.view()
walk.psiY.view()
walk.plot('./test.png')
walk.destroy()
