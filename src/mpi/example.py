#!/usr/bin/env python2.7
import sys, os, errno, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np

import pyctqw_MPI, pyctqw_plots

OptDB = PETSc.Options()
N = OptDB.getInt('N', 100)
t = OptDB.getInt('t', 20.0)
d = [3,4]
amp = [2.0,1.5]

rank =  PETSc.Comm.Get_rank(PETSc.COMM_WORLD)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------- Arbitrary CTQW --------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if rank == 0:	print '\ncycle CTQW\n'

#init_state = [[2,1.0]]

#walk = pyctqw_MPI.ctqwGraph(6)
#walk.createH('out/test.txt','txt',d=[0,5],amp=[100,100])
#walk.createInitState(init_state)

#walk.EigSolver.setEigSolver(tol=1.e-3)
#walk.propagate(t,method='expm')

#walk.plot('out/test-cycle.png')

#walk.destroy()

#sys.exit()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---------------------------------- 2P line ------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if rank == 0:	print '\n2p quantum walk\n'

init_state = [[0,1,1.0/np.sqrt(2.0)], [1,1,1.0j/np.sqrt(2.0)]]

walk = pyctqw_MPI.Line2P(N)
walk.createH(d,amp)
walk.createInitState(init_state)

walk.EigSolver.setEigSolver(tol=1.e-3)
walk.propagate(t,method='chebyshev')

walk.plot('out/testp2.png')
walk.exportState('out/test.txt','bin')

walk.destroy()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---------------------------------- 1P line ------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if rank == 0:	print '\n1p quantum walk\n'

#init_state = [[0.,1.0/np.sqrt(2.0)], [1.,1.0/np.sqrt(2.0)]]

#walk = pyctqw_MPI.Line(N)
#walk.createH(d,amp)
#walk.createInitState(init_state)

#walk.EigSolver.setEigSolver(tol=1.e-3)
#walk.propagate(t,method='chebyshev')

#walk.plot('out/testp1.png')

#walk.destroy()








