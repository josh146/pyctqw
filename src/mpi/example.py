#!/usr/bin/env python2.7
import sys, os, errno, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np

import pyctqw_MPI, pyctqw_plots

OptDB = PETSc.Options()
N = OptDB.getInt('N', 100)
t = OptDB.getReal('t', 20)
d = [3,4]
amp = [2.0,1.5]

rank =  PETSc.Comm.Get_rank(PETSc.COMM_WORLD)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------- Arbitrary CTQW --------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if rank == 0:	print '3-Caley Tree CTQW\n'

init_state = [[0,1,1.0/np.sqrt(2.0)], [1,1,1.0j/np.sqrt(2.0)]]

walk = pyctqw_MPI.ctqwGraph2P(10)
walk.createH('3-caley.txt','txt',d=[0,0],amp=[0,0])

pyctqw_plots.exportMat(walk.H.mat,'3-caley-2p.txt','txt')
walk.createInitState(init_state)

walk.EigSolver.setEigSolver(tol=1.e-2,verbose=True,emin_estimate=0.)
walk.propagate(t,method='chebyshev')

walk.plot('out/3-caley.png')


walk.destroy()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---------------------------------- 2P line ------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if rank == 0:	print '\n2p quantum walk\n'

init_state = [[0,1,1.0/np.sqrt(2.0)], [1,1,1.0j/np.sqrt(2.0)]]

walk = pyctqw_MPI.Line2P(N)
walk.createH(d,amp)
walk.createInitState(init_state)

walk.EigSolver.setEigSolver(tol=1.e-2)
walk.EigSolver.setEigSolver(emin_estimate=0)

for t in range(1,21):
	walk.propagate(1,method='chebyshev')
	walk.plot('out/'+str(t)+'.png')
	walk.psiToInit()

#walk.exportState('out/test.txt','bin')

walk.destroy()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---------------------------------- 1P line ------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if rank == 0:	print '\n1p quantum walk\n'

init_state = [[0.,1.0/np.sqrt(2.0)], [1.,1.0/np.sqrt(2.0)]]

walk = pyctqw_MPI.Line(N)
walk.createH(d,amp)
walk.createInitState(init_state)

walk.EigSolver.setEigSolver(tol=1.e-3)
walk.propagate(t,method='chebyshev')

walk.plot('out/testp1.png')

walk.destroy()








