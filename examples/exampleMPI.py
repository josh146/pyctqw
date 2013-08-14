#!/usr/bin/env python2.7
import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np

import pyCTQW.MPI as qw

OptDB = PETSc.Options()
N = OptDB.getInt('N', 100)
t = OptDB.getReal('t', 20)
d = [3,4]
amp = [2.0,1.5]

rank =  PETSc.Comm.Get_rank(PETSc.COMM_WORLD)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------- Arbitrary CTQW --------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if rank == 0:	print '1P 3-Caley Tree CTQW\n'

init_state = [[0,1.0/np.sqrt(2.0)], [1,1.0j/np.sqrt(2.0)]]

walk = qw.ctqwGraph(10)
walk.createH('../graphs/3-caley.txt','txt',d=[0,0],amp=[0,0])

walk.createInitState(init_state)

walk.EigSolver.setEigSolver(tol=1.e-2,verbose=True,emin_estimate=0.)

walk.propagate(t,method='krylov')

walk.plot('out/3-caley-p1.png')
walk.plotGraph(nodetextbg='blue')
walk.plotGraph(output='out/3-caley-graph-p1.png')


walk.destroy()
sys.exit()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------- Arbitrary CTQW --------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if rank == 0:	print '2P 3-Caley Tree CTQW\n'

init_state = [[0,1,1.0/np.sqrt(2.0)], [1,1,1.0j/np.sqrt(2.0)]]

walk = qw.ctqwGraph2P(10)
walk.createH('../graphs/3-caley.txt','txt',d=[0,0],amp=[0,0],layout='spring')

#qw.func.exportMat(walk.H.mat,'3-caley-2p.txt','txt')
walk.createInitState(init_state)

walk.EigSolver.setEigSolver(tol=1.e-2,verbose=True,emin_estimate=0.)

for t2 in range(1,6):
	walk.propagate(t2,method='chebyshev')
	walk.plotLiveGraph(0.2)
	walk.psiToInit()

walk.clearLiveGraph()

walk.plot('out/3-caley.png')
walk.plotGraph(output='out/plot.png')


walk.destroy()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---------------------------------- 2P line ------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if rank == 0:	print '\n2p quantum walk\n'

init_state = [[0,1,1.0/np.sqrt(2.0)], [1,1,1.0j/np.sqrt(2.0)]]

walk = qw.Line2P(N)
walk.createH(d,amp)
walk.createInitState(init_state)

walk.EigSolver.setEigSolver(tol=1.e-2)
walk.EigSolver.setEigSolver(emin_estimate=0)

for t2 in range(1,21):
	walk.propagate(1,method='chebyshev')
	walk.plot('out/'+str(t2)+'.png')
	walk.psiToInit()

#walk.exportState('out/test.txt','bin')

walk.destroy()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---------------------------------- 1P line ------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if rank == 0:	print '\n1p quantum walk\n'

init_state = [[0.,1.0/np.sqrt(2.0)], [1.,1.0/np.sqrt(2.0)]]

walk = qw.Line(N)
walk.createH(d,amp)
walk.createInitState(init_state)

walk.EigSolver.setEigSolver(tol=1.e-3)
walk.propagate(t,method='chebyshev')

walk.plot('out/testp1.png')

walk.destroy()








