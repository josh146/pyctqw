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
#------------------------------- Strong Regular --------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if rank == 0:	print '1P Strong Regular\n'

init_state = [[0,1]]

walk = qw.Graph(25)
walk.createH('../graphs/strong-regular-25-12-5-6/1.txt','txt',layout='circle')

walk.createInitState(init_state)

walk.EigSolver.setEigSolver(tol=1.e-2,verbose=False,emin_estimate=0.,emax_estimate=15.)

walk.propagate(100.,method='krylov')

walk.plot('out/reg.png')
walk.plotGraph(output='out/reg-graph.png')

walk.destroy()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 2P Strong Regular -----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if rank == 0:	print '2P Strong Regular\n'

init_state = [[0,1,1.0/np.sqrt(2.0)], [1,1,1.0j/np.sqrt(2.0)]]

walk = qw.Graph2P(25)
walk.createH('../graphs/strong-regular-25-12-5-6/1.txt','txt',layout='circle',d=[0],amp=[0.])
#qw.func.exportMat(walk.H.Adj,'../graphs/strong-regular-25-12-5-6/1p.txt','txt',mattype='adj')

walk.createInitState(init_state)

walk.EigSolver.setEigSolver(tol=1.e-2,verbose=False,emin_estimate=0.,emax_estimate=30.)

walk.propagate(100,method='krylov')

walk.plot('out/reg-2p.png')
walk.plotGraph(output='out/reg-graph-2p.png')

walk.destroy()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------- 1P 3-Caley Tree CTQW --------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if rank == 0:	print '1P 3-Caley Tree CTQW\n'

init_state = [[0,1.0/np.sqrt(2.0)], [1,1.0j/np.sqrt(2.0)]]

walk = qw.Graph(10)
walk.createH('../graphs/3-caley.txt','txt',d=[0],amp=[3.])

walk.createInitState(init_state)

walk.EigSolver.setEigSolver(tol=1.e-2,verbose=False,emin_estimate=0.,emax_estimate=10.)

walk.watch([0,1,2,3,4])

for t2 in np.arange(0.01,t+0.01,0.01):
	walk.propagate(t2,method='chebyshev')

walk.plot('out/3-caley-1p.png')
#walk.plotGraph(nodetextbg='blue')
walk.plotGraph(output='out/3-caley-1p-graph.png')
walk.plotNodes('out/3-caley-1p-nodes.png')

walk.destroy()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------- 2P 3-Caley Tree CTQW ------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if rank == 0:	print '2P 3-Caley Tree CTQW\n'

init_state = [[0,1,1.0/np.sqrt(2.0)], [1,1,1.0j/np.sqrt(2.0)]]

walk = qw.Graph2P(10)
walk.createH('../graphs/3-caley.txt','txt',d=d,amp=amp,layout='spring',interaction=0.5)

qw.func.exportMat(walk.H.mat,'out/3-caley-2p-hamiltonian.txt','txt')
walk.createInitState(init_state)

walk.EigSolver.setEigSolver(tol=1.e-2,verbose=False,emin_estimate=0.)

walk.watch([0,1,2,3,4,9])

for t2 in np.arange(0.01,5+0.01,0.01):
	walk.propagate(t2,method='chebyshev')
#	walk.plotLiveGraph(0.2)

#walk.clearLiveGraph()

walk.plot('out/3-caley-2p.png')
walk.plotGraph(output='out/3-caley-2p-graph.png')
walk.plotNode('out/3-caley-2p-node1.png',1)
walk.plotNodes('out/3-caley-2p-nodes-particle1.png',p=1)
walk.plotNodes('out/3-caley-2p-nodes-particle2.png',p=2)


walk.destroy()


# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #--------------------------- 3P 3-Caley Tree CTQW ------------------------------
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if rank == 0:	print '3P 3-Caley Tree CTQW\n'

walk = qw.Graph3P(10)
walk.createH('../graphs/3-caley.txt','txt',d=d,amp=amp,layout='spring',interaction=1.5)

init_state = [[0,1,4,1.0/np.sqrt(2.0)], [1,1,4,1.j/np.sqrt(2.0)]]
walk.createInitState(init_state)

walk.EigSolver.setEigSolver(tol=1.e-2,verbose=False,emin_estimate=0.)

walk.watch([0,1,2,3,4,9])

for t2 in np.arange(0.01,5+0.01,0.01):
	walk.propagate(t2,method='chebyshev')

walk.plot('out/3-caley-3p.png')
walk.plotNode('out/3-caley-3p-node1.png',1)
walk.plotNodes('out/3-caley-3p-nodes-particle1.png',p=1)
walk.plotNodes('out/3-caley-3p-nodes-particle2.png',p=2)
walk.plotNodes('out/3-caley-3p-nodes-particle3.png',p=3)


walk.destroy()

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #---------------------------------- 1P line ------------------------------------
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if rank == 0:	print '1P Line\n'

init_state = [[0.,1.0/np.sqrt(2.0)], [1.,1.0/np.sqrt(2.0)]]

walk = qw.Line(N)
walk.createH(d,amp)
walk.createInitState(init_state)

walk.EigSolver.setEigSolver(tol=1.e-3)

walk.watch([0,1,-5])

for i in range(int(t/0.01)):
	walk.propagate(0.01,method='chebyshev')
	walk.psiToInit()

walk.plot('out/line1p.png')
walk.plotNodes('out/line1p-nodes.png')

walk.destroy()

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #---------------------------------- 2P line ------------------------------------
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if rank == 0:	print '2P Line\n'

init_state = [[0,1,1.0/np.sqrt(2.0)], [1,1,1.0j/np.sqrt(2.0)]]

walk = qw.Line2P(N)
walk.createH(interaction=10.)
walk.createInitState(init_state)

walk.EigSolver.setEigSolver(tol=1.e-2)
walk.EigSolver.setEigSolver(emin_estimate=0)

walk.propagate(t,method='chebyshev')
walk.plot('out/line2.png')

walk.exportState('out/line2p-state.bin','bin')

walk.destroy()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---------------------------------- 3P line ------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if rank == 0:	print '3P Line\n'

init_state = [[0,4,4,1.0]]

walk = qw.Line3P(20)
walk.createH(interaction=10.)
walk.createInitState(init_state)

qw.func.exportVec( walk.psi0, 'text.txt', 'txt')

walk.EigSolver.setEigSolver(tol=1.e-2)
walk.EigSolver.setEigSolver(emin_estimate=0)

walk.propagate(2.,method='chebyshev')
walk.plot('out/line3.png')

walk.destroy()