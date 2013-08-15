#!/usr/bin/python
import sys, os, errno, time

from petsc4py import PETSc
from libpyctqw_MPI import ctqwmpi
import func
import numpy as np

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------ Eigensolver object (PETSc matrix input) --------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class EigSolver(object):
	def __init__(self,mat,**kwargs):		
		self.__default = {'esolver' : 'krylovschur',
			          'workType': 'null',
			          'workSize': '35',
			          'tol'     : 0.,
			          'maxIt'   : 0,
			          'verbose' : False,
			          'emax_estimate' : None,
			          'emin_estimate' : None,
			          'Emin_val': None,
			          'Emax_val': None}
			   
		self.__mat = mat
		self.__rank = PETSc.Comm.Get_rank(PETSc.COMM_WORLD)
			   
		for key,default in self.__default.iteritems():
			setattr(self, key, kwargs.get(key,default))

	def setEigSolver(self,**kwargs):
		for key in kwargs:
			if hasattr(self, key):
				setattr(self, key, kwargs.get(key))
			else:
				'Property type does not exist!'
			
	def getEigSolver(self,*args):
		for key in args:
			if hasattr(self, key):
				return getattr(self, key)
			else:
				'Property type does not exist!'

	def findEmax(self):
		# Calcaulate the eigenvalues
		if self.emax_estimate is not None:
			self.Emax_val = self.emax_estimate
			self.Emax_err = 0.

		else:
			EmaxStage = PETSc.Log.Stage('Emax'); EmaxStage.push()
			Emax,Emax_error,ierr = ctqwmpi.min_max_eigs(self.__mat.fortran,self.__rank,'max',
				self.esolver,self.workType,self.workSize,self.tol,self.maxIt,self.verbose)
			EmaxStage.pop()
			if ierr==0:
				self.Emax_val = Emax
				self.Emax_err = Emax_error

	def findEmin(self):
		# Calcaulate the eigenvalues
		
		if self.emin_estimate is not None:
			self.Emin_val = self.emin_estimate
			self.Emin_err = 0.

		else:
			EminStage = PETSc.Log.Stage('Emin'); EminStage.push()
			Emin,Emin_error,ierr = ctqwmpi.min_max_eigs(self.__mat.fortran,self.__rank,'min',
				self.esolver,self.workType,self.workSize,self.tol,self.maxIt,self.verbose)
			EminStage.pop()
			if ierr==0:
				self.Emin_val = Emin
				self.Emin_err = Emin_error
			
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------ Hamiltonian object (grid size  input) ----------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Hamiltonian(object):

	def __init__(self,N):
		self.rank = PETSc.Comm.Get_rank(PETSc.COMM_WORLD)
		self.N = N

		# define matrices
		self.mat = PETSc.Mat()
		self.mat.create(PETSc.COMM_WORLD)
		
		# create eigenvalue solver
		self.EigSolver = EigSolver(self.mat)
	
	def reinitialize(self):
		self.destroy()
		self.mat = PETSc.Mat()
		self.mat.create(PETSc.COMM_WORLD)
		
	def importAdj(self,filename,filetype,d=[0],amp=[0.],layout='spring'):
		try:
			if self.mat.isAssembled():
				self.reinitialize()
		except: pass
		# create the Hamiltonian
		Hamiltonian = PETSc.Log.Stage('Hamiltonian')
		Hamiltonian.push()
		
		self.Adj = func.loadMat(filename,filetype)
		self.mat = func.adjToH(self.Adj,d=d,amp=amp)
		
		Hamiltonian.pop()

		self.nodePos, self.lineX, self.lineY = func.getGraphNodes(self.Adj,layout=layout)

	def importAdjToH(self,filename,filetype,d=[0],amp=[0.],p='1',layout='spring'):
		try:
			if self.mat.isAssembled():
				self.reinitialize()
		except: pass
		# create the Hamiltonian
		Hamiltonian = PETSc.Log.Stage('Hamiltonian')
		Hamiltonian.push()
		
		ctqwmpi.importadjtoh(self.mat.fortran,filename,p,d,amp)
		
		Hamiltonian.pop()

		# create the adjacency matrix
		self.Adj = func.loadMat(filename,filetype)
		self.nodePos, self.lineX, self.lineY = func.getGraphNodes(self.Adj,layout=layout)
	
	def createLine2P(self,d=[0],amp=[0.]):
		try:
			if self.mat.isAssembled():
				self.reinitialize()
		except: pass
		self.defectNodes = d
		self.defectAmp = amp
		# create the Hamiltonian
		Hamiltonian = PETSc.Log.Stage('Hamiltonian')
		Hamiltonian.push()
		ctqwmpi.hamiltonian_2p_line(self.mat.fortran,self.defectNodes,self.defectAmp,self.N)
		Hamiltonian.pop()
	
	def createLine(self,d=[0],amp=[0.]):
		try:
			if self.mat.isAssembled():
				self.reinitialize()
		except: pass
		self.defectNodes = d
		self.defectAmp = amp
		# create the Hamiltonian
		Hamiltonian = PETSc.Log.Stage('Hamiltonian')
		Hamiltonian.push()
		ctqwmpi.hamiltonian_1p_line(self.mat.fortran,self.defectNodes,self.defectAmp,self.N)
		Hamiltonian.pop()
	
	def Emax(self,**kwargs):
		if self.EigSolver.Emax_val is None:
			self.EigSolver.findEmax()
		return self.EigSolver.Emax_val
	
	def Emin(self,**kwargs):
		if self.EigSolver.Emin_val is None:
			self.EigSolver.findEmin()
		return self.EigSolver.Emin_val
		
	def destroy(self):
		self.mat.destroy()
		try:
			self.Adj.destroy()
		except:
			pass

class nodeHandle(object):
	def __init__(self,nodes,t,psi):
		self.rank = PETSc.COMM_WORLD.Get_rank()
		self.time = [t]

		self.all_nodes = nodes

		self.local_nodes = []
		(Istart,Iend) = psi.getOwnershipRange()

		for i in self.all_nodes:
			if Istart <= i < Iend:
				self.local_nodes.append(i)

		self.local_prob = psi.getValues(self.local_nodes)

	def update(self,t,psi):
		self.time.append(t)

		prob_update = np.square(np.abs(psi.getValues(self.local_nodes)))
		self.local_prob = np.vstack([self.local_prob,prob_update])

	def getLocalNodes(self,t=None):
		if t is not None:
			try:
				indt = self.time.index(t)
			except ValueError:
				if self.rank == 0:	print '\nERROR: time {} was not handled'.format(t)
				return

			return self.local_prob[indt]

		else:
			return np.array(self.time), self.local_prob.T



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------- 1 particle CTQW   -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class QuantumWalkP1(object):
	def __init__(self,N):
		self.rank = PETSc.Comm.Get_rank(PETSc.COMM_WORLD)
		self.N = N
		self.t = 0
		
		# define vectors
		self.psi0 = PETSc.Vec()
		self.psi0.create(PETSc.COMM_WORLD)
		self.psi0.setSizes(self.N)
		self.psi0.setUp()
		self.psi = self.psi0.duplicate()
		self.prob = self.psi0.duplicate()
		
		# define matrices
		self.H = Hamiltonian(self.N)
		self.EigSolver = self.H.EigSolver
	
	def importInitState(self,filename,filetype):
		self.initState = 'file:'+filename
		# create the inital stage
		initStateS = PETSc.Log.Stage('initState')
		initStateS.push()
		try:
			self.psi0 = func.loadVec(filename,filetype)
			self.marginal(self.psi0)
		except:
			print '\nERROR: incorrect state (is it the correct length?'
			sys.exit()
		initStateS.pop()

	def marginal(self,vec):
		# calculate marginal probabilities
		Marginal = PETSc.Log.Stage('Marginal'); Marginal.push()
		ctqwmpi.p1prob(vec.fortran,self.prob.fortran,self.N)
		Marginal.pop()

	def watch(self,nodes,type='prob'):
		self.handle = nodeHandle(nodes,self.t,self.prob)
		
	def propagate(self,t,method='chebyshev',**kwargs):
		self.t = t
		self.EigSolver.setEigSolver(**kwargs)
		
		if method=='krylov':
			# SLEPc matrix exponential
			expmS = PETSc.Log.Stage('SLEPc expm'); expmS.push()
			ctqwmpi.expm(self.H.mat.fortran,self.t,self.psi0.fortran,self.psi.fortran)
			expmS.pop()
			
		elif method=='chebyshev':
			# Chebyshev algorithm
			chebyS = PETSc.Log.Stage('Chebyshev'); chebyS.push()
			ctqwmpi.qw_cheby(self.psi0.fortran,self.psi.fortran,self.t,self.H.mat.fortran,
					self.H.Emin(),self.H.Emax(),self.rank,self.N)
			chebyS.pop()
		
		self.marginal(self.psi)

		try:
			self.handle.update(self.t,self.prob)
		except:
			pass

	def exportState(self,filename,filetype):
		func.exportVec(self.psi,filename,filetype)

	def psiToInit(self):
		self.psi0 = self.psi
	
	def destroy(self):
		self.H.destroy()
		self.psi.destroy()
		self.psi0.destroy()
		self.prob.destroy()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------- 2 particle CTQW   -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class QuantumWalkP2(object):
	def __init__(self,N):
		self.rank = PETSc.Comm.Get_rank(PETSc.COMM_WORLD)
		self.N = N
		self.t = 0
		
		# define vectors
		self.psi0 = PETSc.Vec()
		self.psi0.create(PETSc.COMM_WORLD)
		self.psi0.setSizes(self.N**2)
		self.psi0.setUp()
		self.psi = self.psi0.duplicate()
		
		# create additional vectors
		self.psiX = PETSc.Vec()
		self.psiX.create(PETSc.COMM_WORLD)
		self.psiX.setSizes(self.N)
		self.psiX.setUp()
		self.psiY = self.psiX.duplicate()
		
		# define matrices
		self.H = Hamiltonian(self.N)
		self.EigSolver = self.H.EigSolver
	
	def importInitState(self,filename,filetype):
		self.initState = 'file:'+filename
		# create the inital stage
		initStateS = PETSc.Log.Stage('initState')
		initStateS.push()
		try:
			if filetype == 'txt':
				self.psi0 = func.loadMatToVec(filename,filetype)
			elif filetype == 'bin':
				self.psi0 = func.loadVec(filename,filetype)
		except:
			print '\nERROR: incorrect state (is it the correct length?)'
			sys.exit()
		initStateS.pop()

	def marginal(self):
		# calculate marginal probabilities
		Marginal = PETSc.Log.Stage('Marginal'); Marginal.push()
		ctqwmpi.marginal(self.psi.fortran,self.psiX.fortran,'x',self.N)
		ctqwmpi.marginal(self.psi.fortran,self.psiY.fortran,'y',self.N)
		Marginal.pop()
		
	def propagate(self,t,method='expm',**kwargs):
		self.t = t
		self.EigSolver.setEigSolver(**kwargs)
		
		if method=='krylov':
			# SLEPc matrix exponential
			krylov = PETSc.Log.Stage('SLEPc krylov'); krylov.push()
			ctqwmpi.expm(self.H.mat.fortran,self.t,self.psi0.fortran,self.psi.fortran)
			krylov.pop()
			
		elif method=='chebyshev':
			# Chebyshev algorithm
			chebyS = PETSc.Log.Stage('Chebyshev'); chebyS.push()
			ctqwmpi.qw_cheby(self.psi0.fortran,self.psi.fortran,self.t,self.H.mat.fortran,
					self.H.Emin(),self.H.Emax(),self.rank,self.N)
			chebyS.pop()
		
		self.marginal()
	
	def exportState(self,filename,filetype):
		if filetype == 'txt':
			func.exportVecToMat(self.psi,filename,filetype)
		elif filetype == 'bin':
			func.exportVec(self.psi,filename,filetype)

	def psiToInit(self):
		self.psi0 = self.psi
	
	def destroy(self):
		self.H.destroy()
		self.psi.destroy()
		self.psi0.destroy()
		self.psiX.destroy()
		self.psiY.destroy()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------- Arbitrary CTQW --------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class ctqwGraph(QuantumWalkP1):
	def __init__(self,N,filename=None,filetype=None,d=None,amp=None):
		QuantumWalkP1.__init__(self,N)
		self.liveplot = False
		
		if (filename and filetype) is not None:
			self.createH(filename,filetype,d=d,amp=amp)
		
	def createH(self,filename,filetype,d=None,amp=None,layout='spring'):
		if (d and amp) is None:
			self.defectNodes = [0]
			self.defectAmp = [0.]
		else:
			self.defectNodes = d
			self.defectAmp = amp
		self.H.importAdj(filename,filetype,d=self.defectNodes,amp=self.defectAmp,layout=layout)
		
	def createInitState(self,initState):
		self.initState = np.vstack([np.array(initState).T[0]-self.N/2+1,
			   	 	    np.array(initState).T[1]]).T.tolist()
	
		# create the inital stage
		initStateS = PETSc.Log.Stage('initState')
		initStateS.push()
		ctqwmpi.p1_init(self.psi0.fortran,self.initState,self.N)
		initStateS.pop()

		self.marginal(self.psi0)
		
	def plot(self,filename):
		if os.path.isabs(filename):
			outDir = os.path.dirname(filename)
		else:
			outDir = './'+os.path.dirname(filename)
	
		# create output directory if it doesn't exist
		try:
			os.mkdir(outDir)
		except OSError as exception:
			if exception.errno != errno.EEXIST:
				raise

		initstateLabels = []
		for i in range(len(self.initState)):
			initstateLabels.append([sum(pair).real for pair in zip(self.initState[i], [self.N/2-1,0])])

		plotStage = PETSc.Log.Stage('Plotting'); plotStage.push()		
		func.plot(np.arange(self.N),self.prob,filename,self.t,initstateLabels,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()

	def plotGraph(self,size=(12,8),probX=True,output=None,**kwargs):

		rank = PETSc.COMM_WORLD.Get_rank()

		if probX:
			# scatter prob to process 0
			commX = self.prob.getComm()
			rank = PETSc.COMM_WORLD.getRank()
			scatterX, probX0 = PETSc.Scatter.toZero(self.prob)
			scatterX.scatter(self.prob, probX0, False, PETSc.Scatter.Mode.FORWARD)

		if rank==0:
			from matplotlib import pyplot as plt
			import mpl_toolkits.mplot3d as plt3d
			self.fig = plt.figure(figsize=size)
			self.ax = plt3d.Axes3D(self.fig)
			self.ax.view_init(45, -50)
			self.ax.set_axis_off()

			if probX:
				prob = np.real(np.asarray(probX0))
			else:
				prob = None

			func.plotGraph(self.ax,self.H.nodePos,self.H.lineX,self.H.lineY,
				prob=prob,**kwargs)
		
			self.ax.set_title('$t={}$'.format(self.t))

			if type(output) is str:
				plt.savefig(output)
			else:
				plt.show(block=True)
				plt.close()

		if probX:
			# deallocate	
			commX.barrier()
			scatterX.destroy()
			probX0.destroy()

	def plotLiveGraph(self,dt,size=(12,8),**kwargs):

		rank = PETSc.COMM_WORLD.Get_rank()

		# scatter prob to process 0
		commX = self.prob.getComm()
		rank = PETSc.COMM_WORLD.getRank()
		scatterX, probX0 = PETSc.Scatter.toZero(self.prob)
		scatterX.scatter(self.prob, probX0, False, PETSc.Scatter.Mode.FORWARD)

		if rank==0:
			from matplotlib import pyplot as plt
			import mpl_toolkits.mplot3d as plt3d
			if not self.liveplot:
				plt.ion()
				self.figLive = plt.figure(figsize=size)
				self.liveplot = True
			else:
				plt.cla()

			self.axLive = plt3d.Axes3D(self.figLive)
			self.axLive.view_init(45, -50)
			self.axLive.set_axis_off()

			prob = np.real(np.asarray(probX0))

			func.plotGraph(self.axLive,self.H.nodePos,self.H.lineX,self.H.lineY,
				prob=prob,**kwargs)
		
			self.axLive.set_title('$t={}$'.format(self.t))

			plt.draw()#show(block=True)
			time.sleep(dt)

		# deallocate	
		commX.barrier()
		scatterX.destroy()
		probX0.destroy()

	def clearLiveGraph(self):
		rank = PETSc.COMM_WORLD.Get_rank()

		if rank == 0:
			from matplotlib import pyplot as plt
			plt.show(block=True)
			self.liveplot = False
			plt.close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------- 2P Arbitrary CTQW -----------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class ctqwGraph2P(QuantumWalkP2):
	def __init__(self,N,filename=None,filetype=None,d=None,amp=None):
		QuantumWalkP2.__init__(self,N)
		self.liveplot = False
		
		if (filename and filetype) is not None:
			self.createH(filename,filetype,d=d,amp=amp)
		
	def createH(self,filename,filetype,d=None,amp=None,layout='spring'):
		if (d and amp) is None:
			self.defectNodes = [0]
			self.defectAmp = [0.]
		else:
			self.defectNodes = d
			self.defectAmp = amp

		self.H.importAdjToH(filename,filetype,
			d=self.defectNodes,amp=self.defectAmp,p='2',layout=layout)

		# ctqwmpi.importadjtoh(self.H.mat.fortran,filename,'2',
		# 	d=self.defectNodes,amp=self.defectAmp,layout=layout)
		
	def createInitState(self,initState):
		self.initState = np.vstack([np.array(initState).T[0]-self.N/2+1,
			   	 	    np.array(initState).T[1]-self.N/2+1,np.array(initState).T[2]]).T.tolist()
	
		# create the inital stage
		initStateS = PETSc.Log.Stage('initState')
		initStateS.push()
		ctqwmpi.p2_init(self.psi0.fortran,self.initState,self.N)
		initStateS.pop()
		
	def plot(self,filename):
		if os.path.isabs(filename):
			outDir = os.path.dirname(filename)
		else:
			outDir = './'+os.path.dirname(filename)
	
		# create output directory if it doesn't exist
		try:
			os.mkdir(outDir)
		except OSError as exception:
			if exception.errno != errno.EEXIST:
				raise

		initstateLabels = []
		for i in range(len(self.initState)):
			initstateLabels.append([sum(pair).real for pair in zip(self.initState[i], [self.N/2-1,self.N/2-1,0])])

		plotStage = PETSc.Log.Stage('Plotting'); plotStage.push()		
		func.plot2P(np.arange(self.N),self.psiX,self.psiY,filename,self.t,initstateLabels,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()

	def plotGraph(self,size=(12,8),probX=True,probY=True,output=None,**kwargs):

		rank = PETSc.COMM_WORLD.Get_rank()

		if probX:
			# scatter prob to process 0
			commX = self.psiX.getComm()
			rank = PETSc.COMM_WORLD.getRank()
			scatterX, probX0 = PETSc.Scatter.toZero(self.psiX)
			scatterX.scatter(self.psiX, probX0, False, PETSc.Scatter.Mode.FORWARD)

		if probY:
			# scatter prob to process 0
			commY = self.psiY.getComm()
			rank = PETSc.COMM_WORLD.getRank()
			scatterY, probY0 = PETSc.Scatter.toZero(self.psiY)
			scatterX.scatter(self.psiY, probY0, False, PETSc.Scatter.Mode.FORWARD)

		if rank==0:
			from matplotlib import pyplot as plt
			import mpl_toolkits.mplot3d as plt3d

			self.fig = plt.figure(figsize=size)
			self.ax = plt3d.Axes3D(self.fig)
			self.ax.view_init(45, -50)
			self.ax.set_axis_off()

			if probX:
				prob = np.real(np.asarray(probX0))
			else:
				prob = None

			if probY:
				prob2 = np.real(np.asarray(probY0))
			else:
				prob2 = None

			func.plotGraph(self.ax,self.H.nodePos,self.H.lineX,self.H.lineY,
				prob=prob,prob2=prob2,**kwargs)
		
			self.ax.set_title('$t={}$'.format(self.t))

			if type(output) is str:
				plt.savefig(output)
			else:
				plt.show(block=True)

		if probX:
			# deallocate	
			commX.barrier()
			scatterX.destroy()
			probX0.destroy()

		if probY:
			# deallocate	
			commY.barrier()
			scatterY.destroy()
			probY0.destroy()

	def plotLiveGraph(self,dt,size=(12,8),
				probX=True,probY=True,**kwargs):

		rank = PETSc.COMM_WORLD.Get_rank()

		if probX:
			# scatter prob to process 0
			commX = self.psiX.getComm()
			rank = PETSc.COMM_WORLD.getRank()
			scatterX, probX0 = PETSc.Scatter.toZero(self.psiX)
			scatterX.scatter(self.psiX, probX0, False, PETSc.Scatter.Mode.FORWARD)

		if probY:
			# scatter prob to process 0
			commY = self.psiY.getComm()
			rank = PETSc.COMM_WORLD.getRank()
			scatterY, probY0 = PETSc.Scatter.toZero(self.psiY)
			scatterX.scatter(self.psiY, probY0, False, PETSc.Scatter.Mode.FORWARD)

		if rank==0:
			from matplotlib import pyplot as plt
			import mpl_toolkits.mplot3d as plt3d

			if not self.liveplot:
				plt.ion()
				self.figLive = plt.figure(figsize=size)
				self.liveplot = True
			else:
				plt.cla()

			self.axLive = plt3d.Axes3D(self.figLive)
			self.axLive.view_init(45, -50)
			self.axLive.set_axis_off()

			if probX:
				prob = np.real(np.asarray(probX0))
			else:
				prob = None

			if probY:
				prob2 = np.real(np.asarray(probY0))
			else:
				prob2 = None

			func.plotGraph(self.axLive,self.H.nodePos,self.H.lineX,self.H.lineY,
				prob=prob,prob2=prob2,**kwargs)
		
			self.axLive.set_title('$t={}$'.format(self.t))

			plt.draw()#show(block=True)
			time.sleep(dt)

		if probX:
			# deallocate	
			commX.barrier()
			scatterX.destroy()
			probX0.destroy()

		if probY:
			# deallocate	
			commY.barrier()
			scatterY.destroy()
			probY0.destroy()

	def clearLiveGraph(self):
		rank = PETSc.COMM_WORLD.Get_rank()

		if rank == 0:
			from matplotlib import pyplot as plt
			plt.show(block=True)
			self.liveplot = False
			plt.close()
		
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------- 1 particle CTQW on a line -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Line(QuantumWalkP1):
	def __init__(self,N,d=None,amp=None):
		QuantumWalkP1.__init__(self,N)
		if (d is not None) and (amp is not None):
			self.defectNodes = d
			self.defectAmp = amp
			self.H.createLine(d,amp)
		
	def createH(self,d=None,amp=None):
		if (d and amp) is None:
			self.defectNodes = [0]
			self.defectAmp = [0.]
		else:
			self.defectNodes = d
			self.defectAmp = amp
		self.H.createLine(self.defectNodes,self.defectAmp)
		
	def createInitState(self,initState):
		self.initState = initState
		# create the inital stage
		initStateS = PETSc.Log.Stage('initState')
		initStateS.push()
		ctqwmpi.p1_init(self.psi0.fortran,initState,self.N)
		initStateS.pop()
		
	def plot(self,filename):
		if os.path.isabs(filename):
			outDir = os.path.dirname(filename)
		else:
			outDir = './'+os.path.dirname(filename)
	
		# create output directory if it doesn't exist
		try:
			os.mkdir(outDir)
		except OSError as exception:
			if exception.errno != errno.EEXIST:
				raise

		plotStage = PETSc.Log.Stage('Plotting'); plotStage.push()		
		func.plot(np.arange(1-self.N/2,self.N/2+1),self.prob,filename,self.t,self.initState,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()
			
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------- 2 particle CTQW on a line -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Line2P(QuantumWalkP2):
	def __init__(self,N,d=None,amp=None):
		QuantumWalkP2.__init__(self,N)

		if (d is not None) and (amp is not None):
			self.defectNodes = d
			self.defectAmp = amp
			self.H.createLine2P(d,amp)
		
	def createH(self,d=None,amp=None):
		if (d and amp) is None:
			self.defectNodes = [0]
			self.defectAmp = [0.]
		else:
			self.defectNodes = d
			self.defectAmp = amp
		self.H.createLine2P(self.defectNodes,self.defectAmp)
		
	def createInitState(self,initState):
		self.initState = initState
		# create the inital stage
		initStateS = PETSc.Log.Stage('initState')
		initStateS.push()
		ctqwmpi.p2_init(self.psi0.fortran,initState,self.N)
		initStateS.pop()
		
	def plot(self,filename):
		if os.path.isabs(filename):
			outDir = os.path.dirname(filename)
		else:
			outDir = './'+os.path.dirname(filename)
	
		# create output directory if it doesn't exist
		try:
			os.mkdir(outDir)
		except OSError as exception:
			if exception.errno != errno.EEXIST:
				raise

		plotStage = PETSc.Log.Stage('Plotting'); plotStage.push()		
		func.plot2P(np.arange(1-self.N/2,self.N/2+1),self.psiX,self.psiY,filename,self.t,self.initState,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()

