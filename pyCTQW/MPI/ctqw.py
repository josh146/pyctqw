#!/usr/bin/python
import os as _os
import errno as _errno
import sys as _sys
import time as _time

from petsc4py import PETSc as _PETSc
from libpyctqw_MPI import ctqwmpi
import func
import numpy as _np

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------ Eigensolver object (_PETSc matrix i_nput) --------------------
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
		self.__rank = _PETSc.Comm.Get_rank(_PETSc.COMM_WORLD)
			   
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
			EmaxStage = _PETSc.Log.Stage('Emax'); EmaxStage.push()
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
			EminStage = _PETSc.Log.Stage('Emin'); EminStage.push()
			Emin,Emin_error,ierr = ctqwmpi.min_max_eigs(self.__mat.fortran,self.__rank,'min',
				self.esolver,self.workType,self.workSize,self.tol,self.maxIt,self.verbose)
			EminStage.pop()
			if ierr==0:
				self.Emin_val = Emin
				self.Emin_err = Emin_error
			
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------ Hamiltonian object (grid size  i_nput) ----------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Hamiltonian(object):

	def __init__(self,N):
		self.rank = _PETSc.Comm.Get_rank(_PETSc.COMM_WORLD)
		self.N = N

		# define matrices
		self.mat = _PETSc.Mat()
		self.mat.create(_PETSc.COMM_WORLD)
		
		# create eigenvalue solver
		self.EigSolver = EigSolver(self.mat)
	
	def reinitialize(self):
		self.destroy()
		self.mat = _PETSc.Mat()
		self.mat.create(_PETSc.COMM_WORLD)
		
	def importAdj(self,filename,filetype,d=[0],amp=[0.],layout='spring',delimiter=None):
		try:
			if self.mat.isAssembled():
				self.reinitialize()
		except: pass
		# create the Hamiltonian
		Hamiltonian = _PETSc.Log.Stage('Hamiltonian')
		Hamiltonian.push()
		
		self.Adj = func.loadMat(filename,filetype,delimiter=delimiter)
		self.mat = func.adjToH(self.Adj,d=d,amp=amp)
		
		Hamiltonian.pop()

		self.nodePos, self.lineX, self.lineY = func.getGraphNodes(self.Adj,layout=layout)

	def importAdjToH(self,filename,filetype,d=[0.],amp=[0.],p='1',layout='spring',delimiter=None,interaction=0.):
		try:
			if self.mat.isAssembled():
				self.reinitialize()
		except: pass
		# create the Hamiltonian
		Hamiltonian = _PETSc.Log.Stage('Hamiltonian')
		Hamiltonian.push()
		
		ctqwmpi.importadjtoh(self.mat.fortran,filename,p,d,amp,interaction)
		
		Hamiltonian.pop()

		# create the adjacency matrix
		self.Adj = func.loadMat(filename,filetype,delimiter=delimiter)
		self.nodePos, self.lineX, self.lineY = func.getGraphNodes(self.Adj,layout=layout)
	
	def createLine2P(self,d=[0],amp=[0.],interaction=0.):
		try:
			if self.mat.isAssembled():
				self.reinitialize()
		except: pass
		self.defectNodes = d
		self.defectAmp = amp
		# create the Hamiltonian
		Hamiltonian = _PETSc.Log.Stage('Hamiltonian')
		Hamiltonian.push()
		ctqwmpi.hamiltonian_2p_line(self.mat.fortran,self.defectNodes,self.defectAmp,interaction,self.N)
		Hamiltonian.pop()

	def createLine3P(self,d=[0],amp=[0.],interaction=0.):
		try:
			if self.mat.isAssembled():
				self.reinitialize()
		except: pass
		self.defectNodes = d
		self.defectAmp = amp
		# create the Hamiltonian
		Hamiltonian = _PETSc.Log.Stage('Hamiltonian')
		Hamiltonian.push()
		ctqwmpi.hamiltonian_3p_line(self.mat.fortran,self.defectNodes,self.defectAmp,interaction,self.N)
		Hamiltonian.pop()
	
	def createLine(self,d=[0],amp=[0.]):
		try:
			if self.mat.isAssembled():
				self.reinitialize()
		except: pass
		self.defectNodes = d
		self.defectAmp = amp
		# create the Hamiltonian
		Hamiltonian = _PETSc.Log.Stage('Hamiltonian')
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
	def __init__(self,nodes,t,psi,psi2=None,psi3=None):
		self.rank = _PETSc.COMM_WORLD.Get_rank()
		self.time = [t]

		self.all_nodes = nodes

		self.local_nodes = []
		(Istart,Iend) = psi.getOwnershipRange()

		for i in self.all_nodes:
			if Istart <= i < Iend:
				self.local_nodes.append(i)

		self.local_prob = psi.getValues(self.local_nodes)

		if psi2 is not None:
			self.local_prob2 = psi2.getValues(self.local_nodes)

		if psi3 is not None:
			self.local_prob3 = psi3.getValues(self.local_nodes)

	def update(self,t,psi,psi2=None,psi3=None):
		self.time.append(t)

		prob_update = psi.getValues(self.local_nodes)
		self.local_prob = _np.vstack([self.local_prob,prob_update])

		if psi2 is not None:
			prob_update = psi2.getValues(self.local_nodes)
			self.local_prob2 = _np.vstack([self.local_prob2,prob_update])

		if psi3 is not None:
			prob_update = psi3.getValues(self.local_nodes)
			self.local_prob3 = _np.vstack([self.local_prob3,prob_update])

	def getLocalNodes(self,t=None,p=1):
		if t is not None:
			try:
				indt = self._time.index(t)
			except ValueError:
				if self.rank == 0:	print '\nERROR: time {} was not handled'.format(t)
				return

			if p == 1:
				return self.local_prob[indt]
			elif p == 2:
				return self.local_prob[indt], self.local_prob2[indt]
			elif p == 3:
				return self.local_prob[indt], self.local_prob2[indt], self.local_prob3[indt]

		else:
			if p == 1:
				return _np.array(self.time), self.local_prob.T.tolist()
			elif p == 2:
				return _np.array(self.time), self.local_prob.T.tolist(), self.local_prob2.T.tolist()
			elif p == 3:
				return _np.array(self.time), self.local_prob.T.tolist(), self.local_prob2.T.tolist(), self.local_prob3.T.tolist()

	def getLocalNode(self,node,p=1):
		try:
			indN = self.local_nodes.index(node)
		except ValueError:
			if p == 1:
				return _np.array(self.time), []
			elif p == 2:
				return _np.array(self.time), [], []
			elif p == 3:
				return _np.array(self.time), [], [], []

		if p == 1:
				return _np.array(self.time), self.local_prob.T.tolist()[indN]
		elif p == 2:
			return _np.array(self.time), self.local_prob.T.tolist()[indN], self.local_prob2.T.tolist()[indN]
		elif p == 3:
			return _np.array(self.time), self.local_prob.T.tolist()[indN], self.local_prob2.T.tolist()[indN], self.local_prob3.T.tolist()[indN]




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------- 1 particle CTQW   -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class QuantumWalkP1(object):
	def __init__(self,N):
		self.rank = _PETSc.Comm.Get_rank(_PETSc.COMM_WORLD)
		self.N = N
		self.t = 0
		self._timestep = False
		
		# define vectors
		self.psi0 = _PETSc.Vec()
		self.psi0.create(_PETSc.COMM_WORLD)
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
		initStateS = _PETSc.Log.Stage('initState')
		initStateS.push()
		try:
			self.psi0 = func.loadVec(filename,filetype)
			self.marginal(self.psi0)
		except:
			print '\nERROR: incorrect state (is it the correct length?'
			_sys.exit()
		initStateS.pop()

		self.marginal(self.psi0)

	def marginal(self,vec):
		# calculate marginal probabilities
		Marginal = _PETSc.Log.Stage('Marginal'); Marginal.push()
		ctqwmpi.p1prob(vec.fortran,self.prob.fortran,self.N)
		Marginal.pop()

	def watch(self,nodes,type='prob'):
		self.handle = nodeHandle(nodes,self.t,self.prob)
		
	def propagate(self,t,method='chebyshev',**kwargs):
		if self._timestep:
			self.t = t + self.t
		else:
			self.t = t

		self.EigSolver.setEigSolver(**kwargs)
		
		if method=='krylov':
			# SLEPc matrix exponential
			expmS = _PETSc.Log.Stage('SLEPc expm'); expmS.push()
			ctqwmpi.expm(self.H.mat.fortran,t,self.psi0.fortran,self.psi.fortran)
			expmS.pop()
			
		elif method=='chebyshev':
			# Chebyshev algorithm
			chebyS = _PETSc.Log.Stage('Chebyshev'); chebyS.push()
			ctqwmpi.qw_cheby(self.psi0.fortran,self.psi.fortran,t,self.H.mat.fortran,
					self.H.Emin(),self.H.Emax(),self.rank,self.N)
			chebyS.pop()
		
		self.marginal(self.psi)

		try:
			self.handle.update(self.t,self.prob)
		except:
			pass

		self._timestep = False

	def plotNodes(self,filename,t=None):

		comm = _PETSc.COMM_WORLD
		rank = comm.Get_rank()
		node = self.handle.local_nodes
		timeArray, probArray = self.handle.getLocalNodes(t=t)

		probArray = comm.tompi4py().gather(probArray)
		node = comm.tompi4py().gather(node)

		if rank == 0:
			timeArray = _np.array(timeArray)
			nodeArray = _np.array([item for sublist in node for item in sublist])
			probArray = _np.array([item for sublist in probArray for item in sublist]).real

			func.plotNodes(timeArray,nodeArray,probArray,filename)

	def exportState(self,filename,filetype):
		func.exportVec(self.psi,filename,filetype)

	def psiToInit(self):
		self.psi0 = self.psi
		self._timestep = True
	
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
		self.rank = _PETSc.Comm.Get_rank(_PETSc.COMM_WORLD)
		self.N = N
		self.t = 0
		self.timestep = False
		
		# define vectors
		self.psi0 = _PETSc.Vec()
		self.psi0.create(_PETSc.COMM_WORLD)
		self.psi0.setSizes(self.N**2)
		self.psi0.setUp()
		self.psi = self.psi0.duplicate()
		
		# create additional vectors
		self.psiX = _PETSc.Vec()
		self.psiX.create(_PETSc.COMM_WORLD)
		self.psiX.setSizes(self.N)
		self.psiX.setUp()
		self.psiY = self.psiX.duplicate()
		
		# define matrices
		self.H = Hamiltonian(self.N)
		self.EigSolver = self.H.EigSolver
	
	def importInitState(self,filename,filetype):
		self.initState = 'file:'+filename
		# create the inital stage
		initStateS = _PETSc.Log.Stage('initState')
		initStateS.push()
		try:
			if filetype == 'txt':
				self.psi0 = func.loadMatToVec(filename,filetype)
			elif filetype == 'bin':
				self.psi0 = func.loadVec(filename,filetype)
		except:
			print '\nERROR: incorrect state (is it the correct length?)'
			_sys.exit()
		initStateS.pop()

		self.marginal(self.psi0)

	def marginal(self,vec):
		# calculate marginal probabilities
		Marginal = _PETSc.Log.Stage('Marginal'); Marginal.push()
		ctqwmpi.marginal(vec.fortran,self.psiX.fortran,'x',self.N)
		ctqwmpi.marginal(vec.fortran,self.psiY.fortran,'y',self.N)
		Marginal.pop()

	def watch(self,nodes,type='prob'):
		self.handle = nodeHandle(nodes,self.t,self.psiX,psi2=self.psiY)
		
	def propagate(self,t,method='expm',**kwargs):
		if self.timestep:
			self.t = t + self.t
		else:
			self.t = t

		self.EigSolver.setEigSolver(**kwargs)
		
		if method=='krylov':
			# SLEPc matrix exponential
			krylov = _PETSc.Log.Stage('SLEPc krylov'); krylov.push()
			ctqwmpi.expm(self.H.mat.fortran,t,self.psi0.fortran,self.psi.fortran)
			krylov.pop()
			
		elif method=='chebyshev':
			# Chebyshev algorithm
			chebyS = _PETSc.Log.Stage('Chebyshev'); chebyS.push()
			ctqwmpi.qw_cheby(self.psi0.fortran,self.psi.fortran,t,self.H.mat.fortran,
					self.H.Emin(),self.H.Emax(),self.rank,self.N)
			chebyS.pop()
		
		self.marginal(self.psi)

		try:
			self.handle.update(self.t,self.psiX,psi2=self.psiY)
		except:
			pass

		self.timestep = False

	def plotNode(self,filename,node,t=None):

		comm = _PETSc.COMM_WORLD
		rank = comm.Get_rank()
		timeArray, probXArray, probYArray = self.handle.getLocalNode(node,p=2)

		probXArray = comm.tompi4py().gather(probXArray)
		probYArray = comm.tompi4py().gather(probYArray)

		if rank == 0:
			timeArray = _np.array(timeArray)
			probXArray = _np.array([item for sublist in probXArray for item in sublist]).real
			probYArray = _np.array([item for sublist in probYArray for item in sublist]).real

			func.plotNodes2P(timeArray,node,probXArray,probYArray,filename)

	def plotNodes(self,filename,p=1,t=None):

		comm = _PETSc.COMM_WORLD
		rank = comm.Get_rank()
		node = self.handle.local_nodes
		timeArray, probXArray, probYArray = self.handle.getLocalNodes(t=t,p=2)

		if p == 1:
			probArray = probXArray
		elif p==2:
			probArray = probYArray
		else:
			print 'p must be either 1 or 2'
			return

		probArray = comm.tompi4py().gather(probArray)
		node = comm.tompi4py().gather(node)

		if rank == 0:
			timeArray = _np.array(timeArray)
			nodeArray = _np.array([item for sublist in node for item in sublist])
			probArray = _np.array([item for sublist in probArray for item in sublist]).real

			func.plotNodes(timeArray,nodeArray,probArray,filename,p=p)
	
	def exportState(self,filename,filetype):
		if filetype == 'txt':
			func.exportVecToMat(self.psi,filename,filetype)
		elif filetype == 'bin':
			func.exportVec(self.psi,filename,filetype)

	def psiToInit(self):
		self.psi0 = self.psi
		self.timestep = True
	
	def destroy(self):
		self.H.destroy()
		self.psi.destroy()
		self.psi0.destroy()
		self.psiX.destroy()
		self.psiY.destroy()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------- 3 particle CTQW   -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class QuantumWalkP3(object):
	def __init__(self,N):
		self.rank = _PETSc.Comm.Get_rank(_PETSc.COMM_WORLD)
		self.N = N
		self.t = 0
		self.timestep = False
		
		# define vectors
		self.psi0 = _PETSc.Vec()
		self.psi0.create(_PETSc.COMM_WORLD)
		self.psi0.setSizes(self.N**3)
		self.psi0.setUp()
		self.psi = self.psi0.duplicate()
		
		# create additional vectors
		self.psiX = _PETSc.Vec()
		self.psiX.create(_PETSc.COMM_WORLD)
		self.psiX.setSizes(self.N)
		self.psiX.setUp()
		self.psiY = self.psiX.duplicate()
		self.psiZ = self.psiX.duplicate()
		
		# define matrices
		self.H = Hamiltonian(self.N)
		self.EigSolver = self.H.EigSolver
	
	def importInitState(self,filename,filetype):
		self.initState = 'file:'+filename
		# create the inital stage
		initStateS = _PETSc.Log.Stage('initState')
		initStateS.push()
		try:
			if filetype == 'txt':
				print 'Text file state import not supported for three walkers.'
				print 'Please try again with a binary file state import.'
				return
			elif filetype == 'bin':
				self.psi0 = func.loadVec(filename,filetype)
		except:
			print '\nERROR: incorrect state (is it the correct length?)'
			_sys.exit()
		initStateS.pop()

		self.marginal(self.psi0)

	def marginal(self,vec):
		# calculate marginal probabilities
		Marginal = _PETSc.Log.Stage('Marginal'); Marginal.push()
		ctqwmpi.marginal3(vec.fortran,self.psiX.fortran,'x',self.N)
		ctqwmpi.marginal3(vec.fortran,self.psiY.fortran,'y',self.N)
		ctqwmpi.marginal3(vec.fortran,self.psiZ.fortran,'z',self.N)
		Marginal.pop()

	def watch(self,nodes,type='prob'):
		self.handle = nodeHandle(nodes,self.t,self.psiX,psi2=self.psiY,psi3=self.psiZ)
		
	def propagate(self,t,method='expm',**kwargs):
		if self.timestep:
			self.t = t + self.t
		else:
			self.t = t

		self.EigSolver.setEigSolver(**kwargs)
		
		if method=='krylov':
			# SLEPc matrix exponential
			krylov = _PETSc.Log.Stage('SLEPc krylov'); krylov.push()
			ctqwmpi.expm(self.H.mat.fortran,t,self.psi0.fortran,self.psi.fortran)
			krylov.pop()
			
		elif method=='chebyshev':
			# Chebyshev algorithm
			chebyS = _PETSc.Log.Stage('Chebyshev'); chebyS.push()
			ctqwmpi.qw_cheby(self.psi0.fortran,self.psi.fortran,t,self.H.mat.fortran,
					self.H.Emin(),self.H.Emax(),self.rank,self.N)
			chebyS.pop()
		
		self.marginal(self.psi)

		try:
			self.handle.update(self.t,self.psiX,psi2=self.psiY,psi3=self.psiZ)
		except:
			pass

		self.timestep = False

	def plotNode(self,filename,node,t=None):

		comm = _PETSc.COMM_WORLD
		rank = comm.Get_rank()
		timeArray, probXArray, probYArray, probZArray = self.handle.getLocalNode(node,p=3)

		probXArray = comm.tompi4py().gather(probXArray)
		probYArray = comm.tompi4py().gather(probYArray)
		probZArray = comm.tompi4py().gather(probZArray)

		if rank == 0:
			timeArray = _np.array(timeArray)
			probXArray = _np.array([item for sublist in probXArray for item in sublist]).real
			probYArray = _np.array([item for sublist in probYArray for item in sublist]).real
			probZArray = _np.array([item for sublist in probZArray for item in sublist]).real

			func.plotNodes3P(timeArray,node,probXArray,probYArray,probZArray,filename)

	def plotNodes(self,filename,p=1,t=None):

		comm = _PETSc.COMM_WORLD
		rank = comm.Get_rank()
		node = self.handle.local_nodes
		timeArray, probXArray, probYArray, probZArray = self.handle.getLocalNodes(t=t,p=3)

		if p == 1:
			probArray = probXArray
		elif p==2:
			probArray = probYArray
		elif p==3:
			probArray = probZArray
		else:
			print 'p must be either 1, 2 or 3'
			return

		probArray = comm.tompi4py().gather(probArray)
		node = comm.tompi4py().gather(node)

		if rank == 0:
			timeArray = _np.array(timeArray)
			nodeArray = _np.array([item for sublist in node for item in sublist])
			probArray = _np.array([item for sublist in probArray for item in sublist]).real

			func.plotNodes(timeArray,nodeArray,probArray,filename,p=p)
	
	def exportState(self,filename,filetype):
		if filetype == 'txt':
			func.exportVec(self.psi,filename,filetype)
		elif filetype == 'bin':
			func.exportVec(self.psi,filename,filetype)

	def psiToInit(self):
		self.psi0 = self.psi
		self.timestep = True
	
	def destroy(self):
		self.H.destroy()
		self.psi.destroy()
		self.psi0.destroy()
		self.psiX.destroy()
		self.psiY.destroy()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------- Arbitrary CTQW --------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Graph(QuantumWalkP1):
	def __init__(self,N,filename=None,filetype=None,d=None,amp=None):
		QuantumWalkP1.__init__(self,N)
		self.liveplot = False
		
		if (filename and filetype) is not None:
			self.createH(filename,filetype,d=d,amp=amp)
		
	def createH(self,filename,filetype,d=None,amp=None,layout='spring',delimiter=None):
		if (d and amp) is None:
			self.defectNodes = [0]
			self.defectAmp = [0.]
		else:
			self.defectNodes = d
			self.defectAmp = amp
		self.H.importAdj(filename,filetype,d=self.defectNodes,amp=self.defectAmp,layout=layout,delimiter=delimiter)
		
	def createInitState(self,initState):
		self.initState = _np.vstack([_np.array(initState).T[0]-self.N/2+1,
			   	 	    _np.array(initState).T[1]]).T.tolist()
	
		# create the inital stage
		initStateS = _PETSc.Log.Stage('initState')
		initStateS.push()
		ctqwmpi.p1_init(self.psi0.fortran,self.initState,self.N)
		initStateS.pop()

		self.marginal(self.psi0)
		
	def plot(self,filename):
		if _os.path.isabs(filename):
			outDir = _os.path.dirname(filename)
		else:
			outDir = './'+_os.path.dirname(filename)
	
		# create output directory if it doesn't exist
		try:
			_os.mkdir(outDir)
		except OSError as exception:
			if exception.errno != _errno.EEXIST:
				raise

		initstateLabels = []
		for i in range(len(self.initState)):
			initstateLabels.append([sum(pair) for pair in zip(self.initState[i], [self.N/2-1,0])])

		plotStage = _PETSc.Log.Stage('Plotting'); plotStage.push()		
		func.plot(_np.arange(self.N),self.prob,filename,self.t,initstateLabels,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()

	def plotGraph(self,size=(12,8),probX=True,output=None,**kwargs):

		rank = _PETSc.COMM_WORLD.Get_rank()

		if probX:
			# scatter prob to process 0
			commX = self.prob.getComm()
			rank = _PETSc.COMM_WORLD.getRank()
			scatterX, probX0 = _PETSc.Scatter.toZero(self.prob)
			scatterX.scatter(self.prob, probX0, False, _PETSc.Scatter.Mode.FORWARD)

		if rank==0:
			from matplotlib import pyplot as plt
			import mpl_toolkits.mplot3d as plt3d
			self.fig = plt.figure(figsize=size)
			self.ax = plt3d.Axes3D(self.fig)
			self.ax.view_init(45, -50)
			self.ax.set_axis_off()

			if probX:
				prob = _np.real(_np.asarray(probX0))
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

		rank = _PETSc.COMM_WORLD.Get_rank()

		# scatter prob to process 0
		commX = self.prob.getComm()
		rank = _PETSc.COMM_WORLD.getRank()
		scatterX, probX0 = _PETSc.Scatter.toZero(self.prob)
		scatterX.scatter(self.prob, probX0, False, _PETSc.Scatter.Mode.FORWARD)

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

			prob = _np.real(_np.asarray(probX0))

			func.plotGraph(self.axLive,self.H.nodePos,self.H.lineX,self.H.lineY,
				prob=prob,**kwargs)
		
			self.axLive.set_title('$t={}$'.format(self.t))

			plt.draw()#show(block=True)
			_time.sleep(dt)

		# deallocate	
		commX.barrier()
		scatterX.destroy()
		probX0.destroy()

	def clearLiveGraph(self):
		rank = _PETSc.COMM_WORLD.Get_rank()

		if rank == 0:
			from matplotlib import pyplot as plt
			plt.show(block=True)
			self.liveplot = False
			plt.close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------- 2P Arbitrary CTQW -----------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Graph2P(QuantumWalkP2):
	def __init__(self,N,filename=None,filetype=None,d=None,amp=None,interaction=0.):
		QuantumWalkP2.__init__(self,N)
		self.liveplot = False
		
		if (filename and filetype) is not None:
			self.createH(filename,filetype,d=d,amp=amp,interaction=interaction)
		
	def createH(self,filename,filetype,d=None,amp=None,layout='spring',delimiter=None,interaction=0.):
		if (d or amp) is None:
			self.defectNodes = [0]
			self.defectAmp = [0.]
		else:
			self.defectNodes = d
			self.defectAmp = amp

		self.interaction = interaction

		self.H.importAdjToH(filename,filetype,
			d=self.defectNodes,amp=self.defectAmp,p='2',interaction=self.interaction,
			layout=layout,delimiter=delimiter)

		# ctqwmpi.importadjtoh(self.H.mat.fortran,filename,'2',
		# 	d=self.defectNodes,amp=self.defectAmp,layout=layout)
		
	def createInitState(self,initState):
		self.initState = _np.vstack([_np.array(initState).T[0]-self.N/2+1,
			   	 	    _np.array(initState).T[1]-self.N/2+1,_np.array(initState).T[2]]).T.tolist()
	
		# create the inital stage
		initStateS = _PETSc.Log.Stage('initState')
		initStateS.push()
		ctqwmpi.p2_init(self.psi0.fortran,self.initState,self.N)
		initStateS.pop()
		
	def plot(self,filename):
		if _os.path.isabs(filename):
			outDir = _os.path.dirname(filename)
		else:
			outDir = './'+_os.path.dirname(filename)
	
		# create output directory if it doesn't exist
		try:
			_os.mkdir(outDir)
		except OSError as exception:
			if exception.errno != _errno.EEXIST:
				raise

		initstateLabels = []
		for i in range(len(self.initState)):
			initstateLabels.append([sum(pair) for pair in zip(self.initState[i], [self.N/2-1,self.N/2-1,0])])

		plotStage = _PETSc.Log.Stage('Plotting'); plotStage.push()		
		func.plot2P(_np.arange(self.N),self.psiX,self.psiY,filename,self.t,initstateLabels,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()

	def plotGraph(self,size=(12,8),probX=True,probY=True,output=None,**kwargs):

		rank = _PETSc.COMM_WORLD.Get_rank()

		if probX:
			# scatter prob to process 0
			commX = self.psiX.getComm()
			rank = _PETSc.COMM_WORLD.getRank()
			scatterX, probX0 = _PETSc.Scatter.toZero(self.psiX)
			scatterX.scatter(self.psiX, probX0, False, _PETSc.Scatter.Mode.FORWARD)

		if probY:
			# scatter prob to process 0
			commY = self.psiY.getComm()
			rank = _PETSc.COMM_WORLD.getRank()
			scatterY, probY0 = _PETSc.Scatter.toZero(self.psiY)
			scatterX.scatter(self.psiY, probY0, False, _PETSc.Scatter.Mode.FORWARD)

		if rank==0:
			from matplotlib import pyplot as plt
			import mpl_toolkits.mplot3d as plt3d

			self.fig = plt.figure(figsize=size)
			self.ax = plt3d.Axes3D(self.fig)
			self.ax.view_init(45, -50)
			self.ax.set_axis_off()

			if probX:
				prob = _np.real(_np.asarray(probX0))
			else:
				prob = None

			if probY:
				prob2 = _np.real(_np.asarray(probY0))
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

		rank = _PETSc.COMM_WORLD.Get_rank()

		if probX:
			# scatter prob to process 0
			commX = self.psiX.getComm()
			rank = _PETSc.COMM_WORLD.getRank()
			scatterX, probX0 = _PETSc.Scatter.toZero(self.psiX)
			scatterX.scatter(self.psiX, probX0, False, _PETSc.Scatter.Mode.FORWARD)

		if probY:
			# scatter prob to process 0
			commY = self.psiY.getComm()
			rank = _PETSc.COMM_WORLD.getRank()
			scatterY, probY0 = _PETSc.Scatter.toZero(self.psiY)
			scatterX.scatter(self.psiY, probY0, False, _PETSc.Scatter.Mode.FORWARD)

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
				prob = _np.real(_np.asarray(probX0))
			else:
				prob = None

			if probY:
				prob2 = _np.real(_np.asarray(probY0))
			else:
				prob2 = None

			func.plotGraph(self.axLive,self.H.nodePos,self.H.lineX,self.H.lineY,
				prob=prob,prob2=prob2,**kwargs)
		
			self.axLive.set_title('$t={}$'.format(self.t))

			plt.draw()#show(block=True)
			_time.sleep(dt)

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
		rank = _PETSc.COMM_WORLD.Get_rank()

		if rank == 0:
			from matplotlib import pyplot as plt
			plt.show(block=True)
			self.liveplot = False
			plt.close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------- 3P Arbitrary CTQW -----------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Graph3P(QuantumWalkP3):
	def __init__(self,N,filename=None,filetype=None,d=None,amp=None,interaction=0.):
		QuantumWalkP3.__init__(self,N)
		
		if (filename and filetype) is not None:
			self.createH(filename,filetype,d=d,amp=amp,interaction=interaction)
		
	def createH(self,filename,filetype,d=None,amp=None,layout='spring',delimiter=None,interaction=0.):
		if (d and amp) is None:
			self.defectNodes = [0]
			self.defectAmp = [0.]
		else:
			self.defectNodes = d
			self.defectAmp = amp

		self.interaction = interaction

		self.H.importAdjToH(filename,filetype,
			d=self.defectNodes,amp=self.defectAmp,p='3',interaction=self.interaction,
			layout=layout,delimiter=delimiter)

		# ctqwmpi.importadjtoh(self.H.mat.fortran,filename,'2',
		# 	d=self.defectNodes,amp=self.defectAmp,layout=layout)
		
	def createInitState(self,initState):
		self.initState = _np.vstack(  [ _np.array(initState).T[0],
			   	 	    				_np.array(initState).T[1],
			   	 	    				_np.array(initState).T[2],
			   	 	    				_np.array(initState).T[3] ]).T.tolist()
	
		# create the inital stage
		initStateS = _PETSc.Log.Stage('initState')
		initStateS.push()
		ctqwmpi.p3_init(self.psi0.fortran,self.initState,self.N)
		initStateS.pop()
		
	def plot(self,filename):
		if _os.path.isabs(filename):
			outDir = _os.path.dirname(filename)
		else:
			outDir = './'+_os.path.dirname(filename)
	
		# create output directory if it doesn't exist
		try:
			_os.mkdir(outDir)
		except OSError as exception:
			if exception.errno != _errno.EEXIST:
				raise

		plotStage = _PETSc.Log.Stage('Plotting'); plotStage.push()		
		func.plot3P(_np.arange(self.N),self.psiX,self.psiY,self.psiZ,filename,self.t,self.initState,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()
		
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------- 1 particle CTQW on a line -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Line(QuantumWalkP1):
	def __init__(self,N,d=None,amp=None):
		QuantumWalkP1.__init__(self,N)
		if (d is not None) and (amp is not None):
			self.defectNodes = d
			self.defectAmp = amp
			self.H.createLine(d=d,amp=amp)
		
	def createH(self,d=None,amp=None):
		if (d and amp) is None:
			self.defectNodes = [0]
			self.defectAmp = [0.]
		else:
			self.defectNodes = d
			self.defectAmp = amp
		self.H.createLine(d=self.defectNodes,amp=self.defectAmp)
		
	def createInitState(self,initState):
		self.initState = initState
		# create the inital stage
		initStateS = _PETSc.Log.Stage('initState')
		initStateS.push()
		ctqwmpi.p1_init(self.psi0.fortran,initState,self.N)
		initStateS.pop()

	def watch(self,nodes,type='prob'):
		nodes = [i+self.N/2-1 for i in nodes]
		super(Line,self).watch(nodes,type=type)
		
	def plot(self,filename):
		if _os.path.isabs(filename):
			outDir = _os.path.dirname(filename)
		else:
			outDir = './'+_os.path.dirname(filename)
	
		# create output directory if it doesn't exist
		try:
			_os.mkdir(outDir)
		except OSError as exception:
			if exception.errno != _errno.EEXIST:
				raise

		plotStage = _PETSc.Log.Stage('Plotting'); plotStage.push()		
		func.plot(_np.arange(1-self.N/2,self.N/2+1),self.prob,filename,self.t,self.initState,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()

	def plotNodes(self,filename,t=None):

		comm = _PETSc.COMM_WORLD
		rank = comm.Get_rank()
		node = self.handle.local_nodes
		timeArray, probArray = self.handle.getLocalNodes(t=t)

		probArray = comm.tompi4py().gather(probArray)
		node = comm.tompi4py().gather(node)

		if rank == 0:
			timeArray = _np.array(timeArray)
			nodeArray = _np.array([item-self.N/2+1 for sublist in node for item in sublist])
			probArray = _np.array([item for sublist in probArray for item in sublist]).real

			func.plotNodes(timeArray,nodeArray,probArray,filename)
			
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------- 2 particle CTQW on a line -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Line2P(QuantumWalkP2):
	def __init__(self,N,d=None,amp=None,interaction=0.):
		QuantumWalkP2.__init__(self,N)

		if (d is not None) and (amp is not None):
			self.defectNodes = d
			self.defectAmp = amp
			self.H.createLine2P(d=d,amp=amp,interaction=interaction)
		
	def createH(self,d=None,amp=None,interaction=0.):
		if (d and amp) is None:
			self.defectNodes = [0]
			self.defectAmp = [0.]
		else:
			self.defectNodes = d
			self.defectAmp = amp
			
		self.interaction = interaction
		self.H.createLine2P(d=self.defectNodes,amp=self.defectAmp,interaction=self.interaction)
		
	def createInitState(self,initState):
		self.initState = initState
		# create the inital stage
		initStateS = _PETSc.Log.Stage('initState')
		initStateS.push()
		ctqwmpi.p2_init(self.psi0.fortran,initState,self.N)
		initStateS.pop()

	def watch(self,nodes,type='prob'):
		nodes = [i+self.N/2-1 for i in nodes]
		super(Line2P,self).watch(nodes,type=type)
		
	def plot(self,filename):
		if _os.path.isabs(filename):
			outDir = _os.path.dirname(filename)
		else:
			outDir = './'+_os.path.dirname(filename)
	
		# create output directory if it doesn't exist
		try:
			_os.mkdir(outDir)
		except OSError as exception:
			if exception.errno != _errno.EEXIST:
				raise

		plotStage = _PETSc.Log.Stage('Plotting'); plotStage.push()		
		func.plot2P(_np.arange(1-self.N/2,self.N/2+1),self.psiX,self.psiY,filename,self.t,self.initState,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()

	def plotNode(self,filename,node,t=None):
		node = node+self.N/2-1
		super(Line2P,self).plotNode(filename,node)

	def plotNodes(self,filename,p=1,t=None):

		comm = _PETSc.COMM_WORLD
		rank = comm.Get_rank()
		node = self.handle.local_nodes
		timeArray, probXArray, probYArray = self.handle.getLocalNodes(t=t,p=2)

		if p == 1:
			probArray = probXArray
		elif p==2:
			probArray = probYArray
		else:
			print 'p must be either 1 or 2'
			return

		probArray = comm.tompi4py().gather(probArray)
		node = comm.tompi4py().gather(node)

		if rank == 0:
			timeArray = _np.array(timeArray)
			nodeArray = _np.array([item-self.N/2+1 for sublist in node for item in sublist])
			probArray = _np.array([item for sublist in probArray for item in sublist]).real

			func.plotNodes(timeArray,nodeArray,probArray,filename,p=p)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------- 3 particle CTQW on a line -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Line3P(QuantumWalkP3):
	def __init__(self,N,d=None,amp=None,interaction=0.):
		QuantumWalkP3.__init__(self,N)

		if (d is not None) and (amp is not None):
			self.defectNodes = d
			self.defectAmp = amp
			self.H.createLine3P(d=d,amp=amp,interaction=interaction)
		
	def createH(self,d=None,amp=None,interaction=0.):
		if (d and amp) is None:
			self.defectNodes = [0]
			self.defectAmp = [0.]
		else:
			self.defectNodes = d
			self.defectAmp = amp

		self.interaction = interaction

		self.H.createLine3P(d=self.defectNodes,amp=self.defectAmp,interaction=self.interaction)
		
	def createInitState(self,initState):
		self.initState = _np.vstack(  [ _np.array(initState).T[0]+self.N/2-1,
			   	 	    				_np.array(initState).T[1]+self.N/2-1,
			   	 	    				_np.array(initState).T[2]+self.N/2-1,
			   	 	    				_np.array(initState).T[3] ]).T.tolist()

		# create the inital stage
		initStateS = _PETSc.Log.Stage('initState')
		initStateS.push()
		ctqwmpi.p3_init(self.psi0.fortran,self.initState,self.N)
		initStateS.pop()

	def watch(self,nodes,type='prob'):
		nodes = [i+self.N/2-1 for i in nodes]
		super(Line3P,self).watch(nodes,type=type)
		
	def plot(self,filename):
		if _os.path.isabs(filename):
			outDir = _os.path.dirname(filename)
		else:
			outDir = './'+_os.path.dirname(filename)
	
		# create output directory if it doesn't exist
		try:
			_os.mkdir(outDir)
		except OSError as exception:
			if exception.errno != _errno.EEXIST:
				raise

		initstateLabels = []
		for i in range(len(self.initState)):
			initstateLabels.append([sum(pair) for pair in zip(self.initState[i], [1-self.N/2,1-self.N/2,1-self.N/2,0])])

		plotStage = _PETSc.Log.Stage('Plotting'); plotStage.push()		
		func.plot3P(_np.arange(1-self.N/2,self.N/2+1),self.psiX,self.psiY,self.psiZ,filename,self.t,initstateLabels,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()

	def plotNode(self,filename,node,t=None):
		node = node+self.N/2-1
		super(Line3P,self).plotNode(filename,node)

	def plotNodes(self,filename,p=1,t=None):

		comm = _PETSc.COMM_WORLD
		rank = comm.Get_rank()
		node = self.handle.local_nodes
		timeArray, probXArray, probYArray, probZArray = self.handle.getLocalNodes(t=t,p=3)

		if p == 1:
			probArray = probXArray
		elif p==2:
			probArray = probYArray
		elif p==3:
			probArray = probZArray
		else:
			print 'p must be either 1, 2 or 3'
			return

		probArray = comm.tompi4py().gather(probArray)
		node = comm.tompi4py().gather(node)

		if rank == 0:
			timeArray = _np.array(timeArray)
			nodeArray = _np.array([item-self.N/2+1 for sublist in node for item in sublist])
			probArray = _np.array([item for sublist in probArray for item in sublist]).real

			func.plotNodes(timeArray,nodeArray,probArray,filename,p=p)