#!/usr/bin/python
"""
Additional ctqw related classes.
"""
import os as _os
import errno as _errno
import sys as _sys
import time as _time
import glob as _glob

from petsc4py import PETSc as _PETSc
from libpyctqw_MPI import ctqwmpi as _ctqwmpi
import io as _io
import plot as _plot
import numpy as _np

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------ Eigensolver object (PETSc matrix input) ------------------
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
			Emax,Emax_error,ierr = _ctqwmpi.min_max_eigs(self.__mat.fortran,self.__rank,'max',
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
			Emin,Emin_error,ierr = _ctqwmpi.min_max_eigs(self.__mat.fortran,self.__rank,'min',
				self.esolver,self.workType,self.workSize,self.tol,self.maxIt,self.verbose)
			EminStage.pop()
			if ierr==0:
				self.Emin_val = Emin
				self.Emin_err = Emin_error
			
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------ Hamiltonian object (grid size input) -----------------------
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
		
		self.Adj = _io.loadMat(filename,filetype,delimiter=delimiter)
		self.mat = _io.adjToH(self.Adj,d=d,amp=amp)
		
		Hamiltonian.pop()

		self.nodePos, self.lineX, self.lineY = _plot.getGraphNodes(self.Adj,layout=layout)

	def importAdjToH(self,filename,filetype,d=[0.],amp=[0.],p='1',layout='spring',delimiter=None,interaction=0.):
		try:
			if self.mat.isAssembled():
				self.reinitialize()
		except: pass
		# create the Hamiltonian
		Hamiltonian = _PETSc.Log.Stage('Hamiltonian')
		Hamiltonian.push()

		_ctqwmpi.importadjtoh(self.mat.fortran,filename,p,d,amp,interaction)
		
		Hamiltonian.pop()

		# create the adjacency matrix
		self.Adj = _io.loadMat(filename,filetype,delimiter=delimiter)
		self.nodePos, self.lineX, self.lineY = _plot.getGraphNodes(self.Adj,layout=layout)
	
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
		_ctqwmpi.hamiltonian_2p_line(self.mat.fortran,self.defectNodes,self.defectAmp,interaction,self.N)
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
		_ctqwmpi.hamiltonian_3p_line(self.mat.fortran,self.defectNodes,self.defectAmp,interaction,self.N)
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
		_ctqwmpi.hamiltonian_1p_line(self.mat.fortran,self.defectNodes,self.defectAmp,self.N)
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

class entanglementHandle(object):
	def __init__(self,t,psi,N,**kwargs):
		self.__default = {'esolver' : 'krylovschur',
			          'workType': 'null',
			          'workSize': '35',
			          'tol'     : 0.,
			          'maxIt'   : 0,
			          'verbose' : False,
			          'p'		: 2}

		self._rank = _PETSc.COMM_WORLD.Get_rank()
		self.time = [t]
		self.N = N

		for key,default in self.__default.iteritems():
			setattr(self, key, kwargs.get(key,default))

		if self.p==2:
			entInit, ierr = _ctqwmpi.entanglement(psi.fortran,self.N,
				self.esolver,self.workType,self.workSize,self.tol,self.maxIt,self.verbose)

		self.entanglement = [entInit]

	def update(self,t,psi):
		self.time.append(t)

		if self.p==2:
			entUpdate, ierr = _ctqwmpi.entanglement(psi.fortran,self.N,
					self.esolver,self.workType,self.workSize,self.tol,self.maxIt,self.verbose)

		self.entanglement.append(entUpdate)

	def getEntanglement(self,t=None):
		if t is not None:
			try:
				indt = self._time.index(t)
			except ValueError:
				if self._rank == 0:	print '\nERROR: time {} was not handled'.format(t)
				return

			return self.entanglement[indt]

		else:
			return _np.array(self.time), _np.array(self.entanglement)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------- 1 particle CTQW   ---------------------------------
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
		"""	Imports the initial state of the quantum walk from a file.

			:param str filename: path to the file containing the input state.

							.. note::
								The number of elements in the imported input state vector
								**must** match the number of nodes the Graph object 
								is initialized with.

			:param str filetype: the filetype of the imported adjacency matrix.

								* ``'txt'`` - an :math:`N` element column vector in text format.
								* ``'bin'`` - an :math:`N` element PETSc binary vector.

			:returns:	this creates a PETSc vector containing the initial state, accessed via

						>>> walk.psi0

			:rtype: ``petsc4py.PETSc.Vec()``
			"""		
		self.initState = 'file:'+filename
		# create the inital stage
		initStateS = _PETSc.Log.Stage('initState')
		initStateS.push()
		try:
			self.psi0 = _io.loadVec(filename,filetype)
			self.marginal(self.psi0)
		except:
			print '\nERROR: incorrect state (is it the correct length?'
			_sys.exit()
		initStateS.pop()

		self.marginal(self.psi0)

	def marginal(self,vec):
		# calculate marginal probabilities
		Marginal = _PETSc.Log.Stage('Marginal'); Marginal.push()
		_ctqwmpi.p1prob(vec.fortran,self.prob.fortran,self.N)
		Marginal.pop()

	def watch(self,nodes):
		"""	Creates a handle that watches node probability during propagation.

			:param array nodes:	the nodes to watch (e.g. ``[0,1,4]``).

			:returns:	creates a handle that can be
						accessed to retrieve node probabilities for various :math:`t`

						For example,

						>>> walk.watch([0,1,2,3,4])
						>>> walk.propagate(5.,method='chebyshev')
						>>> timeArray, probArray = walk.handle.getLocalNodes()

						.. warning::
							note that `walk.handle` attributes are **not** collective;
							if running on multiple nodes, only *local* values will be
							returned.

			:rtype: ``pyCTQW.MPI.ctqw.nodeHandle()``

			"""				
		self.handle = nodeHandle(nodes,self.t,self.prob)
		
	def propagate(self,t,method='chebyshev',**kwargs):
		"""	Propagates the quantum walk for time t.

			:param float t:	the timestep over which to propagate the state ``walk.psi0``

							.. math::
								\\left|\\psi\\right\\rangle = e^{-iHt}\\left|\\psi0\\right\\rangle

			:param str method:	the propagation algorithm to use
								(``'chebyshev'`` *(default)* or ``'krylov'``)

			:param kwargs:	EigSolver keywords can also be passed to the propagator;
							for more details of the available EigSolver properties,
							see :class:`pyCTQW.MPI.ctqw.EigSolver`.

							.. note::
								these EigSolver properties only take effect if
								`method='chebyshev'`.

			:returns:	this creates a PETSc vector containing the propagated state, accessed via

						>>> walk.psi

						as well as a PETSc vector containing the marginal probabilities,

						>>> walk.prob

			:rtype: ``petsc4py.PETSc.Vec()``
			"""				
		if self._timestep:
			self.t = t + self.t
		else:
			self.t = t

		self.EigSolver.setEigSolver(**kwargs)
		
		if method=='krylov':
			# SLEPc matrix exponential
			expmS = _PETSc.Log.Stage('SLEPc expm'); expmS.push()
			_ctqwmpi.expm(self.H.mat.fortran,t,self.psi0.fortran,self.psi.fortran)
			expmS.pop()
			
		elif method=='chebyshev':
			# Chebyshev algorithm
			chebyS = _PETSc.Log.Stage('Chebyshev'); chebyS.push()
			_ctqwmpi.qw_cheby(self.psi0.fortran,self.psi.fortran,t,self.H.mat.fortran,
					self.H.Emin(),self.H.Emax(),self.rank,self.N)
			chebyS.pop()
		
		self.marginal(self.psi)

		try:
			self.handle.update(self.t,self.prob)
		except:
			pass

		self._timestep = False

	def plotNodes(self,filename,t=None):
		"""	Creates a plot of the node probablities over time.

			:param str filename: the absolute/relative path to the desired output file.

			.. note::
				* :func:`watch` **must** be called prior to propagation in order for\
				the probabilities at the specified nodes to be stored.
				* if multiple nodes are watched, they will **all** be plotted.
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* Timesteps are recorded when propagations occur - to increase the number\
				of timesteps, increase the number of propagations.

			"""				

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

			_plot.plotNodes(timeArray,nodeArray,probArray,filename)

	def exportState(self,filename,filetype):
		"""	Exports the final state of the quantum walk to a file.

			:param str filename: path to desired output file.

			:param str filetype: the filetype of the exported adjacency matrix.

								* ``'txt'`` - an :math:`N` element column vector in text format.
								* ``'bin'`` - an :math:`N` element PETSc binary vector.
			"""		
		_io.exportVec(self.psi,filename,filetype)

	def psiToInit(self):
		"""	Copies the state :attr:`Graph.psi` to :attr:`Graph.psi0`.

			For example, this can be used to iterate the propagation
			over discrete timesteps: ::

				for i in range(10):
					walk.propagate(0.01,method='chebyshev')
					walk.psiToInit()

			"""				
		self.psi0 = self.psi
		self._timestep = True
	
	def destroy(self):
		"""	Destroys the 1 particle :class:`Graph` object, and all \
			associated PETSc matrices/vectors.

			>>> walk.destroy()

			"""	
		self.H.destroy()
		self.psi.destroy()
		self.psi0.destroy()
		self.prob.destroy()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------- 2 particle CTQW   ---------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class QuantumWalkP2(object):
	def __init__(self,N):
		self.rank = _PETSc.Comm.Get_rank(_PETSc.COMM_WORLD)
		self.N = N
		self.t = 0
		self.timestep = False
		self._watchProb = False
		self._watchEnt = False
		
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
				self.psi0 = _io.loadMatToVec(filename,filetype)
			elif filetype == 'bin':
				self.psi0 = _io.loadVec(filename,filetype)
		except:
			print '\nERROR: incorrect state (is it the correct length?)'
			_sys.exit()
		initStateS.pop()

		self.marginal(self.psi0)

	def marginal(self,vec):
		# calculate marginal probabilities
		Marginal = _PETSc.Log.Stage('Marginal'); Marginal.push()
		_ctqwmpi.marginal(vec.fortran,self.psiX.fortran,'x',self.N)
		_ctqwmpi.marginal(vec.fortran,self.psiY.fortran,'y',self.N)
		Marginal.pop()

	def watch(self,nodes,watchtype='prob',**kwargs):
		if watchtype == 'prob':
			self._watchProb = True
			self.handle = nodeHandle(nodes,self.t,self.psiX,psi2=self.psiY)
		elif watchtype == 'entanglement':
			self._watchEnt = True
			self.entanglementHandle = entanglementHandle(self.t,self.psi0,self.N,**kwargs)
		
	def propagate(self,t,method='expm',**kwargs):
		if self.timestep:
			self.t = t + self.t
		else:
			self.t = t

		self.EigSolver.setEigSolver(**kwargs)
		
		if method=='krylov':
			# SLEPc matrix exponential
			krylov = _PETSc.Log.Stage('SLEPc krylov'); krylov.push()
			_ctqwmpi.expm(self.H.mat.fortran,t,self.psi0.fortran,self.psi.fortran)
			krylov.pop()
			
		elif method=='chebyshev':
			# Chebyshev algorithm
			chebyS = _PETSc.Log.Stage('Chebyshev'); chebyS.push()
			_ctqwmpi.qw_cheby(self.psi0.fortran,self.psi.fortran,t,self.H.mat.fortran,
					self.H.Emin(),self.H.Emax(),self.rank,self.N)
			chebyS.pop()
		
		self.marginal(self.psi)

		if self._watchProb:
			self.handle.update(self.t,self.psiX,psi2=self.psiY)

		if self._watchEnt:
			self.entanglementHandle.update(self.t,self.psi)

		self.timestep = False

	def plotEntanglement(self,filename):

		comm = _PETSc.COMM_WORLD
		rank = comm.Get_rank()

		if rank == 0:
			timeArray, entArray = self.entanglementHandle.getEntanglement()
			_plot.plotEntanglement(timeArray,entArray,filename,self.initState,self.defectNodes,self.defectAmp)

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

			_plot.plotNodes2P(timeArray,node,probXArray,probYArray,filename)

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

			_plot.plotNodes(timeArray,nodeArray,probArray,filename,p=p)
	
	def exportState(self,filename,filetype):
		if filetype == 'txt':
			_io.exportVecToMat(self.psi,filename,filetype)
		elif filetype == 'bin':
			_io.exportVec(self.psi,filename,filetype)

	def exportPartialTrace(self,filename,filetype,p=1):
		if p == 1:
			if filetype == 'bin':
				rho = _PETSc.Mat().create(_PETSc.COMM_WORLD)
				_ctqwmpi.partial_trace_mat(self.psi.fortran,rho.fortran,self.N)
				_io.exportMat(rho, filename, filetype)
				rho.destroy()
			else:
				rho = _ctqwmpi.partial_trace_array(self.psi.fortran,self.N)
				if self.rank == 0:	_np.savetxt(filename, rho.real)
		if p == 2:
			#to do
			pass

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
#--------------------------- 3 particle CTQW   ---------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class QuantumWalkP3(object):
	def __init__(self,N):
		self.rank = _PETSc.Comm.Get_rank(_PETSc.COMM_WORLD)
		self.N = N
		self.t = 0
		self.timestep = False
		self._watchProb = False
		self._watchEnt = False
		
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
				self.psi0 = _io.loadVec(filename,filetype)
		except:
			print '\nERROR: incorrect state (is it the correct length?)'
			_sys.exit()
		initStateS.pop()

		self.marginal(self.psi0)

	def marginal(self,vec):
		# calculate marginal probabilities
		Marginal = _PETSc.Log.Stage('Marginal'); Marginal.push()
		_ctqwmpi.marginal3(vec.fortran,self.psiX.fortran,'x',self.N)
		_ctqwmpi.marginal3(vec.fortran,self.psiY.fortran,'y',self.N)
		_ctqwmpi.marginal3(vec.fortran,self.psiZ.fortran,'z',self.N)
		Marginal.pop()

	def watch(self,nodes,watchtype='prob',**kwargs):
		if watchtype == 'prob':
			self._watchProb = True
			self.handle = nodeHandle(nodes,self.t,self.psiX,psi2=self.psiY,psi3=self.psiZ)
		elif watchtype == 'entanglement':
			self._watchEnt = True
			self.entanglementHandle = entanglementHandle(self.t,self.psi0,self.N,**kwargs)
		
	def propagate(self,t,method='expm',**kwargs):
		if self.timestep:
			self.t = t + self.t
		else:
			self.t = t

		self.EigSolver.setEigSolver(**kwargs)
		
		if method=='krylov':
			# SLEPc matrix exponential
			krylov = _PETSc.Log.Stage('SLEPc krylov'); krylov.push()
			_ctqwmpi.expm(self.H.mat.fortran,t,self.psi0.fortran,self.psi.fortran)
			krylov.pop()
			
		elif method=='chebyshev':
			# Chebyshev algorithm
			chebyS = _PETSc.Log.Stage('Chebyshev'); chebyS.push()
			_ctqwmpi.qw_cheby(self.psi0.fortran,self.psi.fortran,t,self.H.mat.fortran,
					self.H.Emin(),self.H.Emax(),self.rank,self.N)
			chebyS.pop()
		
		self.marginal(self.psi)

		if self._watchProb:
			self.handle.update(self.t,self.psiX,psi2=self.psiY,psi3=self.psiZ)

		if self._watchEnt:
			self.entanglementHandle.update(self.t,self.psi)

		self.timestep = False

	def plotEntanglement(self,filename):

		comm = _PETSc.COMM_WORLD
		rank = comm.Get_rank()

		if rank == 0:
			timeArray, entArray = self.entanglementHandle.getEntanglement()
			_plot.plotEntanglement(timeArray,entArray,filename,self.initState,self.defectNodes,self.defectAmp)

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

			_plot.plotNodes3P(timeArray,node,probXArray,probYArray,probZArray,filename)

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

			_plot.plotNodes(timeArray,nodeArray,probArray,filename,p=p)
	
	def exportState(self,filename,filetype):
		if filetype == 'txt':
			_io.exportVec(self.psi,filename,filetype)
		elif filetype == 'bin':
			_io.exportVec(self.psi,filename,filetype)

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
	"""	Performs and analyses 1 particle continuous-time quantum walks on graphs

		:param int N:	number of nodes to initialize the walker with. Nodes are \
						labeled :math:`j\\in\\{0,N-1\\}`.

		For example, to create a CTQW Graph object for a 10 node graph,

		>>> walk = pyCTQW.MPI.Graph(10)

		.. note::	**if filename *and* filetype are provided**, this automatically
					creates a PETSc Hamiltonian matrix, neglecting the need to run
					:func:`createH`. For details on the other keyword arguments,
					see :func:`createH`.
		"""

	def __init__(self,N,filename=None,filetype=None,d=None,amp=None,layout='spring',delimiter=None):
		QuantumWalkP1.__init__(self,N)
		self.liveplot = False
		
		if (filename and filetype) is not None:
			self.createH(filename,filetype,d=d,amp=amp,layout=layout,delimiter=delimiter)
		
	def createH(self,filename,filetype,d=None,amp=None,layout='spring',delimiter=None):
		"""	Generate the Hamiltonian of the graph.

			.. note::	this needs to be called **only** if the filename and filetype
						of the graph were not already called when the Graph object
						was initialized.

			:param str filename: path to the file containing the adjacency matrix of the graph

							.. note::
								The number of nodes in the imported adjacency matrix
								**must** match the number of nodes the Graph object 
								is initialized with.

			:param str filetype: the filetype of the imported adjacency matrix.

								* ``'txt'`` - an :math:`N\\times N` dense 2D array in text format.
								* ``'bin'`` - an :math:`N\\times N` PETSc binary matrix.

			:param array d:	an array containing *integers* indicating the nodes
							where diagonal defects are to be placed.

							>>> d = [0,1,4]

			:param array amp:	an array containing *floats* indicating the diagonal defect
							amplitudes corresponding to each element in ``d``.

							>>> amp = [0.5,-1.,4.2]
							>>> len(d) == len(amp)
							True

							.. warning:: the size of ``a`` and ``d`` must be identical

			:param str layout:	the format to store the position of the nodes (only used
								when running :func:`plotGraph`).

								* ``spring`` *(default)* - spring layout.
								* ``circle`` - nodes are arranged in a circle.
								* ``spectral`` - nodes are laid out according to the \
												spectrum of the graph.
								* ``random`` - nodes are arranged in a random pattern.

			:param str delimiter: this is passed to `numpy.genfromtxt\
					<http://docs.scipy.org/doc/numpy/reference/generated/numpy.genfromtxt.html>`_
					in the case of strange delimiters in an imported ``txt`` file.


			:returns:	this creates a PETSc Hamiltonian matrix, accessed via

						>>> walk.H.mat

			:rtype: ``petsc4py.PETSc.Mat()``
			"""
		if (d and amp) is None:
			self.defectNodes = [0]
			self.defectAmp = [0.]
		else:
			self.defectNodes = d
			self.defectAmp = amp
		self.H.importAdj(filename,filetype,d=self.defectNodes,amp=self.defectAmp,layout=layout,delimiter=delimiter)
		
	def createInitState(self,initState):
		"""	Generate the initial state of the quantum walk.

			:param array initState:	an :math:`n\\times 2` array, containing the initial
									state of the quantum walker in the format::

										[[j1,amp1],[j2,amp2],...]

									For example, for a CTQW initially located in a
									superposition of nodes 1 and 2, e.g.

									.. math::
										\\left|\\psi(0)\\right\\rangle = \\frac{1}{\\sqrt{2}}\\left|1\\right\\rangle
										- \\frac{1}{\\sqrt{2}} \\left|2\\right\\rangle

									the initial state would be created like so:

									>>> import numpy as np
									>>> init_state = [[1,1./np.sqrt(2.)], [2,-1./np.sqrt(2.)]]
									>>> walk.createInitState(init_state)

			:returns:	this creates a PETSc vector containing the initial state, accessed via

						>>> walk.psi0

			:rtype: ``petsc4py.PETSc.Vec()``
			"""		
		self.initState = _np.vstack([_np.array(initState).T[0]-self.N/2+1,
			   	 	    _np.array(initState).T[1]]).T.tolist()
	
		# create the inital stage
		initStateS = _PETSc.Log.Stage('initState')
		initStateS.push()
		_ctqwmpi.p1_init(self.psi0.fortran,self.initState,self.N)
		initStateS.pop()

		self.marginal(self.psi0)
		
	def plot(self,filename):
		"""	Creates a plot of probability vs node.

			:param str filename: the absolute/relative path to the desired output file.

			.. note::
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* if :func:`propagate` has not been called, the probability of the \
				initial state is plotted (:math:`t=0`). Otherwise, the last propagated \
				time (``walk.t``) is used.
			"""		
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
		_plot.plot(_np.arange(self.N),self.prob,filename,self.t,initstateLabels,
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

			_plot.plotGraph(self.ax,self.H.nodePos,self.H.lineX,self.H.lineY,
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

			_plot.plotGraph(self.axLive,self.H.nodePos,self.H.lineX,self.H.lineY,
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

		# _ctqwmpi.importadjtoh(self.H.mat.fortran,filename,'2',
		# 	d=self.defectNodes,amp=self.defectAmp,layout=layout)
		
	def createInitState(self,initState):
		self.initState = _np.vstack([_np.array(initState).T[0]-self.N/2+1,
			   	 	    _np.array(initState).T[1]-self.N/2+1,_np.array(initState).T[2]]).T.tolist()
	
		# create the inital stage
		initStateS = _PETSc.Log.Stage('initState')
		initStateS.push()
		_ctqwmpi.p2_init(self.psi0.fortran,self.initState,self.N)
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
		_plot.plot2P(_np.arange(self.N),self.psiX,self.psiY,filename,self.t,initstateLabels,
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

			_plot.plotGraph(self.ax,self.H.nodePos,self.H.lineX,self.H.lineY,
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

			_plot.plotGraph(self.axLive,self.H.nodePos,self.H.lineX,self.H.lineY,
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

		# _ctqwmpi.importadjtoh(self.H.mat.fortran,filename,'2',
		# 	d=self.defectNodes,amp=self.defectAmp,layout=layout)
		
	def createInitState(self,initState):
		self.initState = _np.vstack(  [ _np.array(initState).T[0],
			   	 	    				_np.array(initState).T[1],
			   	 	    				_np.array(initState).T[2],
			   	 	    				_np.array(initState).T[3] ]).T.tolist()
	
		# create the inital stage
		initStateS = _PETSc.Log.Stage('initState')
		initStateS.push()
		_ctqwmpi.p3_init(self.psi0.fortran,self.initState,self.N)
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
		_plot.plot3P(_np.arange(self.N),self.psiX,self.psiY,self.psiZ,filename,self.t,self.initState,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()
		
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------- 1 particle CTQW on a line -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Line(QuantumWalkP1):
	"""	Performs and analyses 1 particle continuous-time quantum walks on an infinite line

		:param int N:	an **even** number of nodes to initialize the walker with. Nodes are \
						labeled :math:`j\\in\\{1-N/2,N/2\\}`.

		For example, to create a CTQW Line object for a 10 node line,

		>>> walk = pyCTQW.MPI.Line(10)

		"""
	def __init__(self,N):
		QuantumWalkP1.__init__(self,N)
		
	def createH(self,d=None,amp=None):
		"""	Generate the Hamiltonian of the graph.

			:param array d:	an array containing *integers* indicating the nodes
							where diagonal defects are to be placed.

							>>> d = [0,1,4]

			:param array amp:	an array containing *floats* indicating the diagonal defect
							amplitudes corresponding to each element in ``d``.

							>>> amp = [0.5,-1.,4.2]
							>>> len(d) == len(amp)
							True

							.. warning:: the size of ``a`` and ``d`` must be identical

			:returns:	this creates a PETSc Hamiltonian matrix, accessed via

						>>> walk.H.mat

			:rtype: ``petsc4py.PETSc.Mat()``
			"""
		if (d and amp) is None:
			self.defectNodes = [0]
			self.defectAmp = [0.]
		else:
			self.defectNodes = d
			self.defectAmp = amp
		self.H.createLine(d=self.defectNodes,amp=self.defectAmp)
		
	def createInitState(self,initState):
		"""	Generate the initial state of the quantum walk.

			:param array initState:	an :math:`n\\times 2` array, containing the initial
									state of the quantum walker in the format::

										[[j1,amp1],[j2,amp2],...]

									For example, for a CTQW initially located in a
									superposition of nodes -4 and 2, e.g.

									.. math::
										\\left|\\psi(0)\\right\\rangle = \\frac{1}{\\sqrt{2}}\\left|-4\\right\\rangle
										- \\frac{1}{\\sqrt{2}} \\left|2\\right\\rangle

									the initial state would be created like so:

									>>> import numpy as np
									>>> init_state = [[-4,1./np.sqrt(2.)], [1,-1./np.sqrt(2.)]]
									>>> walk.createInitState(init_state)

			:returns:	this creates a PETSc vector containing the initial state, accessed via

						>>> walk.psi0

			:rtype: ``petsc4py.PETSc.Vec()``
			"""		
		self.initState = initState
		# create the inital stage
		initStateS = _PETSc.Log.Stage('initState')
		initStateS.push()
		_ctqwmpi.p1_init(self.psi0.fortran,initState,self.N)
		initStateS.pop()

	def watch(self,nodes,type='prob'):
		"""	Creates a handle that watches node probability during propagation.

			:param array nodes:	the nodes to watch (e.g. ``[-5,1,4]``).

			:returns:	creates a handle that can be
						accessed to retrieve node probabilities for various :math:`t`

						For example,

						>>> walk.watch([0,1,2,3,4])
						>>> walk.propagate(5.,method='chebyshev')
						>>> timeArray, probArray = walk.handle.getLocalNodes()

						.. warning::
							note that `walk.handle` attributes are **not** collective;
							if running on multiple nodes, only *local* values will be
							returned.

			:rtype: ``pyCTQW.MPI.ctqw.nodeHandle()``

			"""				
		nodes = [i+self.N/2-1 for i in nodes]
		super(Line,self).watch(nodes,type=type)
		
	def plot(self,filename):
		"""	Creates a plot of probability vs node.

			:param str filename: the absolute/relative path to the desired output file.

			.. note::
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* if :func:`propagate` has not been called, the probability of the \
				initial state is plotted (:math:`t=0`). Otherwise, the last propagated \
				time (``walk.t``) is used.
			"""		
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
		_plot.plot(_np.arange(1-self.N/2,self.N/2+1),self.prob,filename,self.t,self.initState,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()

	def plotNodes(self,filename,t=None):
		"""	Creates a plot of the node probablities over time.

			:param str filename: the absolute/relative path to the desired output file.

			.. note::
				* :func:`watch` **must** be called prior to propagation in order for\
				the probabilities at the specified nodes to be stored.
				* if multiple nodes are watched, they will **all** be plotted.
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* Timesteps are recorded when propagations occur - to increase the number\
				of timesteps, increase the number of propagations.

			"""	
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

			_plot.plotNodes(timeArray,nodeArray,probArray,filename)
			
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
		_ctqwmpi.p2_init(self.psi0.fortran,initState,self.N)
		initStateS.pop()

	def watch(self,nodes,watchtype='prob',**kwargs):
		if watchtype == 'prob':
			nodes = [i+self.N/2-1 for i in nodes]
		super(Line2P,self).watch(nodes,watchtype=watchtype,**kwargs)
		
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
		_plot.plot2P(_np.arange(1-self.N/2,self.N/2+1),self.psiX,self.psiY,filename,self.t,self.initState,
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

			_plot.plotNodes(timeArray,nodeArray,probArray,filename,p=p)

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
		_ctqwmpi.p3_init(self.psi0.fortran,self.initState,self.N)
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
		_plot.plot3P(_np.arange(1-self.N/2,self.N/2+1),self.psiX,self.psiY,self.psiZ,filename,self.t,initstateLabels,
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

			_plot.plotNodes(timeArray,nodeArray,probArray,filename,p=p)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------------------- Graph Isomorphism -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class GraphISO(object):
	"""	A graph isomorphism solver, containing functions for 
		creating graph certificates and checking isomorphism
		of adjacency matrices.

		>>> gi = pyCTQW.MPI.GraphISO(p=2,propagator='krylov')

		:param int p:	number of particles, 2 *(default)* or 3,
					to use in constructing the graph certificate.

		:param float freqTol: the tolerance to use when constructing the
						frequency table (*default* ``1.e-2``).

						.. note::
							For ``freqTol=1.e-2``, all decimal places
							below 0.01 are discarded from the probability
							distribution.
						.. seealso:: :py:func:`GIcert`

		:param float compareTol:	the tolerance used when comparing two Graph
							certificates (*default* ``1.e-10``).

							.. note::
								Two isomorphic certificates satisfy
								:math:`\max(|cert_1 - cert_2|) < compareTol`
							.. seealso:: :py:func:`isomorphicQ`

		:param str propagator:	the CTQW propagator algorithm to use
							when calculating the graph certificate
							(``'chebyshev'`` *(default)* or ``'krylov'``).

		:returns: GraphISO solver
		"""

	def __init__(self,**kwargs):		
		self.__default = {
						'p'				: 2,
	                    'freqTol'		: 1.e-2,
	                    'compareTol'    : 1.e-10,
	                    'propagator'	: 'chebyshev'
						}

		self.__eigDefault = {
						'esolver'		: 'krylovschur',
						'emax_estimate' : 0.,
						'workType'		: 'null',
						'workSize'		: '35',
						'tol'     		: 0.,
						'maxIt'   		: 0,
						'verbose' 		: False
						}
			   
		self.__rank = _PETSc.Comm.Get_rank(_PETSc.COMM_WORLD)
			   
		for key, default in self.__default.iteritems():
			setattr(self, key, kwargs.get(key,default))

		for key, default in self.__eigDefault.iteritems():
			setattr(self, key, kwargs.get(key,default))

	def setProperties(self,**kwargs):
		"""
		Set some or all of the GraphISO properties.
		For the list of the properties see :class:`GraphISO`.
		"""
		for key in kwargs:
			if self.__default.has_key(key):
				setattr(self, key, kwargs.get(key))
			else:
				print 'Property {} does not exist!'.format(key)
			
	def getProperty(self,*args):
		"""
		Get some or all of the GraphISO properties.
		For the list of the properties see :class:`GraphISO`.
		"""
		for key in args:
			if self.__default.has_key(key):
				return getattr(self, key)
			else:
				print 'Property {} does not exist!'.format(key)

	def setEigSolver(self,**kwargs):
		"""Set some or all of the eigenvalue solver properties.

			:param str esolver: the eigensolver algorithm to use. 

							* ``'krylovschur'`` *(default)* - Krylov-Schur
							* ``'arnoldi'`` - Arnoldi Method
							* ``'lanczos'`` - Lanczos Method
							* ``'power'`` - Power Iteration/Rayleigh Quotient Iteration
							* ``'gd'`` - Generalized Davidson
							* ``'jd'`` - Jacobi-Davidson,
							* ``'lapack'`` - Uses LAPACK eigenvalue solver routines
							* ``'arpack'`` - *only available if SLEPc is\
												compiled with ARPACK linking*

			:param str workType:	can be used to set the eigensolver worktype
								(either ``'ncv'`` or ``'mpd'``). The default
								is to let SLEPc decide.

			:param int workSize:	sets the work size **if** ``workType`` is set.

			:param float tolIn:	tolerance of the eigenvalue solver
							(*default* ``0.`` (SLEPc decides)).
			
			:param int maxIt:	maximum number of iterations of the eigenvalue solver
							(*default* ``0`` (SLEPc decides)).
			
			:param bool verbose: if ``True``, writes eigensolver information to the console

			:param float emax_estimate:	used to override the calculation
									of the graphs maximum eigenvalue.

						.. warning::	the supplied :math:`\hat{\lambda}_{\max}`
										**must** satisfy :math:`\hat{\lambda}
										_{\max}\geq\lambda_{\max}`. The greater
										the value of :math:`\hat{\lambda}_{\max}
										-\lambda_{\max}`, the longer the convergence
										time of the ``chebyshev`` propagator

			.. note::
				* These properties only apply if ``propagator='chebyshev'``
				* For more information regarding these properties,refer to \
				the `SLEPc documentation\
				<http://www.grycap.upv.es/slepc/documentation/manual.htm>`_

			"""
		for key in kwargs:
			if self.__eigDefault.has_key(key):
				setattr(self, key, kwargs.get(key))
			else:
				print 'Eigsolver property {} does not exist!'.format(key)
			
	def getEigSolver(self,*args):
		"""
		Get some or all of the GraphISO properties.
		For the list of the properties see :func:`setEigSolver`.
		"""
		for key in args:
			if self.__eigDefault.has_key(key):
				return getattr(self, key)
			else:
				print 'Eigsolver property {} does not exist!'.format(key)

	def GIcert(self,adj):
		"""Generate the GI certificate of a graph.

			:param adj:	symmetric adjacency matrix in the form
						of a dense array of numpy array.
			:type adj:	array or numpy.array

			:returns: graph certificate
			:rtype: numpy.array
			"""

		cert, certSize = _ctqwmpi.GraphISCert(adj,self.p,self.freqTol,self.propagator,
			self.esolver,self.emax_estimate,self.workType,self.workSize,self.tol,self.maxIt,self.verbose)

		return _np.array(cert).T[_np.lexsort(_np.array(cert)[:,0:certSize])[::-1]]

	def isomorphicQ(self,adj1,adj2):
		"""Returns ``True`` if two graphs are isomorphic.

			:param adj1,adj2:	symmetric adjacency matrices in the form
								of a dense array of numpy array.
			:type adj1,adj2:	array or numpy.array

			:rtype: bool
			"""

		GIcert1 = self.GIcert(adj1)
		GIcert2 = self.GIcert(adj2)

		comm = _PETSc.COMM_WORLD
		size = comm.getSize()
		result = None

		if self.__rank == 0:
			if GIcert1.shape[0] == GIcert2.shape[0]:
				if _np.abs(_np.subtract(GIcert1,GIcert2)).max() < self.compareTol:
					result = _np.array([True for i in range(size)])
				else:
					result = _np.array([False for i in range(size)])
			else:
				result = _np.array([False for i in range(size)])

		result = comm.tompi4py().scatter(result)
		return result

	def AllIsomorphicQ(self,folder,graphRange=None,info=True,checkSelf=False):
		"""Calculates whether each pair of graphs (in a specified set of graphs)
			are isomorphic, returning an array :math:`R` with :math:`R_{ij} = 1` if
			graphs :math:`i` and :math:`j` are isomorphic, and :math:`R_{ij} = 0` otherwise.

			:param str folder:	path to a folder containing a collection of adjacency
								matrices in dense text file format.

								.. note::
									The text files must have filenames of the form \*X.txt
									where X represents a number (of any number of digits).
									These are used to order the graphs.

			:param array graphRange:	an array containing graph numbers to test.
										By default, all graphs in a folder are tested.

			:param bool info:	if ``True``, information on each :math:`R_{ij}` comparison
								is printed to the console.

			:param bool info:	if ``True``, each graph is also tested against itself.

			:rtype: array
			"""

		if not _os.path.isdir(folder):
			if self.__rank == 0:
				print 'Directory does not exist!'
				return

		numbers = _glob.re.compile(r'(\d+)')

		def numericalSort(value):
		    parts = numbers.split(value)
		    parts[1::2] = map(int, parts[1::2])
		    return parts

		filelist = sorted(_glob.glob(folder+'/*[0-9].txt'), key=numericalSort)

		adj = []
		for graph in filelist:
			if graphRange is not None:
				if int(_glob.re.findall(r'\d+',graph)[-1]) in graphRange:
					adj.append(_np.genfromtxt(graph))
					if (self.__rank==0 and info):	print 'Adding graph ' + graph
			else:
				adj.append(_np.genfromtxt(graph))
				if (self.__rank==0 and info):	print 'Adding graph ' + graph

		NG = len(adj)
		comparisonTable = _np.zeros([NG,NG])

		for i in range(NG):
			for j in range(i if checkSelf else i+1,NG):
				if (self.__rank==0 and info):	print 'Testing graphs ' + str(i) + ',' + str(j)
				comparisonTable[i,j] = 1 if self.isomorphicQ(adj[i],adj[j]) else 0
				if (self.__rank==0 and info):	print '\tIsomorphic' if comparisonTable[i,j] == 1 else '\tNon-isomorphic'

		if checkSelf:
			return comparisonTable + comparisonTable.T - _np.diag(comparisonTable.diagonal())
		else:
			return comparisonTable + comparisonTable.T + _np.diag(_np.ones(NG))