#!/usr/bin/python
import sys, os, errno

from petsc4py import PETSc
from libpyctqw_MPI import ctqwmpi
import pyctqw_plots
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
		
	def importAdj(self,filename,filetype,d=[0],amp=[0.]):
		try:
			if self.mat.isAssembled():
				self.reinitialize()
		except: pass
		# create the Hamiltonian
		Hamiltonian = PETSc.Log.Stage('Hamiltonian')
		Hamiltonian.push()
		
		self.Adj = pyctqw_plots.loadMat(filename,filetype)
		self.mat = pyctqw_plots.adjToH(self.Adj,d=d,amp=amp)
		
		Hamiltonian.pop()

	def importAdjToH(self,filename,filetype,d=[0],amp=[0.],p='1'):
		try:
			if self.mat.isAssembled():
				self.reinitialize()
		except: pass
		# create the Hamiltonian
		Hamiltonian = PETSc.Log.Stage('Hamiltonian')
		Hamiltonian.push()
		
		ctqwmpi.importAdjToH(self.mat.fortran,filename,p,d,amp)
		
		Hamiltonian.pop()
	
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------- 1 particle CTQW   -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class QuantumWalkP1(object):
	def __init__(self,N):
		self.rank = PETSc.Comm.Get_rank(PETSc.COMM_WORLD)
		self.N = N
		
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
			self.psi0 = pyctqw_plots.loadVec(filename,filetype)
		except:
			print '\nERROR: incorrect state (is it the correct length?'
			sys.exit()
		initStateS.pop()

	def marginal(self):
		# calculate marginal probabilities
		Marginal = PETSc.Log.Stage('Marginal'); Marginal.push()
		ctqwmpi.p1prob(self.psi.fortran,self.prob.fortran,self.N)
		Marginal.pop()
		
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
		
		self.marginal()
	
	def exportState(self,filename,filetype):
		pyctqw_plots.exportVec(self.psi,filename,filetype)

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
				self.psi0 = pyctqw_plots.loadMatToVec(filename,filetype)
			elif filetype == 'bin':
				self.psi0 = pyctqw_plots.loadVec(filename,filetype)
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
			pyctqw_plots.exportVecToMat(self.psi,filename,filetype)
		elif filetype == 'bin':
			pyctqw_plots.exportVec(self.psi,filename,filetype)

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
		
		if (filename and filetype) is not None:
			self.createH(filename,filetype,d=d,amp=amp)
		
	def createH(self,filename,filetype,d=None,amp=None):
		if (d and amp) is None:
			self.defectNodes = [0]
			self.defectAmp = [0.]
		else:
			self.defectNodes = d
			self.defectAmp = amp
		self.H.importAdj(filename,filetype,d=self.defectNodes,amp=self.defectAmp)
		
	def createInitState(self,initState):
		self.initState = np.vstack([np.array(initState).T[0]-self.N/2+1,
			   	 	    np.array(initState).T[1]]).T.tolist()
	
		# create the inital stage
		initStateS = PETSc.Log.Stage('initState')
		initStateS.push()
		ctqwmpi.p1_init(self.psi0.fortran,self.initState,self.N)
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
			initstateLabels.append([sum(pair).real for pair in zip(self.initState[i], [self.N/2-1,0])])

		plotStage = PETSc.Log.Stage('Plotting'); plotStage.push()		
		pyctqw_plots.plot(np.arange(self.N),self.prob,filename,self.t,initstateLabels,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------- 2P Arbitrary CTQW -----------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class ctqwGraph2P(QuantumWalkP2):
	def __init__(self,N,filename=None,filetype=None,d=None,amp=None):
		QuantumWalkP2.__init__(self,N)
		
		if (filename and filetype) is not None:
			self.createH(filename,filetype,d=d,amp=amp)
		
	def createH(self,filename,filetype,d=None,amp=None):
		if (d and amp) is None:
			self.defectNodes = [0]
			self.defectAmp = [0.]
		else:
			self.defectNodes = d
			self.defectAmp = amp
		ctqwmpi.importadjtoh(self.H.mat.fortran,filename,'2',d=self.defectNodes,amp=self.defectAmp)
		
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
		pyctqw_plots.plot2P(np.arange(self.N),self.psiX,self.psiY,filename,self.t,initstateLabels,
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
		pyctqw_plots.plot(np.arange(1-self.N/2,self.N/2+1),self.prob,filename,self.t,self.initState,
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
		pyctqw_plots.plot2P(np.arange(1-self.N/2,self.N/2+1),self.psiX,self.psiY,filename,self.t,self.initState,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()

