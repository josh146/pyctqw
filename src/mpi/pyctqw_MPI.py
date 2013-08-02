#!/usr/bin/python
import sys, os, errno

from petsc4py import PETSc
from libpyctqw_MPI import ctqwmpi
import pyctqw_plots
import numpy as np

class Hamiltonian(object):
	def __init__(self,N):
		self.rank = PETSc.Comm.Get_rank(PETSc.COMM_WORLD)
		self.N = N
		self.Emax_val = None
		self.Emin_val = None

		# define matrices
		self.mat = PETSc.Mat()
		self.mat.create(PETSc.COMM_WORLD)
	
	def reinitialize(self):
		self.destroy()
		self.mat = PETSc.Mat()
		self.mat.create(PETSc.COMM_WORLD)
	
	def createLine2P(self,d=[0],amp=[0.]):
		self.reinitialize()
		self.defectNodes = d
		self.defectAmp = amp
		# create the Hamiltonian
		Hamiltonian = PETSc.Log.Stage('Hamiltonian')
		Hamiltonian.push()
		ctqwmpi.hamiltonian_2p_line(self.mat.fortran,self.defectNodes,self.defectAmp,self.N)
		Hamiltonian.pop()
	
	def Emax(self,**kwargs):
		if self.Emax_val is None:
			self.findEmax(**kwargs)
		return self.Emax_val
	
	def Emin(self,**kwargs):
		if self.Emin_val is None:
			self.findEmin(**kwargs)
		return self.Emin_val
	
	def findEmax(self,esolver='krylovschur',workType='null',
	    workSize=35,tol=0.,maxIt=0,verbose=False):
		# Calcaulate the eigenvalues
		EmaxStage = PETSc.Log.Stage('Emax'); EmaxStage.push()
		Emax,Emax_error,ierr = ctqwmpi.min_max_eigs(self.mat.fortran,self.rank,
					'max',esolver,workType,workSize,tol,maxIt,verbose)
		EmaxStage.pop()
		if ierr==0:
			self.Emax_val = Emax
			self.Emax_err = Emax_error
			return Emax, Emax_error
	
	def findEmin(self,esolver='krylovschur',workType='null',
	    workSize=35,tol=0.,maxIt=0,verbose=False):
		# Calcaulate the eigenvalues
		EminStage = PETSc.Log.Stage('Emin'); EminStage.push()
		Emin,Emin_error,ierr = ctqwmpi.min_max_eigs(self.mat.fortran,self.rank,
					'min',esolver,workType,workSize,tol,maxIt,verbose)
		EminStage.pop()
		if ierr==0:
			self.Emin_val = Emin
			self.Emin_err = Emin_error
			return Emin, Emin_error
		
	def destroy(self):
		self.mat.destroy()
		
	

class Line2P(object):
	def __init__(self,N,d=None,amp=None):
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
		if d or amp is None:
			self.H.createLine2P()
		else:
			self.defectNodes = d
			self.defectAmp = amp
			self.H.createLine2P(d,amp)
		
	def createH(self,d=[0],amp=[0.]):
		self.defectNodes = d
		self.defectAmp = amp
		self.H.createLine2P(d,amp)
		
	def createInitState(self,initState):
		self.initState = initState
		# create the inital stage
		initStateS = PETSc.Log.Stage('initState')
		initStateS.push()
		ctqwmpi.p2_init(self.psi0.fortran,initState,self.N)
		initStateS.pop()

	def marginal(self):
		# calculate marginal probabilities
		Marginal = PETSc.Log.Stage('Marginal'); Marginal.push()
		ctqwmpi.marginal(self.psi.fortran,self.psiX.fortran,'x',self.N)
		ctqwmpi.marginal(self.psi.fortran,self.psiY.fortran,'y',self.N)
		Marginal.pop()
		
	def plot(self,filename):
		(outDir,plotName) = os.path.split(filename)
	
		# create output directory if it doesn't exist
		try:
			os.mkdir(outDir)
		except OSError as exception:
			if exception.errno != errno.EEXIST:
				raise

		plotStage = PETSc.Log.Stage('Plotting'); plotStage.push()		
		pyctqw_plots.plot_marginal(self.psiX,self.psiY,plotName,self.t,self.initState,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()
		
	def propagate(self,t,method='expm',**kwargs):
		self.t = t
		
		if method=='expm':
			# SLEPc matrix exponential
			expmS = PETSc.Log.Stage('SLEPc expm'); expmS.push()
			ctqwmpi.expm(self.H.mat.fortran,self.t,self.psi0.fortran,self.psi.fortran)
			expmS.pop()
			
		elif method=='chebyshev':
			# Chebyshev algorithm
			chebyS = PETSc.Log.Stage('Chebyshev'); chebyS.push()
			ctqwmpi.qw_cheby(self.psi0.fortran,self.psi.fortran,self.t,self.H.mat.fortran,
					self.H.Emin(**kwargs),self.H.Emax(**kwargs),self.rank,self.N)
			chebyS.pop()
		
		self.marginal()
	
	def destroy(self):
		self.H.destroy()
		self.psi.destroy()
		self.psi0.destroy()
		self.psiX.destroy()
		self.psiY.destroy()
