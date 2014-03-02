#!/usr/bin/python
#
#  This file provides an interface from the pyCTQW.MPI python module
#  to the libctqwMPI.so FORTRAN library contained in source code src/ctqwMPI.F90
#  ----------------------------------------------------------------------------
#  pyCTQW - Distributed memory CTQW Fortran library and Python module
#  Copyright (C) 2013-2014, Joshua Izaac
#
#  pyCTQW is free software: you can redistribute it and/or modify it under the
#  terms of version 3 of the GNU General Public License as published by
#  the Free Software Foundation.
#
#  pyCTQW is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
#  more details.
#
#  You should have received a copy of the GNU General Public License
#  along with pyCTQW. If not, see <http://www.gnu.org/licenses/>.
#  ----------------------------------------------------------------------------
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
	""" Contains methods for setting up and solving for the boundary eigenvalues of a PETSc matrix.

		Args:
			mat (petsc4py.PETSc.Mat) : input matrix

		Keyword Args:
				esolver (str): the eigensolver algorithm to use. 

							* ``'krylovschur'`` *(default)* - Krylov-Schur
							* ``'arnoldi'`` - Arnoldi Method
							* ``'lanczos'`` - Lanczos Method
							* ``'power'`` - Power Iteration/Rayleigh Quotient Iteration
							* ``'gd'`` - Generalized Davidson
							* ``'jd'`` - Jacobi-Davidson,
							* ``'lapack'`` - Uses LAPACK eigenvalue solver routines
							* ``'arpack'`` - *only available if SLEPc is\
												compiled with ARPACK linking*

				workType (str):    can be used to set the eigensolver worktype \
								(either ``'ncv'`` or ``'mpd'``). The default \
								is to let SLEPc decide.

				workSize (int):    sets the work size **if** :attr:`workType` is set.

				tolIn (float): tolerance of the eigenvalue solver \
								(*default* ``0.`` (SLEPc decides)).
				
				maxIt (int):   maximum number of iterations of the eigenvalue solver \
								(*default* ``0`` (SLEPc decides)).
				
				verbose (bool): if ``True``, writes eigensolver information to the console

				emax_estimate (float): used to override the calculation \
										of the graphs maximum eigenvalue.

				emin_estimate (float): used to override the calculation \
										of the graphs minimum eigenvalue. For a *finite* graph, \
										``emin_estimate=0`` is set by default.

					.. caution::
						* If supplied, the value of :attr:`emax_estimate` (:math:`\hat{\lambda}_{\max}`) \
						**must** satisfy :math:`\hat{\lambda}_{\max}\geq\lambda_{\max}`, \
						where :math:`\lambda_{\max}` is the actual maximum eigenvalue of the graph.
						* Similarly, the value of :attr:`emin_estimate` (:math:`\hat{\lambda}_{\min}`) \
						**must** satisfy :math:`\hat{\lambda}_{\min}\leq\lambda_{\min}`, \
						where :math:`\lambda_{\min}` is the actual minimum eigenvalue of the graph.
						* The greater the difference value of :math:`|\hat{\lambda}_{\max} \
						-\lambda_{\max}|` and :math:`|\hat{\lambda}_{\min} \
						-\lambda_{\min}|`, the longer the convergence \
						time of the ``chebyshev`` propagator when simulating a CTQW.

					Note:
						* For more information regarding these properties,refer to \
						the `SLEPc documentation\
						<http://www.grycap.upv.es/slepc/documentation/manual.htm>`_		
		"""

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
		"""Set some or all of the eigenvalue solver properties.

			Keyword Args:
				esolver (str): the eigensolver algorithm to use. 

							* ``'krylovschur'`` *(default)* - Krylov-Schur
							* ``'arnoldi'`` - Arnoldi Method
							* ``'lanczos'`` - Lanczos Method
							* ``'power'`` - Power Iteration/Rayleigh Quotient Iteration
							* ``'gd'`` - Generalized Davidson
							* ``'jd'`` - Jacobi-Davidson,
							* ``'lapack'`` - Uses LAPACK eigenvalue solver routines
							* ``'arpack'`` - *only available if SLEPc is\
												compiled with ARPACK linking*

				workType (str):    can be used to set the eigensolver worktype
								(either ``'ncv'`` or ``'mpd'``). The default
								is to let SLEPc decide.

				workSize (int):    sets the work size **if** ``workType`` is set.

				tolIn (float): tolerance of the eigenvalue solver
								(*default* ``0.`` (SLEPc decides)).
				
				maxIt (int):   maximum number of iterations of the eigenvalue solver
								(*default* ``0`` (SLEPc decides)).
				
				verbose (bool): if ``True``, writes eigensolver information to the console

				emax_estimate (float): used to override the calculation
										of the graphs maximum eigenvalue.

			.. caution::
				* If supplied, the value of :attr:`emax_estimate`:math:`\hat{\lambda}_{\max}` \
				**must** satisfy :math:`\hat{\lambda}_{\max}\geq\lambda_{\max}`, \
				where :math:`\lambda_{\max}` is the actual maximum eigenvalue of the graph.
				* The greater the value of :math:`\hat{\lambda}_{\max} \
				-\lambda_{\max}`, the longer the convergence \
				time of the ``chebyshev`` propagator when simulating a CTQW.

			Note:
				* For more information regarding these properties,refer to \
				the `SLEPc documentation\
				<http://www.grycap.upv.es/slepc/documentation/manual.htm>`_

			"""
		for key in kwargs:
			if hasattr(self, key):
				setattr(self, key, kwargs.get(key))
			else:
				'Property type does not exist!'
			
	def getEigSolver(self,*args):
		"""	Get some or all of the GraphISO properties.
			For the list of the properties see :func:`setEigSolver`.
			"""
		for key in args:
			if hasattr(self, key):
				return getattr(self, key)
			else:
				'Property type does not exist!'

	def findEmax(self):
		""" Returns the maximum Eigenvalue of the matrix, along with associated error.

			Example:

				>>> Emax, error = EigSolver.findEmax()

			.. admonition:: Fortran interface
			
				This function calls the Fortran function :f:func:`min_max_eigs`.
			"""
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
		""" Returns the minimum Eigenvalue of the matrix, along with associated error.

			Example:

				>>> Emin, error = EigSolver.findEmin()

			.. admonition:: Fortran interface
			
				This function calls the Fortran function :f:func:`min_max_eigs`.
			"""
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
	""" Contains methods for initializing, creating and manipulating Hamiltonian matrices.

		Args:
			N (int) : number of nodes to initialize the Hamiltonian object with.

		Returns:
				:	creates an unallocated PETSc matrix on communicator \
					:attr:`petsc4py.PETSc.COMM_WORLD`. This is accessed \
					via the attribute :attr:`Hamiltonian.mat`.
		:rtype: petsc4py.PETSc.Mat

		Returns:
				:	an associated :class:`EigSolver` object is also created \
					at :attr:`Hamiltonian.EigSolver`.
		:rtype: :class:`EigSolver`
		"""

	def __init__(self,N):
		self.rank = _PETSc.Comm.Get_rank(_PETSc.COMM_WORLD)
		self.N = N

		# define matrices
		self.mat = _PETSc.Mat()
		self.mat.create(_PETSc.COMM_WORLD)
		
		# create eigenvalue solver
		self.EigSolver = EigSolver(self.mat)
	
	def reinitialize(self):
		""" Destroy and reinitialize the PETSc matrix :attr:`Hamiltonian.mat`."""
		self.destroy()
		self.mat = _PETSc.Mat()
		self.mat.create(_PETSc.COMM_WORLD)
		
	def importAdj(self,filename,filetype,d=[0],amp=[0.],layout='spring',delimiter=None):
		""" Create a Hamiltonian from an imported adjacency matrix.

			Args:
				filename (str): path to the file containing the adjacency matrix of the graph

				filetype (str): the filetype of the imported adjacency matrix.

								* ``'txt'`` - an :math:`N\\times N` dense 2D array in text format.
								* ``'bin'`` - an :math:`N\\times N` PETSc binary matrix.

				d (array of ints): an array containing *integers* indicating the nodes
							where diagonal defects are to be placed (e.g. ``d=[0,1,4]``).

				amp (array of floats):   an array containing *floats* indicating the diagonal defect
							amplitudes corresponding to each element in ``d`` (e.g. ``amp=[0.5,-1,4.2]``).

				layout (str):  the format to store the position of the nodes.

								* ``spring`` *(default)* - spring layout.
								* ``circle`` - nodes are arranged in a circle.
								* ``spectral`` - nodes are laid out according to the \
												spectrum of the graph.
								* ``random`` - nodes are arranged in a random pattern.

				delimiter (str): this is passed to `numpy.genfromtxt\
					<http://docs.scipy.org/doc/numpy/reference/generated/numpy.genfromtxt.html>`_
					in the case of strange delimiters in an imported ``txt`` file.

			Returns:
				:this allocates and builds a Hamiltonian matrix at :attr:`Hamiltonian.mat`.
			:rtype: :func:`petsc4py.PETSc.Mat`

			Returns:
				: the original adjacency matrix is also stored as a PETSc matrix at :attr:`Hamiltonian.Adj`.
			:rtype: :func:`petsc4py.PETSc.Mat`
			
			.. important::
				* The number of nodes in the imported adjacency matrix \
				**must** match the number of nodes the Hamiltonian object \
				is initialized with.
				* The length of ``amp`` and ``d`` must be identical.

			.. note::
				* This function calls Python functions located at :mod:`pyCTQW.MPI.io`, \
				and is somewhat slower than :func:`importAdjToH`, but with the advantage of \
				a more flexible text file importer.
				* For the ability to generate 2 or 3 particle Hamiltonians with or without \
				interactions, please see :func:`importAdjToH`.
			"""
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
		""" Create a Hamiltonian from an imported adjacency matrix.

			Args:
				filename (str): path to the file containing the adjacency matrix of the graph

				filetype (str): the filetype of the imported adjacency matrix.

								* ``'txt'`` - an :math:`N\\times N` dense 2D array in text format.
								* ``'bin'`` - an :math:`N\\times N` PETSc binary matrix.

				d (array of ints): an array containing *integers* indicating the nodes
							where diagonal defects are to be placed (e.g. ``d=[0,1,4]``).

				amp (array of floats):   an array containing *floats* indicating the diagonal defect
							amplitudes corresponding to each element in ``d`` (e.g. ``amp=[0.5,-1,4.2]``).

				p (str) : (``'1'``|``'2'``|``'3'``) - specify whether to create a 1, 2 or 3 \
							particle Hamiltonian.

				layout (str):  the format to store the position of the nodes.

								* ``spring`` *(default)* - spring layout.
								* ``circle`` - nodes are arranged in a circle.
								* ``spectral`` - nodes are laid out according to the \
												spectrum of the graph.
								* ``random`` - nodes are arranged in a random pattern.

				delimiter (str): this is passed to `numpy.genfromtxt \
					<http://docs.scipy.org/doc/numpy/reference/generated/numpy.genfromtxt.html>`_
					in the case of strange delimiters in an imported ``txt`` file.

				interaction (float): the amplitude of interaction between the walkers \
									 when located on the same vertex.

			Returns:
				 :	* :attr:`Hamiltonian.mat` *(petsc4py.PETSc.Mat)* - the Hamiltonian matrix.
				 	* :attr:`Hamiltonian.Adj` *(petsc4py.PETSc.Mat)* - the adjacency matrix stored as a PETSc matrix.
					* :attr:`Hamiltonian.nodePos` *(array)* - the :math:`(x,y)` coordinates of the graph vertices.
					* :attr:`Hamiltonian.lineX` *(array)* - the :math:`x` coordinates of the edges connecting graph vertices.
					* :attr:`Hamiltonian.lineY` *(array)* - the :math:`y` coordinates of the edges connecting graph vertices.
			
			.. important::
				* The number of nodes in the imported adjacency matrix \
				 **must** match the number of nodes the Hamiltonian object \
				 is initialized with.
				* The length of ``amp`` and ``d`` must be identical.

			.. admonition:: Fortran interface
				
				This function calls the Fortran function :f:func:`importadjtoh` and is \
				faster than :func:`importAdjToH`, but with a 'pickier' text file importer.
			"""
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
		""" Create a 2 particle infinite line Hamiltonian.

			Args:
				d (array of ints): an array containing *integers* indicating the nodes
							where diagonal defects are to be placed (e.g. ``d=[0,1,4]``).

				amp (array of floats):   an array containing *floats* indicating the diagonal defect
							amplitudes corresponding to each element in ``d`` (e.g. ``amp=[0.5,-1,4.2]``).

				interaction (float): the amplitude of interaction between the walkers \
									 when located on the same vertex.

			Returns:
				 :	:attr:`Hamiltonian.mat`, the Hamiltonian matrix.
			:rtype: :func:`petsc4py.PETSc.Mat`
			
			.. important::
				* The length of ``amp`` and ``d`` must be identical.

			.. admonition:: Fortran interface

				This function calls the Fortran function :f:func:`hamiltonian_p2_line`.
			"""
		try:
			if self.mat.isAssembled():
				self.reinitialize()
		except: pass
		self.defectNodes = d
		self.defectAmp = amp
		# create the Hamiltonian
		Hamiltonian = _PETSc.Log.Stage('Hamiltonian')
		Hamiltonian.push()
		_ctqwmpi.hamiltonian_p2_line(self.mat.fortran,self.defectNodes,self.defectAmp,interaction,self.N)
		Hamiltonian.pop()

	def createLine3P(self,d=[0],amp=[0.],interaction=0.):
		""" Create a 3 particle infinite line Hamiltonian.

			Args:
				d (array of ints): an array containing *integers* indicating the nodes
							where diagonal defects are to be placed (e.g. ``d=[0,1,4]``).

				amp (array of floats):   an array containing *floats* indicating the diagonal defect
							amplitudes corresponding to each element in ``d`` (e.g. ``amp=[0.5,-1,4.2]``).

				interaction (float): the amplitude of interaction between the walkers \
									 when located on the same vertex.

			Returns:
				 :	:attr:`Hamiltonian.mat`, the Hamiltonian matrix.
			:rtype: :func:`petsc4py.PETSc.Mat`
			
			.. important::
				* The length of ``amp`` and ``d`` must be identical.

			.. admonition:: Fortran interface
				
				This function calls the Fortran function :f:func:`hamiltonian_p3_line`.
			"""
		try:
			if self.mat.isAssembled():
				self.reinitialize()
		except: pass
		self.defectNodes = d
		self.defectAmp = amp
		# create the Hamiltonian
		Hamiltonian = _PETSc.Log.Stage('Hamiltonian')
		Hamiltonian.push()
		_ctqwmpi.hamiltonian_p3_line(self.mat.fortran,self.defectNodes,self.defectAmp,interaction,self.N)
		Hamiltonian.pop()
	
	def createLine(self,d=[0],amp=[0.]):
		""" Create a 1 particle infinite line Hamiltonian.

			Args:
				d (array of ints): an array containing *integers* indicating the nodes
							where diagonal defects are to be placed (e.g. ``d=[0,1,4]``).

				amp (array of floats):   an array containing *floats* indicating the diagonal defect
							amplitudes corresponding to each element in ``d`` (e.g. ``amp=[0.5,-1,4.2]``).

				interaction (float): the amplitude of interaction between the walkers \
									 when located on the same vertex.

			Returns:
				 :	:attr:`Hamiltonian.mat`, the Hamiltonian matrix.
			:rtype: :func:`petsc4py.PETSc.Mat`
			
			.. important::
				* The length of ``amp`` and ``d`` must be identical.

			.. admonition:: Fortran interface

				This function calls the Fortran function :f:func:`hamiltonian_p1_line`.
			"""
		try:
			if self.mat.isAssembled():
				self.reinitialize()
		except: pass
		self.defectNodes = d
		self.defectAmp = amp
		# create the Hamiltonian
		Hamiltonian = _PETSc.Log.Stage('Hamiltonian')
		Hamiltonian.push()
		_ctqwmpi.hamiltonian_p1_line(self.mat.fortran,self.defectNodes,self.defectAmp,self.N)
		Hamiltonian.pop()
	
	def Emax(self,**kwargs):
		""" Returns the maximum Eigenvalue of the matrix.

			Example:

				>>> Emax = Hamiltonian.Emax()
			"""
		if self.EigSolver.Emax_val is None:
			self.EigSolver.findEmax()
		return self.EigSolver.Emax_val
	
	def Emin(self,**kwargs):
		""" Returns the minimum Eigenvalue of the matrix.

			Example:

				>>> Emin = Hamiltonian.Emin()
			"""
		if self.EigSolver.Emin_val is None:
			self.EigSolver.findEmin()
		return self.EigSolver.Emin_val
		
	def destroy(self):
		""" Destroys all associated PETSc matrices."""
		self.mat.destroy()
		try:
			self.Adj.destroy()
		except:
			pass

class nodeHandle(object):
	""" Creates a handle with the ability to store and updates node probability.

		Args:
			nodes (array of ints): the nodes to watch (e.g. ``[0,1,4]``).
			t (float): the elapsed time at handle initialization.
			psi (petsc4py.PETSc.Vec) : vector containing the global probability for particle 1. 
			psi2 (None or petsc4py.PETSc.Vec) : vector containing the global probability for particle 2.
			psi3 (None or petsc4py.PETSc.Vec) : vector containing the global probability for particle 3.

		Warning:
			Note that the handle attributes/methods are **not** collective;
			if running on multiple nodes, only *local* values will be
			stored/updated on each node.
		"""
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
		""" Update the handle with node probability for an additional timestep.

			Args:
				t (float): the elapsed time at handle update.
				psi (petsc4py.PETSc.Vec) : vector containing the global probability for particle 1. 
				psi2 (None or petsc4py.PETSc.Vec) : vector containing the global probability for particle 2.
				psi3 (None or petsc4py.PETSc.Vec) : vector containing the global probability for particle 3.
			"""        
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
		""" Get thet stored node probability and timestep array for particle :attr:`p` over all watched nodes.

			Keyword Args:
				t (None or float): the elapsed time at which to return the node probabilities. \
										If ``None``, **all** timesteps will be returned.
				p (int) : (``1``,``2``,``3``) - the particle to return probability data for.

			Returns:
				tuple of arrays : returns arrays corresponding to the times recorded (if ``t=None``) \
									and the probabilities over watched nodes.

			Example:

				>>> timeArray, probXArray, probYArray = nodeHandle.getLocalNodes()

				or

				>>> probXArray, probYArray = nodeHandle.getLocalNodes(t=4)

			"""        
		if t is not None:
			try:
				indt = self._time.index(t)
			except ValueError:
				if self.rank == 0:  print '\nERROR: time {} was not handled'.format(t)
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
		""" Get thet stored probability and timestep array for particle :attr:`p` at a specified node.

			Args:
				node (int) : the node to retrieve data from.

			Keyword Args:
				p (int) : (``1``,``2``,``3``) - the particle to return probability data for.

			Returns:
				tuple of arrays : returns arrays corresponding to the times recorded (if ``t=None``) \
									and the probabilities over the specified node.

			Example:

					>>> timeArray, probXArray, probYArray = nodeHandle.getLocalNodes(3)

			"""   
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
	""" Creates a handle with the ability to store and updates entanglement.

		Args:
			t (float): the elapsed time at handle initialization.
			psi (petsc4py.PETSc.Vec) : vector containing the statespace.
			N (int) : the size of the quantum walker's position space (i.e. number \
						of graph vertices).

		Keyword Args:
			: EigSolver keywords can also be passed;
				for more details of the available EigSolver properties,
				see :class:`EigSolver`.

				.. important::
					The Eigsolver *must* be set to lapack, since **all**
					the eigenvalues must be found in the calculation of entanglement.

		.. admonition:: Fortran interface
		
			This function calls the Fortran function :f:func:`entanglement`.
		"""
	def __init__(self,t,psi,N,**kwargs):
		self.__default = {'esolver' : 'lapack',
					  'workType': 'null',
					  'workSize': '35',
					  'tol'     : 0.,
					  'maxIt'   : 0,
					  'verbose' : False,
					  'p'       : 2}

		self._rank = _PETSc.COMM_WORLD.Get_rank()
		self.time = [t]
		self.N = N

		for key,default in self.__default.iteritems():
			setattr(self, key, kwargs.get(key,default))

		entanglementS = _PETSc.Log.Stage('Entanglement')
		entanglementS.push()
		if self.p==2:
			entInit, ierr = _ctqwmpi.entanglement(psi.fortran,self.N,
				self.esolver,self.workType,self.workSize,self.tol,self.maxIt,self.verbose)
		entanglementS.pop()

		self.entanglement = [entInit]

	def update(self,t,psi):
		""" Creates a handle with the ability to store and updates entanglement.

			Args:
				t (float): the timestep of the handle update.
				psi (petsc4py.PETSc.Vec) : vector containing the statespace.

			.. admonition:: Fortran interface
			
				This function calls the Fortran function :f:func:`entanglement`.
			"""
		self.time.append(t)

		entanglementS = _PETSc.Log.Stage('Entanglement')
		entanglementS.push()
		if self.p==2:
			entUpdate, ierr = _ctqwmpi.entanglement(psi.fortran,self.N,
					self.esolver,self.workType,self.workSize,self.tol,self.maxIt,self.verbose)
		entanglementS.pop()

		self.entanglement.append(entUpdate)

	def getEntanglement(self,t=None):
		""" Returns the stored entanglement.

			Keyword Args:
				t (None or float): the timestep to retrieve the entanglement. \
									If ``None``, the the entanglement at **all** \
									recorded timesteps is returned.
			Returns:
				tuple of arrays : returns arrays corresponding to the times recorded (if ``t=None``) \
									and the entanglement over the recorded times.

			Example:

				>>> timeArray, entArray = entanglementHandle.getEntanglement()

				or

				>>> entArray = entanglementHandle.getEntanglement(t=4)
			"""
		if t is not None:
			try:
				indt = self._time.index(t)
			except ValueError:
				if self._rank == 0: print '\nERROR: time {} was not handled'.format(t)
				return

			return self.entanglement[indt]

		else:
			return _np.array(self.time), _np.array(self.entanglement)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------- 1 particle CTQW   ---------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class QuantumWalkP1(object):
	""" Methods providing a very basic framework of a 1 particle CTQW.

		.. caution::
			This object constructs some of the *framework* of a 1 particle
			quantum walk, but is not feature complete. For a feature complete
			construction, see :func:`pyCTQW.MPI.Graph` or :func:`pyCTQW.MPI.Line`
			for which this is a base.
		"""
	def __init__(self,N):
		self.rank = _PETSc.Comm.Get_rank(_PETSc.COMM_WORLD)
		self.N = N
		self.t = 0
		self._timestep = False
		
		# define vectors

		#: Initial state vector
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
		""" Imports the initial state of the quantum walk from a file.

			Args:
				filename (str): path to the file containing the input state.
				filetype (str): the filetype of the imported adjacency matrix.

									* ``'txt'`` - an :math:`N` element column vector in text format.
									* ``'bin'`` - an :math:`N` element PETSc binary vector.
			Returns:
				: this creates a PETSc vector containing the initial state, \
											accessed via :attr:`Graph.psi0`.
			:rtype: :func:`petsc4py.PETSc.Vec`

			Warning:
				The number of elements in the imported input state vector
				**must** match the number of nodes the 1P CTQW object
				is initialized with. 

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
		_ctqwmpi.marginal1(vec.fortran,self.prob.fortran,self.N)
		Marginal.pop()

	def watch(self,nodes):
		""" Creates a handle that watches node probability during propagation.

			Args:
				nodes (array of ints): the nodes to watch (e.g. ``[0,1,4]``).

			Returns:
				:	creates a handle that can be
					accessed to retrieve node probabilities for various :math:`t`
			:rtype: :func:`pyCTQW.MPI.ctqw.nodeHandle`

			Example:
				>>> walk.watch([0,1,2,3,4])
				>>> walk.propagate(5.,method='chebyshev')
				>>> timeArray, probArray = walk.handle.getLocalNodes()

			Warning:
				Note that `walk.handle` attributes are **not** collective;
				if running on multiple nodes, only *local* values will be
				returned.
			"""             
		self.handle = nodeHandle(nodes,self.t,self.prob)
		
	def propagate(self,t,method='chebyshev',**kwargs):
		""" Propagates the quantum walk for time t.

			Args:
				t (float): the timestep over which to propagate the state ``walk.psi0`` via the relationship
					:math:`\\left|\\psi\\right\\rangle = e^{-iHt}\\left|\\psi0\\right\\rangle`.

				method (str):  the propagation algorithm to use
								(``'chebyshev'`` *(default)* or ``'krylov'``).

			Keyword Args:
				esolver (str): the eigensolver algorithm to use. 

							* ``'krylovschur'`` *(default)* - Krylov-Schur
							* ``'arnoldi'`` - Arnoldi Method
							* ``'lanczos'`` - Lanczos Method
							* ``'power'`` - Power Iteration/Rayleigh Quotient Iteration
							* ``'gd'`` - Generalized Davidson
							* ``'jd'`` - Jacobi-Davidson,
							* ``'lapack'`` - Uses LAPACK eigenvalue solver routines
							* ``'arpack'`` - *only available if SLEPc is\
												compiled with ARPACK linking*

				workType (str):    can be used to set the eigensolver worktype \
								(either ``'ncv'`` or ``'mpd'``). The default \
								is to let SLEPc decide.

				workSize (int):    sets the work size **if** :attr:`workType` is set.

				tolIn (float): tolerance of the eigenvalue solver \
								(*default* ``0.`` (SLEPc decides)).
				
				maxIt (int):   maximum number of iterations of the eigenvalue solver \
								(*default* ``0`` (SLEPc decides)).
				
				verbose (bool): if ``True``, writes eigensolver information to the console

				emax_estimate (float): used to override the calculation \
										of the graphs maximum eigenvalue.

				emin_estimate (float): used to override the calculation \
										of the graphs minimum eigenvalue. For a *finite* graph, \
										``emin_estimate=0`` is set by default.

					.. caution::
						* If supplied, the value of :attr:`emax_estimate` (:math:`\hat{\lambda}_{\max}`) \
						**must** satisfy :math:`\hat{\lambda}_{\max}\geq\lambda_{\max}`, \
						where :math:`\lambda_{\max}` is the actual maximum eigenvalue of the graph.
						* Similarly, the value of :attr:`emin_estimate` (:math:`\hat{\lambda}_{\min}`) \
						**must** satisfy :math:`\hat{\lambda}_{\min}\leq\lambda_{\min}`, \
						where :math:`\lambda_{\min}` is the actual minimum eigenvalue of the graph.
						* The greater the difference value of :math:`|\hat{\lambda}_{\max} \
						-\lambda_{\max}|` and :math:`|\hat{\lambda}_{\min} \
						-\lambda_{\min}|`, the longer the convergence \
						time of the ``chebyshev`` propagator.

					Note:
						* The keyword arguments properties only apply if ``propagator='chebyshev'``
						* For more information regarding these properties,refer to \
						the `SLEPc documentation\
						<http://www.grycap.upv.es/slepc/documentation/manual.htm>`_

			Returns:
				:	this creates a PETSc vector containing the
					propagated state, accessed via the attribute
					:attr:`psi`, as well as a PETSc vector containing
					the marginal probabilities, :attr:`prob`.
			:rtype: :func:`petsc4py.PETSc.Vec`

			Note:
				The EigSolver properties only take effect if
				`method='chebyshev'`.

			.. admonition:: Fortran interface
			
				This function calls the Fortran function :f:func:`expm`. for the Krylov
				method, and the function :f:func:`qw_cheby` for the Chebyshev method.
			"""             
		if self._timestep:
			self.t = t + self.t
		else:
			self.t = t

		self.EigSolver.setEigSolver(**kwargs)
		
		if method=='krylov':
			# SLEPc matrix exponential
			expmS = _PETSc.Log.Stage('SLEPc expm'); expmS.push()
			_ctqwmpi.qw_krylov(self.H.mat.fortran,t,self.psi0.fortran,self.psi.fortran)
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
		""" Creates a plot of the node probablities over time.

			Args:
				filename (str): the absolute/relative path to the desired output file.

			.. important::
				:func:`watch` **must** be called prior to propagation in order for\
				the probabilities at the specified nodes to be stored.

			Note:
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
		""" Exports the final state of the quantum walk to a file.

			Args:
				filename (str): path to desired output file.
				filetype (str): the filetype of the exported adjacency matrix.

								* ``'txt'`` - an :math:`N` element column vector in text format.
								* ``'bin'`` - an :math:`N` element PETSc binary vector.
			"""     
		_io.exportVec(self.psi,filename,filetype)

	def psiToInit(self):
		""" Copies the state :attr:`self.psi` to :attr:`self.psi0`.

			Example:
				This can be used to iterate the propagation
				over discrete timesteps: ::

					for i in range(10):
						walk.propagate(0.01,method='chebyshev')
						walk.psiToInit()
			"""             
		self.psi0 = self.psi
		self._timestep = True
	
	def destroy(self):
		""" Destroys the 1 particle quantum walk object, and all \
			associated PETSc matrices/vectors.
			""" 
		self.H.destroy()
		self.psi.destroy()
		self.psi0.destroy()
		self.prob.destroy()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------- 2 particle CTQW   ---------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class QuantumWalkP2(object):
	""" Methods providing a very basic framework of a 2 particle CTQW.

		.. caution::
			This object constructs some of the *framework* of a 2 particle
			quantum walk, but is not feature complete. For a feature complete
			construction, see :func:`pyCTQW.MPI.Graph2P` or :func:`pyCTQW.MPI.Line2P`
			for which this is a base.
		"""
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
		""" Imports the initial state of the quantum walk from a file.

			Args:
				filename (str): path to the file containing the input state.
				filetype (str): the filetype of the imported adjacency matrix.

									* ``'txt'`` - an :math:`N\\times N` *2D array* in text format.
									* ``'bin'`` - an :math:`N^2` element PETSc binary *vector*.
			Returns:
				: this creates a PETSc vector containing the initial statespace, \
					accessed via :attr:`self.psi0`.
			:rtype: :func:`petsc4py.PETSc.Vec`

			Warning:
				The number of elements in the imported input statespace
				**must** equal :math:`N^2` (either as a :math:`N\\times N`
				text file or a :math:`N^2` PETSc binary vector),
				where :math:`N` is the number of nodes the 2P CTQW object
				is initialized with. 

			"""
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
		_ctqwmpi.marginal2(vec.fortran,self.psiX.fortran,'x',self.N)
		_ctqwmpi.marginal2(vec.fortran,self.psiY.fortran,'y',self.N)
		Marginal.pop()

	def watch(self,nodes,watchtype='prob',**kwargs):
		""" Creates a handle that watches either marginal probability \
			or walker entanglement during propagation.

			Args:
				nodes (array of ints): the nodes to watch marginal probability (e.g. ``[0,1,4]``).
				watchtype (str):	(`'prob'` , `'entanglement'`).
									
									the type of watch handle to produce.

			Keyword Args:
				: If ``watchtype='entanglement'``, EigSolver keywords can also be passed;
					for more details of the available EigSolver properties,
					see :func:`propagate`.

			Returns:
				:	* if ``watchtype='prob'``, creates a handle that can be \
					accessed to retrieve marginal node probabilities for various :math:`t`
					* if ``watchtype='entanglment'``, creates a handle that can be \
					accessed to retrieve entanglement values for various :math:`t`
			:rtype: * if ``watchtype='prob'``: :func:`pyCTQW.MPI.ctqw.nodeHandle`
					* if ``watchtype='entanglement'``: :func:`pyCTQW.MPI.ctqw.entanglementHandle`

			Examples:
				To watch the entanglement,
					>>> walk.watch(None,watchtype='entanglement')
					>>> walk.propagate(5.,method='chebyshev')
					>>> timeArray, entArray = walk.entanglementHandle.getEntanglement()

					.. note::
						* The entanglement measure used in Von Neumann entropy, calculated \
							via :math:`S=-\\sum_{i}\\lambda_i\log_{2}\\lambda_i`, where :math:`\lambda_i` \
							are the eigenvalues of the reduced density matrix \
							:math:`\\rho_2 = \\text{Tr}_1(|\\psi(t)\\rangle\\langle\\psi(t)|)`
						* Nodes do not need to be specified, as entanglement is a global measurement.
						* As it is a global measurement, there is a large amount of node communication \
						which may increase overall program run time.

				To watch the probabilities,
					>>> walk.watch([0,1,4])
					>>> walk.propagate(2.,method='chebyshev')					
					>>> timeArray, probXArray, probYArray = walk.handle.getLocalNode(4,p=2)

					.. warning::
						Note that `walk.handle` attributes are **not** collective;
						if running on multiple nodes, only *local* values will be
						returned.
			"""      
		if watchtype == 'prob':
			self._watchProb = True
			self.handle = nodeHandle(nodes,self.t,self.psiX,psi2=self.psiY)
		elif watchtype == 'entanglement':
			self._watchEnt = True
			self.entanglementHandle = entanglementHandle(self.t,self.psi0,self.N,**kwargs)
		
	def propagate(self,t,method='chebyshev',**kwargs):
		""" Propagates the quantum walk for time t.

			Args:
				t (float): the timestep over which to propagate the state ``walk.psi0`` via the relationship
					:math:`\\left|\\psi\\right\\rangle = e^{-iHt}\\left|\\psi0\\right\\rangle`.

				method (str):  the propagation algorithm to use
								(``'chebyshev'`` *(default)* or ``'krylov'``).

			Keyword Args:
				esolver (str): the eigensolver algorithm to use. 

							* ``'krylovschur'`` *(default)* - Krylov-Schur
							* ``'arnoldi'`` - Arnoldi Method
							* ``'lanczos'`` - Lanczos Method
							* ``'power'`` - Power Iteration/Rayleigh Quotient Iteration
							* ``'gd'`` - Generalized Davidson
							* ``'jd'`` - Jacobi-Davidson,
							* ``'lapack'`` - Uses LAPACK eigenvalue solver routines
							* ``'arpack'`` - *only available if SLEPc is\
												compiled with ARPACK linking*

				workType (str):    can be used to set the eigensolver worktype \
								(either ``'ncv'`` or ``'mpd'``). The default \
								is to let SLEPc decide.

				workSize (int):    sets the work size **if** :attr:`workType` is set.

				tolIn (float): tolerance of the eigenvalue solver \
								(*default* ``0.`` (SLEPc decides)).
				
				maxIt (int):   maximum number of iterations of the eigenvalue solver \
								(*default* ``0`` (SLEPc decides)).
				
				verbose (bool): if ``True``, writes eigensolver information to the console

				emax_estimate (float): used to override the calculation \
										of the graphs maximum eigenvalue.

				emin_estimate (float): used to override the calculation \
										of the graphs minimum eigenvalue. For a *finite* graph, \
										``emin_estimate=0`` is set by default.

					.. caution::
						* If supplied, the value of :attr:`emax_estimate` (:math:`\hat{\lambda}_{\max}`) \
						**must** satisfy :math:`\hat{\lambda}_{\max}\geq\lambda_{\max}`, \
						where :math:`\lambda_{\max}` is the actual maximum eigenvalue of the graph.
						* Similarly, the value of :attr:`emin_estimate` (:math:`\hat{\lambda}_{\min}`) \
						**must** satisfy :math:`\hat{\lambda}_{\min}\leq\lambda_{\min}`, \
						where :math:`\lambda_{\min}` is the actual minimum eigenvalue of the graph.
						* The greater the difference value of :math:`|\hat{\lambda}_{\max} \
						-\lambda_{\max}|` and :math:`|\hat{\lambda}_{\min} \
						-\lambda_{\min}|`, the longer the convergence \
						time of the ``chebyshev`` propagator.

					Note:
						* The keyword arguments properties only apply if ``propagator='chebyshev'``
						* For more information regarding these properties,refer to \
						the `SLEPc documentation\
						<http://www.grycap.upv.es/slepc/documentation/manual.htm>`_

			Returns:
				:	this creates a PETSc vector containing the
					propagated state, accessed via the attribute
					:attr:`psi`, as well as a PETSc vector containing
					the marginal probabilities, :attr:`prob`.
			:rtype: :func:`petsc4py.PETSc.Vec`

			Note:
				The EigSolver properties only take effect if
				`method='chebyshev'`.

			.. admonition:: Fortran interface
			
				This function calls the Fortran function :f:func:`expm`. for the Krylov
				method, and the function :f:func:`qw_cheby` for the Chebyshev method.
			"""             
		if self.timestep:
			self.t = t + self.t
		else:
			self.t = t

		self.EigSolver.setEigSolver(**kwargs)
		
		if method=='krylov':
			# SLEPc matrix exponential
			krylov = _PETSc.Log.Stage('SLEPc krylov'); krylov.push()
			_ctqwmpi.qw_krylov(self.H.mat.fortran,t,self.psi0.fortran,self.psi.fortran)
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
		""" Creates a plot of the entanglement over time.

			Args:
				filename (str): the absolute/relative path to the desired output file.

			.. important::
				:func:`watch` **must** be called prior to propagation in order for\
				the entanglement to be stored during propagation.

			Note:
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* Timesteps are recorded when propagations occur - to increase the number\
				of timesteps, increase the number of propagations prior to plotting.
			"""             

		comm = _PETSc.COMM_WORLD
		rank = comm.Get_rank()

		if rank == 0:
			timeArray, entArray = self.entanglementHandle.getEntanglement()
			_plot.plotEntanglement(timeArray,entArray,filename,self.initState,self.defectNodes,self.defectAmp)

	def plotNode(self,filename,node,t=None):
		""" Creates a plot of the marginal probablities on a *specified node* over time.

			Args:
				filename (str): the absolute/relative path to the desired output file.
				node (int): the node to plot.

			.. important::
				:func:`watch` **must** be called prior to propagation in order for\
				the probabilities at the specified nodes to be stored.

			Note:
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* Timesteps are recorded when propagations occur - to increase the number\
				of timesteps, increase the number of propagations.
				* See :func:`plotNodes` to plot multiple nodes.

			"""             

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
		""" Creates a plot of the node probablities over time.

			Args:
				filename (str): the absolute/relative path to the desired output file.
				p (int): (1|2) - choose whether to plot the marginal probability of particle 1 or 2.

			.. important::
				:func:`watch` **must** be called prior to propagation in order for\
				the probabilities at the specified nodes to be stored.

			Note:
				* if multiple nodes are watched, they will **all** be plotted.
				* if you wish to plot marginal probabilities for both particle \
					1 and 2, see :func:`plotNode`.
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* Timesteps are recorded when propagations occur - to increase the number\
				of timesteps, increase the number of propagations.

			"""             

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
		""" Exports the final state :attr:`psi` of the quantum walk to a file.

			Args:
				filename (str): path to desired output file.
				filetype (str): the filetype of the exported adjacency matrix.

								* ``'txt'`` - an :math:`N\\times N` 2D array in text format.
								* ``'bin'`` - an :math:`N^2` element PETSc binary vector.
			"""     
		if filetype == 'txt':
			_io.exportVecToMat(self.psi,filename,filetype)
		elif filetype == 'bin':
			_io.exportVec(self.psi,filename,filetype)

	def exportPartialTrace(self,filename,filetype,p=1):
		""" Exports the partial trace of the density matrix to a file.

			Args:
				filename (str): path to desired output file.
				filetype (str): the filetype of the exported adjacency matrix.

								* ``'txt'`` - an :math:`N\\times N` 2D array in text format.
								* ``'bin'`` - an :math:`N\\times N` PETSc binary matrix.

				p (int): (1|2) - the particle to trace over.

			For example, when ``p=1``, :math:`\\rho_2 = \\text{Tr}_1(\\rho_{12})`,
			where :math:`\\rho_{12} = |\\psi(t)\\rangle\\langle\\psi(t)|`, is exported.
				
			.. admonition:: Fortran interface
			
				This function calls the Fortran function :f:func:`partial_trace_mat` for ``bin``
				filetypes, and the function :f:func:`partial_trace_array` otherwise.
			""" 
		if p == 1:
			if filetype == 'bin':
				rho = _PETSc.Mat().create(_PETSc.COMM_WORLD)
				_ctqwmpi.partial_trace_mat(self.psi.fortran,rho.fortran,self.N)
				_io.exportMat(rho, filename, filetype)
				rho.destroy()
			else:
				rho = _ctqwmpi.partial_trace_array(self.psi.fortran,self.N)
				if self.rank == 0:  _np.savetxt(filename, rho.real)
		if p == 2:
			#to do
			pass

	def psiToInit(self):
		""" Copies the state :attr:`self.psi` to :attr:`self.psi0`.

			Example:
				This can be used to iterate the propagation
				over discrete timesteps: ::

					for i in range(10):
						walk.propagate(0.01,method='chebyshev')
						walk.psiToInit()
			"""             
		self.psi0 = self.psi
		self.timestep = True
	
	def destroy(self):
		""" Destroys the 2 particle quantum walk object, and all \
			associated PETSc matrices/vectors.
			""" 
		self.H.destroy()
		self.psi.destroy()
		self.psi0.destroy()
		self.psiX.destroy()
		self.psiY.destroy()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------- 3 particle CTQW   ---------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class QuantumWalkP3(object):
	""" Methods providing a very basic framework of a 3 particle CTQW.

		.. caution::
			This object constructs some of the *framework* of a 3 particle
			quantum walk, but is not feature complete. For a feature complete
			construction, see :func:`pyCTQW.MPI.Graph3P` or :func:`pyCTQW.MPI.Line3P`
			for which this is a base.
		"""
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
		""" Imports the initial state of the quantum walk from a file.

			Args:
				filename (str): path to the file containing the input state.
				filetype (str): the filetype of the imported adjacency matrix.

									* ``'txt'`` - TODO: not yet implemented! Please use ``'bin'`` for now.
									* ``'bin'`` - an :math:`N^3` element PETSc binary *vector*.
			Returns:
				: this creates a PETSc vector containing the initial statespace, \
					accessed via :attr:`self.psi0`.
			:rtype: :func:`petsc4py.PETSc.Vec`

			Warning:
				The number of elements in the imported input statespace
				**must** equal :math:`N^3` where :math:`N` is the number
				of nodes the 3P CTQW object	is initialized with. 
			"""
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
		""" Creates a handle that watches either marginal probability \
			or walker entanglement during propagation.

			Args:
				nodes (array of ints): the nodes to watch marginal probability (e.g. ``[0,1,4]``).
				watchtype (str):	(`'prob'` , `'entanglement'`).
									
									the type of watch handle to produce.

			Keyword Args:
				: If ``watchtype='entanglement'``, EigSolver keywords can also be passed;
					for more details of the available EigSolver properties,
					see :func:`propagate`.

			Returns:
				:	* if ``watchtype='prob'``, creates a handle that can be \
					accessed to retrieve marginal node probabilities for various :math:`t`
					* if ``watchtype='entanglment'``, creates a handle that can be \
					accessed to retrieve entanglement values for various :math:`t`
			:rtype: * if ``watchtype='prob'``: :func:`pyCTQW.MPI.ctqw.nodeHandle`
					* if ``watchtype='entanglement'``: :func:`pyCTQW.MPI.ctqw.entanglementHandle`

			Examples:
				To watch the entanglement,
					>>> walk.watch(None,watchtype='entanglement')
					>>> walk.propagate(5.,method='chebyshev')
					>>> timeArray, entArray = walk.entanglementHandle.getEntanglement()

					.. note::
						* The entanglement measure used in Von Neumann entropy, calculated \
							via :math:`S=-\\sum_{i}\\lambda_i\log_{2}\\lambda_i`, where :math:`\lambda_i` \
							are the eigenvalues of the reduced density matrix \
							:math:`\\rho_2 = \\text{Tr}_1(|\\psi(t)\\rangle\\langle\\psi(t)|)`
						* Nodes do not need to be specified, as entanglement is a global measurement.
						* As it is a global measurement, there is a large amount of node communication \
						which may increase overall program run time.

				To watch the probabilities,
					>>> walk.watch([0,1,4])
					>>> walk.propagate(2.,method='chebyshev')					
					>>> timeArray, probXArray, probYArray, probZArray = walk.handle.getLocalNodes(p=2)

					.. warning::
						Note that `walk.handle` attributes are **not** collective;
						if running on multiple nodes, only *local* values will be
						returned.
			"""      
		if watchtype == 'prob':
			self._watchProb = True
			self.handle = nodeHandle(nodes,self.t,self.psiX,psi2=self.psiY,psi3=self.psiZ)
		elif watchtype == 'entanglement':
			self._watchEnt = True
			self.entanglementHandle = entanglementHandle(self.t,self.psi0,self.N,**kwargs)
		
	def propagate(self,t,method='expm',**kwargs):
		""" Propagates the quantum walk for time t.

			Args:
				t (float): the timestep over which to propagate the state ``walk.psi0`` via the relationship
					:math:`\\left|\\psi\\right\\rangle = e^{-iHt}\\left|\\psi0\\right\\rangle`.

				method (str):  the propagation algorithm to use
								(``'chebyshev'`` *(default)* or ``'krylov'``).

			Keyword Args:
				esolver (str): the eigensolver algorithm to use. 

							* ``'krylovschur'`` *(default)* - Krylov-Schur
							* ``'arnoldi'`` - Arnoldi Method
							* ``'lanczos'`` - Lanczos Method
							* ``'power'`` - Power Iteration/Rayleigh Quotient Iteration
							* ``'gd'`` - Generalized Davidson
							* ``'jd'`` - Jacobi-Davidson,
							* ``'lapack'`` - Uses LAPACK eigenvalue solver routines
							* ``'arpack'`` - *only available if SLEPc is\
												compiled with ARPACK linking*

				workType (str):    can be used to set the eigensolver worktype \
								(either ``'ncv'`` or ``'mpd'``). The default \
								is to let SLEPc decide.

				workSize (int):    sets the work size **if** :attr:`workType` is set.

				tolIn (float): tolerance of the eigenvalue solver \
								(*default* ``0.`` (SLEPc decides)).
				
				maxIt (int):   maximum number of iterations of the eigenvalue solver \
								(*default* ``0`` (SLEPc decides)).
				
				verbose (bool): if ``True``, writes eigensolver information to the console

				emax_estimate (float): used to override the calculation \
										of the graphs maximum eigenvalue.

				emin_estimate (float): used to override the calculation \
										of the graphs minimum eigenvalue. For a *finite* graph, \
										``emin_estimate=0`` is set by default.

					.. caution::
						* If supplied, the value of :attr:`emax_estimate` (:math:`\hat{\lambda}_{\max}`) \
						**must** satisfy :math:`\hat{\lambda}_{\max}\geq\lambda_{\max}`, \
						where :math:`\lambda_{\max}` is the actual maximum eigenvalue of the graph.
						* Similarly, the value of :attr:`emin_estimate` (:math:`\hat{\lambda}_{\min}`) \
						**must** satisfy :math:`\hat{\lambda}_{\min}\leq\lambda_{\min}`, \
						where :math:`\lambda_{\min}` is the actual minimum eigenvalue of the graph.
						* The greater the difference value of :math:`|\hat{\lambda}_{\max} \
						-\lambda_{\max}|` and :math:`|\hat{\lambda}_{\min} \
						-\lambda_{\min}|`, the longer the convergence \
						time of the ``chebyshev`` propagator.

					Note:
						* The keyword arguments properties only apply if ``propagator='chebyshev'``
						* For more information regarding these properties,refer to \
						the `SLEPc documentation\
						<http://www.grycap.upv.es/slepc/documentation/manual.htm>`_

			Returns:
				:	this creates a PETSc vector containing the
					propagated state, accessed via the attribute
					:attr:`psi`, as well as a PETSc vector containing
					the marginal probabilities, :attr:`prob`.
			:rtype: :func:`petsc4py.PETSc.Vec`

			Note:
				The EigSolver properties only take effect if
				`method='chebyshev'`.

			.. admonition:: Fortran interface
			
				This function calls the Fortran function :f:func:`expm`. for the Krylov
				method, and the function :f:func:`qw_cheby` for the Chebyshev method.
			"""             
		if self.timestep:
			self.t = t + self.t
		else:
			self.t = t

		self.EigSolver.setEigSolver(**kwargs)
		
		if method=='krylov':
			# SLEPc matrix exponential
			krylov = _PETSc.Log.Stage('SLEPc krylov'); krylov.push()
			_ctqwmpi.qw_krylov(self.H.mat.fortran,t,self.psi0.fortran,self.psi.fortran)
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
		""" Creates a plot of the entanglement over time.

			Args:
				filename (str): the absolute/relative path to the desired output file.

			.. important::
				:func:`watch` **must** be called prior to propagation in order for\
				the entanglement to be stored during propagation.

			Note:
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* Timesteps are recorded when propagations occur - to increase the number\
				of timesteps, increase the number of propagations prior to plotting.
			"""		

		comm = _PETSc.COMM_WORLD
		rank = comm.Get_rank()

		if rank == 0:
			timeArray, entArray = self.entanglementHandle.getEntanglement()
			_plot.plotEntanglement(timeArray,entArray,filename,self.initState,self.defectNodes,self.defectAmp)

	def plotNode(self,filename,node,t=None):
		""" Creates a plot of the marginal probablities on a *specified node* over time.

			Args:
				filename (str): the absolute/relative path to the desired output file.
				node (int): the node to plot.

			.. important::
				:func:`watch` **must** be called prior to propagation in order for\
				the probabilities at the specified nodes to be stored.

			Note:
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* Timesteps are recorded when propagations occur - to increase the number\
				of timesteps, increase the number of propagations.
				* See :func:`plotNodes` to plot multiple nodes.

			"""

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
		""" Creates a plot of the node probablities over time.

			Args:
				filename (str): the absolute/relative path to the desired output file.
				p (int): (1|2|3) - choose whether to plot the marginal probability of particle 1, 2 or 3.

			.. important::
				:func:`watch` **must** be called prior to propagation in order for\
				the probabilities at the specified nodes to be stored.

			Note:
				* if multiple nodes are watched, they will **all** be plotted.
				* if you wish to plot marginal probabilities for all particles, see :func:`plotNode`.
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* Timesteps are recorded when propagations occur - to increase the number\
				of timesteps, increase the number of propagations.

			""" 

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
		""" Exports the final state :attr:`psi` of the quantum walk to a file.

			Args:
				filename (str): path to desired output file.
				filetype (str): the filetype of the exported adjacency matrix.

								* ``'txt'`` - an :math:`N^3` element column vector in text format.
								* ``'bin'`` - an :math:`N^3` element PETSc binary vector.
			"""     
		if filetype == 'txt':
			_io.exportVec(self.psi,filename,filetype)
		elif filetype == 'bin':
			_io.exportVec(self.psi,filename,filetype)

	def psiToInit(self):
		""" Copies the state :attr:`psi` to :attr:`psi0`.

			Example:
				This can be used to iterate the propagation
				over discrete timesteps: ::

					for i in range(10):
						walk.propagate(0.01,method='chebyshev')
						walk.psiToInit()
			"""             
		self.psi0 = self.psi
		self.timestep = True
	
	def destroy(self):
		""" Destroys the 3 particle quantum walk object, and all \
			associated PETSc matrices/vectors.
			""" 
		self.H.destroy()
		self.psi.destroy()
		self.psi0.destroy()
		self.psiX.destroy()
		self.psiY.destroy()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------- Arbitrary CTQW --------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Graph(QuantumWalkP1):
	""" Performs and analyses 1 particle continuous-time quantum walks on graphs

		Args:
			N (int) :	number of nodes to initialize the walker with. Nodes are \
						labeled :math:`j\\in\\{0,N-1\\}`.

		Example:
			To create a CTQW Graph object for a 10 node graph,

			>>> walk = pyCTQW.MPI.Graph(10)

		Note:
			**if filename AND filetype are provided**, this automatically
			creates a PETSc Hamiltonian matrix, neglecting the need to run
			:func:`createH`. For details on the other keyword arguments,
			see :func:`createH`.
		"""

	def __init__(self,N,filename=None,filetype=None,d=None,amp=None,layout='spring',delimiter=None):
		QuantumWalkP1.__init__(self,N)
		self.liveplot = False
		self.EigSolver.setEigSolver(emin_estimate=0.)
		
		if (filename and filetype) is not None:
			self.createH(filename,filetype,d=d,amp=amp,layout=layout,delimiter=delimiter)
		
	def createH(self,filename,filetype,d=None,amp=None,layout='spring',delimiter=None):
		""" Generate the Hamiltonian of the graph.

			Args:
				filename (str): path to the file containing the adjacency matrix of the graph

				filetype (str): the filetype of the imported adjacency matrix.

								* ``'txt'`` - an :math:`N\\times N` dense 2D array in text format.
								* ``'bin'`` - an :math:`N\\times N` PETSc binary matrix.

				d (array of ints): an array containing *integers* indicating the nodes
							where diagonal defects are to be placed (e.g. ``d=[0,1,4]``).

				amp (array of floats):   an array containing *floats* indicating the diagonal defect
							amplitudes corresponding to each element in ``d`` (e.g. ``amp=[0.5,-1,4.2]``).


				layout (str):  the format to store the position of the nodes (only used
								when running :func:`plotGraph`).

								* ``spring`` *(default)* - spring layout.
								* ``circle`` - nodes are arranged in a circle.
								* ``spectral`` - nodes are laid out according to the \
												spectrum of the graph.
								* ``random`` - nodes are arranged in a random pattern.

				delimiter (str): this is passed to `numpy.genfromtxt\
					<http://docs.scipy.org/doc/numpy/reference/generated/numpy.genfromtxt.html>`_
					in the case of strange delimiters in an imported ``txt`` file.

			Returns:
				:this creates a Hamiltonian object, accessed via the attibute
						:attr:`Graph.H`.
			:rtype: :func:`pyCTQW.MPI.ctqw.Hamiltonian`

			Note:
				This needs to be called **only** if the filename and filetype
				of the graph were not already called when the Graph object
				was initialized.
			
			.. important::
				* The number of nodes in the imported adjacency matrix \
				**must** match the number of nodes the Graph object \
				is initialized with.
				* The size of ``amp`` and ``d`` must be identical

					>>> amp = [0.5,-1.,4.2]
					>>> len(d) == len(amp)
					True
			"""
		if (d and amp) is None:
			self.defectNodes = [0]
			self.defectAmp = [0.]
		else:
			self.defectNodes = d
			self.defectAmp = amp
		self.H.importAdj(filename,filetype,d=self.defectNodes,amp=self.defectAmp,layout=layout,delimiter=delimiter)
		
	def createInitState(self,initState):
		""" Generate the initial state of the quantum walk.

			Args:
				initState (array) : an :math:`n\\times 2` array, containing the initial \
									state of the quantum walker in the format ``[[j1,amp1],[j2,amp2],...]``.

			Returns:
				:	this creates a PETSc vector containing the initial state,
					accessed via the attribute :attr:`Graph.ps0`.
			:rtype: :func:`petsc4py.PETSc.Vec`

			Example:
				For a CTQW initially located in a
				superposition of nodes 1 and 2, e.g.
				:math:`\\left|\\psi(0)\\right\\rangle = \\frac{1}{\\sqrt{2}}\\left|1\\right\\rangle \
					- \\frac{1}{\\sqrt{2}} \\left|2\\right\\rangle`,
				the initial state would be created like so:

				>>> import numpy as np
				>>> init_state = [[1,1./np.sqrt(2.)], [2,-1./np.sqrt(2.)]]
				>>> walk.createInitState(init_state)

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
		""" Creates a plot of probability vs node.

			Args:
				filename (str): the absolute/relative path to the desired output file.

			Note:
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* if :func:`propagate` has not been called, the probability of the \
				initial state is plotted (:math:`t=0`). Otherwise, the last propagated \
				time (:attr:`Graph.t`) is used.
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

	def plotGraph(self,output=None,probX=True,size=(12,8),**kwargs):
		""" Creates a plot of probability vs node superimposed on a 3D visualisation of the graph vertices.

			Args:
				output (str) :  the absolute/relative path to the desired output file.
				probX (bool) :  if set to ``True`` *(default)*, the probability is represented
								as bars placed on each graph vertex; if ``False``, the graph
								is plotted without any probability.
				size (tuple) :  ``size=(x,y)`` sets the horizontal and vertical size of the output figure.
			
			Keyword Args:
				nodesize (float) :  size of the vertices in the plot. If left blank, this is
									determined automatically.
				nodecolor (str)  :  vertex color (*default* ``'red'``).

									For more details on how to specify a color,
									see the `matplotlib documentation
									<http://matplotlib.org/api/colors_api.html#module-matplotlib.colors>`_.

				nodealpha (float)  :  value between 0 and 1 specifying the vertex opacity (*default* 0.25)

				nodetext (bool) :  if set ``True``, the vertices are labelled by number.
				nodetextcolor (str) : vertex label text color (*default* ``'black'``).
				nodetextbg (str) : vertex label background color (*default* ``'None'``).                
				ntofffset (array of floats): the :math:`(x,y,z)` vertex label offset relative \
											 to the vertex (*default* ``[0.,0.,-0.15]``).

				barscale (float) :  scaled height of the probability bars (*default* ``1``).
				barcolor (str)   :  probability bar color (*default* ``'green'``).
				baralpha (float) : value between 0 and 1 specifying the opacity (*default* 0.25)
				bartext (bool) :  if set ``True``, the probability bars are labelled with their value.
				bartextcolor (str) : probability label text color (*default* ``'black'``).
				bartextbg (str) : probability label background color (*default* ``'None'``).
				btoffset (array of floats): the :math:`(x,y,z)` probability label offset relative to the top of
											the probability bars (*default* ``[-0.025,-0.025,0.05]``)

			Note:
				* ensure a file extension is present so that file type is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* if :func:`propagate` has not been called, the probability of the \
				initial state is plotted (:math:`t=0`). Otherwise, the last propagated \
				time (``walk.t``) is used.
			"""

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
		""" Creates a *live*, updated plot of probability vs node superimposed \
			on a 3D visualisation of the graph vertices.

			Args:
				dt (str) :  the amount of time to 'sleep' before resuming the program.
							This can be used to 'slow down' propagation, providing
							time to view small changes on the graph.
				size (tuple) :  ``size=(x,y)`` sets the horizontal and vertical size of the output figure.
			
			Keyword Args:
				: For available keyword arguments, see :func:`Graph.plotGraph`.

			Note:
				* Once the live plot is no longer needed, :func:`Graph.clearLiveGraph()` \
				should be called in order to destroy the live plot object (this will *not* \
				close the plot window).
				* if :func:`propagate` has not been called, the probability of the \
				initial state is plotted (:math:`t=0`). Otherwise, the last propagated \
				time (:attr:`Graph.t`) is used.

			Warning:
				This feature attempts to uses ``matplotlib`` in a sort of `hackish` way
				and is not well tested - your mileage may vary.
			"""

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
		""" Destroys the live plot object previously created by :func:`Graph.plotLiveGraph()`.
			"""
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
	""" Performs and analyses 2 particle continuous-time quantum walks on graphs

		Args:
			N (int) :	number of nodes to initialize the walker with. Nodes are \
						labeled :math:`j\\in\\{0,N-1\\}`.

		Example:
			To create a 2P CTQW Graph object for a 10 node graph,

			>>> walk = pyCTQW.MPI.Graph2P(10)

		Note:
			**if filename AND filetype are provided**, this automatically
			creates a PETSc Hamiltonian matrix, neglecting the need to run
			:func:`createH`. For details on the other keyword arguments,
			see :func:`createH`.
		"""
	def __init__(self,N,filename=None,filetype=None,d=None,amp=None,interaction=0.):
		QuantumWalkP2.__init__(self,N)
		self.liveplot = False
		self.EigSolver.setEigSolver(emin_estimate=0.)
		
		if (filename and filetype) is not None:
			self.createH(filename,filetype,d=d,amp=amp,interaction=interaction)
		
	def createH(self,filename,filetype,d=None,amp=None,layout='spring',delimiter=None,interaction=0.):
		""" Generate the Hamiltonian of the graph.

			Args:
				filename (str): path to the file containing the adjacency matrix of the graph

				filetype (str): the filetype of the imported adjacency matrix.

								* ``'txt'`` - an :math:`N\\times N` dense 2D array in text format.
								* ``'bin'`` - an :math:`N\\times N` PETSc binary matrix.

				d (array of ints): an array containing *integers* indicating the nodes
							where diagonal defects are to be placed (e.g. ``d=[0,1,4]``).

				amp (array of floats):   an array containing *floats* indicating the diagonal defect
							amplitudes corresponding to each element in ``d`` (e.g. ``amp=[0.5,-1,4.2]``).


				layout (str):  the format to store the position of the nodes (only used
								when running :func:`plotGraph`).

								* ``spring`` *(default)* - spring layout.
								* ``circle`` - nodes are arranged in a circle.
								* ``spectral`` - nodes are laid out according to the \
												spectrum of the graph.
								* ``random`` - nodes are arranged in a random pattern.

				delimiter (str): this is passed to `numpy.genfromtxt\
					<http://docs.scipy.org/doc/numpy/reference/generated/numpy.genfromtxt.html>`_
					in the case of strange delimiters in an imported ``txt`` file.

				interaction (float): the amplitude of interaction between the two walkers \
									 when located on the same vertex.

			Returns:
				:this creates a Hamiltonian object, accessed via the attibute
						:attr:`Graph2P.H`.
			:rtype: :func:`pyCTQW.MPI.ctqw.Hamiltonian`

			Note:
				This needs to be called **only** if the filename and filetype
				of the graph were not already called when the Graph object
				was initialized.
			
			.. important::
				* The number of nodes in the imported adjacency matrix \
				**must** match the number of nodes the :class:`Graph2P` object \
				is initialized with.
				* The size of ``amp`` and ``d`` must be identical

					>>> amp = [0.5,-1.,4.2]
					>>> len(d) == len(amp)
					True
			"""
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
		#   d=self.defectNodes,amp=self.defectAmp,layout=layout)
		
	def createInitState(self,initState):
		""" Generate the initial state of the quantum walk.

			Args:
				initState (array) : an :math:`n\\times 3` array, containing the initial \
									state of the quantum walker in the format ``[[x1,y1,amp1],[x2,y2,amp2],...]``.

			Returns:
				:	this creates a PETSc vector containing the initial state,
					accessed via the attribute :attr:`Graph2P.ps0`.
			:rtype: :func:`petsc4py.PETSc.Vec`

			Example:
				For a CTQW initially located in state
				:math:`\\left|\\psi(0)\\right\\rangle = \\frac{1}{\\sqrt{2}}\\left|0\\right\\rangle\\otimes\\left|1\\right\\rangle \
					- \\frac{1}{\\sqrt{2}} \\left|1\\right\\rangle\\otimes\\left|0\\right\\rangle`,
				the initial state would be created like so:

				>>> import numpy as np
				>>> init_state = [[0,1,1./np.sqrt(2.)], [1,0,-1./np.sqrt(2.)]]
				>>> walk.createInitState(init_state)

			.. admonition:: Fortran interface
			
				This function calls the Fortran function :f:func:`p2_init`.

			"""     
		self.initState = _np.vstack([_np.array(initState).T[0]-self.N/2+1,
						_np.array(initState).T[1]-self.N/2+1,_np.array(initState).T[2]]).T.tolist()
	
		# create the inital stage
		initStateS = _PETSc.Log.Stage('initState')
		initStateS.push()
		_ctqwmpi.p2_init(self.psi0.fortran,self.initState,self.N)
		initStateS.pop()
		
	def plot(self,filename):
		""" Creates a plot of probability vs node.

			Args:
				filename (str): the absolute/relative path to the desired output file.

			Note:
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* if :func:`propagate` has not been called, the probability of the \
				initial state is plotted (:math:`t=0`). Otherwise, the last propagated \
				time (:attr:`Graph2P.t`) is used.
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
			initstateLabels.append([sum(pair) for pair in zip(self.initState[i], [self.N/2-1,self.N/2-1,0])])

		plotStage = _PETSc.Log.Stage('Plotting'); plotStage.push()      
		_plot.plot2P(_np.arange(self.N),self.psiX,self.psiY,filename,self.t,initstateLabels,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()

	def plotGraph(self,size=(12,8),probX=True,probY=True,output=None,**kwargs):
		""" Creates a plot of probability vs node superimposed on a 3D visualisation of the graph vertices.

			Args:
				output (str) :  the absolute/relative path to the desired output file.
				probX (bool) :  if set to ``True`` *(default)*, the particle 1 marginale probability is represented
								as bars placed on each graph vertex.
				probY (bool) :  if set to ``True`` *(default)*, the particle 2 marginal probability is represented
								as bars placed below each graph vertex.
				size (tuple) :  ``size=(x,y)`` sets the horizontal and vertical size of the output figure.
			
			Keyword Args:
				nodesize (float) :  size of the vertices in the plot. If left blank, this is
									determined automatically.
				nodecolor (str)  :  vertex color (*default* ``'red'``).

									For more details on how to specify a color,
									see the `matplotlib documentation
									<http://matplotlib.org/api/colors_api.html#module-matplotlib.colors>`_.

				nodealpha (float)  :  value between 0 and 1 specifying the vertex opacity (*default* 0.25)

				nodetext (bool) :  if set ``True``, the vertices are labelled by number.
				nodetextcolor (str) : vertex label text color (*default* ``'black'``).
				nodetextbg (str) : vertex label background color (*default* ``'None'``).                
				ntofffset (array of floats): the :math:`(x,y,z)` vertex label offset relative \
											 to the vertex (*default* ``[0.,0.,-0.15]``).

				barscaleP (float) :  scaled height of the probability bars (*default* ``1``).
				barcolorP (str)   :  probability bar color (*default* ``'green'``).
				baralphaP (float) : value between 0 and 1 specifying the opacity (*default* 0.25)
				bartext (bool) :  if set ``True``, the probability bars are labelled with their value.
				bartextcolorP (str) : probability label text color (*default* ``'black'``).
				bartextbgP (str) : probability label background color (*default* ``'None'``).
				btoffsetP (array of floats): the :math:`(x,y,z)` probability label offset relative to the top of \
											the probability bars (*default* ``[-0.025,-0.025,+-0.05]``)

					.. important:
						where a keyword argument ends with ``P`` above, substitute ``P=(1|2)`` to
						access the probability bar property for particle 1 and 2 respectively.

			Note:
				* ensure a file extension is present so that file type is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* if :func:`propagate` has not been called, the probability of the \
				initial state is plotted (:math:`t=0`). Otherwise, the last propagated \
				time (``walk.t``) is used.
			"""

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
		""" Creates a *live*, updated plot of probability vs node superimposed \
			on a 3D visualisation of the graph vertices.

			Args:
				dt (str) :  the amount of time to 'sleep' before resuming the program.
							This can be used to 'slow down' propagation, providing
							time to view small changes on the graph.
				size (tuple) :  ``size=(x,y)`` sets the horizontal and vertical size of the output figure.
			
			Keyword Args:
				: For available keyword arguments, see :func:`Graph2P.plotGraph`.

			Note:
				* Once the live plot is no longer needed, :func:`Graph.clearLiveGraph()` \
				should be called in order to destroy the live plot object (this will *not* \
				close the plot window).
				* if :func:`propagate` has not been called, the probability of the \
				initial state is plotted (:math:`t=0`). Otherwise, the last propagated \
				time (:attr:`Graph2P.t`) is used.

			Warning:
				This feature attempts to uses ``matplotlib`` in a sort of `hackish` way
				and is not well tested - your mileage may vary.
			"""

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
		""" Destroys the live plot object previously created by :func:`Graph2P.plotLiveGraph()`.
			"""
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
	""" Performs and analyses 3 particle continuous-time quantum walks on graphs

		Args:
			N (int) :	number of nodes to initialize the walker with. Nodes are \
						labeled :math:`j\\in\\{0,N-1\\}`.

		Example:
			To create a 3P CTQW Graph object for a 10 node graph,

			>>> walk = pyCTQW.MPI.Graph3P(10)

		Note:
			**if filename AND filetype are provided**, this automatically
			creates a PETSc Hamiltonian matrix, neglecting the need to run
			:func:`createH`. For details on the other keyword arguments,
			see :func:`createH`.
		"""
	def __init__(self,N,filename=None,filetype=None,d=None,amp=None,interaction=0.):
		QuantumWalkP3.__init__(self,N)
		self.EigSolver.setEigSolver(emin_estimate=0.)
		
		if (filename and filetype) is not None:
			self.createH(filename,filetype,d=d,amp=amp,interaction=interaction)
		
	def createH(self,filename,filetype,d=None,amp=None,layout='spring',delimiter=None,interaction=0.):
		""" Generate the Hamiltonian of the graph.

			Args:
				filename (str): path to the file containing the adjacency matrix of the graph

				filetype (str): the filetype of the imported adjacency matrix.

								* ``'txt'`` - an :math:`N\\times N` dense 2D array in text format.
								* ``'bin'`` - an :math:`N\\times N` PETSc binary matrix.

				d (array of ints): an array containing *integers* indicating the nodes
							where diagonal defects are to be placed (e.g. ``d=[0,1,4]``).

				amp (array of floats):   an array containing *floats* indicating the diagonal defect
							amplitudes corresponding to each element in ``d`` (e.g. ``amp=[0.5,-1,4.2]``).


				layout (str):  the format to store the position of the nodes (only used
								when running :func:`plotGraph`).

								* ``spring`` *(default)* - spring layout.
								* ``circle`` - nodes are arranged in a circle.
								* ``spectral`` - nodes are laid out according to the \
												spectrum of the graph.
								* ``random`` - nodes are arranged in a random pattern.

				delimiter (str): this is passed to `numpy.genfromtxt\
					<http://docs.scipy.org/doc/numpy/reference/generated/numpy.genfromtxt.html>`_
					in the case of strange delimiters in an imported ``txt`` file.

				interaction (float): the amplitude of interaction between the two walkers \
									 when located on the same vertex.

			Returns:
				:this creates a Hamiltonian object, accessed via the attibute
						:attr:`Graph3P.H`.
			:rtype: :func:`pyCTQW.MPI.ctqw.Hamiltonian`

			Note:
				This needs to be called **only** if the filename and filetype
				of the graph were not already called when the Graph object
				was initialized.
			
			.. important::
				* The number of nodes in the imported adjacency matrix \
				**must** match the number of nodes the :class:`Graph2P` object \
				is initialized with.
				* The size of ``amp`` and ``d`` must be identical

					>>> amp = [0.5,-1.,4.2]
					>>> len(d) == len(amp)
					True
			"""		
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
		#   d=self.defectNodes,amp=self.defectAmp,layout=layout)
		
	def createInitState(self,initState):
		""" Generate the initial state of the quantum walk.

			Args:
				initState (array) : an :math:`n\\times 4` array, containing the initial \
									state of the quantum walker in the format ``[[x1,y1,z1,amp1],[x2,y2,z2,amp2],...]``.

			Returns:
				:	this creates a PETSc vector containing the initial state,
					accessed via the attribute :attr:`Graph3P.ps0`.
			:rtype: :func:`petsc4py.PETSc.Vec`

			Example:
				For a CTQW initially located in state
				:math:`\\left|\\psi(0)\\right\\rangle = \\frac{1}{\\sqrt{2}}\\left|0\\right\\rangle \
					\\otimes\\left|1\\right\\rangle \\otimes\\left|5\\right\\rangle \
					- \\frac{i}{\\sqrt{2}} \\left|1\\right\\rangle\\otimes\\left|0\\right\\rangle \
					\\otimes\\left|2\\right\\rangle`,
				the initial state would be created like so:

				>>> import numpy as np
				>>> init_state = [[0,1,5,1./np.sqrt(2.)], [1,0,2,-1j./np.sqrt(2.)]]
				>>> walk.createInitState(init_state)

			.. admonition:: Fortran interface
			
				This function calls the Fortran function :f:func:`p3_init`.

			"""     
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
		""" Creates a plot of probability vs node.

			Args:
				filename (str): the absolute/relative path to the desired output file.

			Note:
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* if :func:`propagate` has not been called, the probability of the \
				initial state is plotted (:math:`t=0`). Otherwise, the last propagated \
				time (:attr:`Graph3P.t`) is used.
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
		_plot.plot3P(_np.arange(self.N),self.psiX,self.psiY,self.psiZ,filename,self.t,self.initState,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()
		
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------- 1 particle CTQW on a line -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Line(QuantumWalkP1):
	""" Performs and analyses 1 particle continuous-time \
		quantum walks on an infinite line

		Args:
			N (int) : an **even** number of nodes to initialize the walker with. Nodes are \
						labeled :math:`j\\in\\{1-N/2,N/2\\}`.

		Example:
			To create a CTQW Line object for a 10 node line,

				>>> walk = pyCTQW.MPI.Line(10)

		"""
	def __init__(self,N):
		QuantumWalkP1.__init__(self,N)
		
	def createH(self,d=None,amp=None):
		""" Generate the Hamiltonian of the graph.

			Args:
				d (array of ints): an array containing *integers* indicating the nodes
							where diagonal defects are to be placed (e.g. ``d=[0,1,4]``).

				amp (array of floats):   an array containing *floats* indicating the diagonal defect
							amplitudes corresponding to each element in ``d`` (e.g. ``amp=[0.5,-1,4.2]``).

			Returns:
				: this creates a Hamiltonian matrix, accessed via the attribute
					:attr:`Line.H`.
			:rtype: :func:`pyCTQW.MPI.ctqw.Hamiltonian`

			Warning:
				The size of ``amp`` and ``d`` must be identical

					>>> amp = [0.5,-1.,4.2]
					>>> len(d) == len(amp)
					True			
			"""
		if (d and amp) is None:
			self.defectNodes = [0]
			self.defectAmp = [0.]
		else:
			self.defectNodes = d
			self.defectAmp = amp
		self.H.createLine(d=self.defectNodes,amp=self.defectAmp)
		
	def createInitState(self,initState):
		""" Generate the initial state of the quantum walk.

			Args:
				initState (array) : an :math:`n\\times 2` array, containing the initial \
									state of the quantum walker in the format ``[[j1,amp1],[j2,amp2],...]``.

			Returns:
				:	this creates a PETSc vector containing the initial state,
					accessed via the attribute :attr:`Line.psi0`.
			:rtype: :func:`petsc4py.PETSc.Vec`

			Example:
				For a CTQW initially located in a
				superposition of nodes -4 and 2, e.g.
				:math:`\\left|\\psi(0)\\right\\rangle = \\frac{1}{\\sqrt{2}}\\left|-4\\right\\rangle \
					- \\frac{1}{\\sqrt{2}} \\left|2\\right\\rangle`,
				the initial state would be created like so:

				>>> import numpy as np
				>>> init_state = [[-4,1./np.sqrt(2.)], [2,-1./np.sqrt(2.)]]
				>>> walk.createInitState(init_state)

			.. admonition:: Fortran interface
			
				This function calls the Fortran function :f:func:`p1_init`.

			"""  
		self.initState = initState
		# create the inital stage
		initStateS = _PETSc.Log.Stage('initState')
		initStateS.push()
		_ctqwmpi.p1_init(self.psi0.fortran,initState,self.N)
		initStateS.pop()

	def watch(self,nodes):
		""" Creates a handle that watches node probability during propagation.

			Args:
				nodes (array of ints): the nodes to watch (e.g. ``[0,1,4]``).

			Returns:
				:	creates a handle that can be
					accessed to retrieve node probabilities for various :math:`t`
			:rtype: :func:`pyCTQW.MPI.ctqw.nodeHandle`

			Example:
				>>> walk.watch([0,1,2,3,4])
				>>> walk.propagate(5.,method='chebyshev')
				>>> timeArray, probArray = walk.handle.getLocalNodes()

			Warning:
				Note that `walk.handle` attributes are **not** collective;
				if running on multiple nodes, only *local* values will be
				returned.
			"""            
		nodes = [i+self.N/2-1 for i in nodes]
		super(Line,self).watch(nodes)
		
	def plot(self,filename):
		""" Creates a plot of probability vs node.

			Args:
				filename (str): the absolute/relative path to the desired output file.

			Note:
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
		""" Creates a plot of the node probablities over time.

			Args:
				filename (str): the absolute/relative path to the desired output file.

			.. important::
				:func:`watch` **must** be called prior to propagation in order for\
				the probabilities at the specified nodes to be stored.

			Note:
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
	""" Performs and analyses 2 particle continuous-time \
		quantum walks on an infinite line

		Args:
			N (int) : an **even** number of nodes to initialize the walker with. Nodes are \
						labeled :math:`j\\in\\{1-N/2,N/2\\}`.

		Example:
			To create a CTQW Line object for a 10 node line,

				>>> walk = pyCTQW.MPI.Line2P(10)

		"""
	def __init__(self,N,d=None,amp=None,interaction=0.):
		QuantumWalkP2.__init__(self,N)

		if (d is not None) and (amp is not None):
			self.defectNodes = d
			self.defectAmp = amp
			self.H.createLine2P(d=d,amp=amp,interaction=interaction)
		
	def createH(self,d=None,amp=None,interaction=0.):
		""" Generate the Hamiltonian of the graph.

			Args:
				d (array of ints): an array containing *integers* indicating the nodes
							where diagonal defects are to be placed (e.g. ``d=[0,1,4]``).

				amp (array of floats):   an array containing *floats* indicating the diagonal defect
							amplitudes corresponding to each element in ``d`` (e.g. ``amp=[0.5,-1,4.2]``).

				interaction (float): the amplitude of interaction between the two walkers \
									 when located on the same vertex.

			Returns:
				: this creates a Hamiltonian matrix, accessed via the attribute
					:attr:`Line2P.H`.
			:rtype: :func:`pyCTQW.MPI.ctqw.Hamiltonian`

			Warning:
				The size of ``amp`` and ``d`` must be identical

					>>> amp = [0.5,-1.,4.2]
					>>> len(d) == len(amp)
					True
			"""
		if (d and amp) is None:
			self.defectNodes = [0]
			self.defectAmp = [0.]
		else:
			self.defectNodes = d
			self.defectAmp = amp
			
		self.interaction = interaction
		self.H.createLine2P(d=self.defectNodes,amp=self.defectAmp,interaction=self.interaction)
		
	def createInitState(self,initState):
		""" Generate the initial state of the quantum walk.

			Args:
				initState (array) : an :math:`n\\times 3` array, containing the initial \
									state of the quantum walker in the format ``[[x1,y1,amp1],[x2,y2,amp2],...]``.

			Returns:
				:	this creates a PETSc vector containing the initial state,
					accessed via the attribute :attr:`Line2P.psi0`.
			:rtype: :func:`petsc4py.PETSc.Vec`

			Example:
				For a CTQW initially located in state
				:math:`\\left|\\psi(0)\\right\\rangle = \\frac{1}{\\sqrt{2}}\\left|0\\right\\rangle \
				\otimes \\left|0\\right\\rangle \
					- \\frac{i}{\\sqrt{2}} \\left|2\\right\\rangle\otimes \\left|2\\right\\rangle`,
				the initial state would be created like so:

				>>> import numpy as np
				>>> init_state = [[0,0,1./np.sqrt(2.)], [2,2,-1j./np.sqrt(2.)]]
				>>> walk.createInitState(init_state)

			.. admonition:: Fortran interface
			
				This function calls the Fortran function :f:func:`p2_init`.

			"""  
		self.initState = initState
		# create the inital stage
		initStateS = _PETSc.Log.Stage('initState')
		initStateS.push()
		_ctqwmpi.p2_init(self.psi0.fortran,initState,self.N)
		initStateS.pop()

	def watch(self,nodes,watchtype='prob',**kwargs):
		""" Creates a handle that watches either marginal probability \
			or walker entanglement during propagation.

			Args:
				nodes (array of ints): the nodes to watch marginal probability (e.g. ``[0,1,4]``).
				watchtype (str):	(`'prob'` , `'entanglement'`).
									
									the type of watch handle to produce.

			Keyword Args:
				: If ``watchtype='entanglement'``, EigSolver keywords can also be passed;
					for more details of the available EigSolver properties,
					see :func:`propagate`.

			Returns:
				:	* if ``watchtype='prob'``, creates a handle that can be \
					accessed to retrieve marginal node probabilities for various :math:`t`
					* if ``watchtype='entanglment'``, creates a handle that can be \
					accessed to retrieve entanglement values for various :math:`t`
			:rtype: * if ``watchtype='prob'``: :func:`pyCTQW.MPI.ctqw.nodeHandle`
					* if ``watchtype='entanglement'``: :func:`pyCTQW.MPI.ctqw.entanglementHandle`

			Examples:
				To watch the entanglement,
					>>> walk.watch(None,watchtype='entanglement')
					>>> walk.propagate(5.,method='chebyshev')
					>>> timeArray, entArray = walk.entanglementHandle.getEntanglement()

					.. note::
						* The entanglement measure used in Von Neumann entropy, calculated \
							via :math:`S=-\\sum_{i}\\lambda_i\log_{2}\\lambda_i`, where :math:`\lambda_i` \
							are the eigenvalues of the reduced density matrix \
							:math:`\\rho_2 = \\text{Tr}_1(|\\psi(t)\\rangle\\langle\\psi(t)|)`
						* Nodes do not need to be specified, as entanglement is a global measurement.
						* As it is a global measurement, there is a large amount of node communication \
						which may increase overall program run time.

				To watch the probabilities,
					>>> walk.watch([0,1,4])
					>>> walk.propagate(2.,method='chebyshev')					
					>>>timeArray, probXArray, probYArray = self.handle.getLocalNode(0,p=2)

					.. warning::
						Note that `walk.handle` attributes are **not** collective;
						if running on multiple nodes, only *local* values will be
						returned.
			"""      
		if watchtype == 'prob':
			nodes = [i+self.N/2-1 for i in nodes]
		super(Line2P,self).watch(nodes,watchtype=watchtype,**kwargs)
		
	def plot(self,filename):
		""" Creates a plot of probability vs node.

			Args:
				filename (str): the absolute/relative path to the desired output file.

			Note:
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* if :func:`propagate` has not been called, the probability of the \
				initial state is plotted (:math:`t=0`). Otherwise, the last propagated \
				time (:attr:`Line2P.t`) is used.
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
		_plot.plot2P(_np.arange(1-self.N/2,self.N/2+1),self.psiX,self.psiY,filename,self.t,self.initState,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()

	def plotNode(self,filename,node,t=None):
		""" Creates a plot of the marginal probablities on a *specified node* over time.

			Args:
				filename (str): the absolute/relative path to the desired output file.
				node (int): the node to plot.

			.. important::
				:func:`watch` **must** be called prior to propagation in order for\
				the probabilities at the specified nodes to be stored.

			Note:
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* Timesteps are recorded when propagations occur - to increase the number\
				of timesteps, increase the number of propagations.
				* See :func:`plotNodes` to plot multiple nodes.

			"""             
		node = node+self.N/2-1
		super(Line2P,self).plotNode(filename,node)

	def plotNodes(self,filename,p=1,t=None):
		""" Creates a plot of the node probablities over time.

			Args:
				filename (str): the absolute/relative path to the desired output file.
				p (int): (1|2) - choose whether to plot the marginal probability of particle 1 or 2.

			.. important::
				:func:`watch` **must** be called prior to propagation in order for\
				the probabilities at the specified nodes to be stored.

			Note:
				* if multiple nodes are watched, they will **all** be plotted.
				* if you wish to plot marginal probabilities for both particle \
					1 and 2, see :func:`plotNode`.
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* Timesteps are recorded when propagations occur - to increase the number\
				of timesteps, increase the number of propagations.

			""" 

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
	""" Performs and analyses 3 particle continuous-time \
		quantum walks on an infinite line

		Args:
			N (int) : an **even** number of nodes to initialize the walker with. Nodes are \
						labeled :math:`j\\in\\{1-N/2,N/2\\}`.

		Example:
			To create a CTQW Line object for a 10 node line,

				>>> walk = pyCTQW.MPI.Line3P(10)

		"""
	def __init__(self,N,d=None,amp=None,interaction=0.):
		QuantumWalkP3.__init__(self,N)

		if (d is not None) and (amp is not None):
			self.defectNodes = d
			self.defectAmp = amp
			self.H.createLine3P(d=d,amp=amp,interaction=interaction)
		
	def createH(self,d=None,amp=None,interaction=0.):
		""" Generate the Hamiltonian of the graph.

			Args:
				d (array of ints): an array containing *integers* indicating the nodes
							where diagonal defects are to be placed (e.g. ``d=[0,1,4]``).

				amp (array of floats):   an array containing *floats* indicating the diagonal defect
							amplitudes corresponding to each element in ``d`` (e.g. ``amp=[0.5,-1,4.2]``).

				interaction (float): the amplitude of interaction between the two walkers \
									 when located on the same vertex.

			Returns:
				: this creates a Hamiltonian matrix, accessed via the attribute
					:attr:`Line3P.H`.
			:rtype: :func:`pyCTQW.MPI.ctqw.Hamiltonian`

			Warning:
				The size of ``amp`` and ``d`` must be identical

					>>> amp = [0.5,-1.,4.2]
					>>> len(d) == len(amp)
					True
			"""
		if (d and amp) is None:
			self.defectNodes = [0]
			self.defectAmp = [0.]
		else:
			self.defectNodes = d
			self.defectAmp = amp

		self.interaction = interaction

		self.H.createLine3P(d=self.defectNodes,amp=self.defectAmp,interaction=self.interaction)
		
	def createInitState(self,initState):
		""" Generate the initial state of the quantum walk.

			Args:
				initState (array) : an :math:`n\\times 4` array, containing the initial \
									state of the quantum walker in the format ``[[x1,y1,z1,amp1],[x2,y2,z2,amp2],...]``.

			Returns:
				:	this creates a PETSc vector containing the initial state,
					accessed via the attribute :attr:`Line3P.ps0`.
			:rtype: :func:`petsc4py.PETSc.Vec`

			Example:
				For a CTQW initially located in state
				:math:`\\left|\\psi(0)\\right\\rangle = \\frac{1}{\\sqrt{2}}\\left|-4\\right\\rangle \
					\\otimes\\left|1\\right\\rangle \\otimes\\left|5\\right\\rangle \
					- \\frac{i}{\\sqrt{2}} \\left|1\\right\\rangle\\otimes\\left|0\\right\\rangle \
					\\otimes\\left|2\\right\\rangle`,
				the initial state would be created like so:

				>>> import numpy as np
				>>> init_state = [[-4,1,5,1./np.sqrt(2.)], [1,0,2,-1j./np.sqrt(2.)]]
				>>> walk.createInitState(init_state)

			.. admonition:: Fortran interface
			
				This function calls the Fortran function :f:func:`p3_init`.

			"""     
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
		""" Creates a handle that watches either marginal probability \
			or walker entanglement during propagation.

			Args:
				nodes (array of ints): the nodes to watch marginal probability (e.g. ``[0,1,4]``).
				watchtype (str):	(`'prob'` , `'entanglement'`).
									
									the type of watch handle to produce.

			Keyword Args:
				: If ``watchtype='entanglement'``, EigSolver keywords can also be passed;
					for more details of the available EigSolver properties,
					see :func:`propagate`.

			Returns:
				:	* if ``watchtype='prob'``, creates a handle that can be \
					accessed to retrieve marginal node probabilities for various :math:`t`
					* if ``watchtype='entanglment'``, creates a handle that can be \
					accessed to retrieve entanglement values for various :math:`t`
			:rtype: * if ``watchtype='prob'``: :func:`pyCTQW.MPI.ctqw.nodeHandle`
					* if ``watchtype='entanglement'``: :func:`pyCTQW.MPI.ctqw.entanglementHandle`

			Examples:
				To watch the entanglement,
					>>> walk.watch(None,watchtype='entanglement')
					>>> walk.propagate(5.,method='chebyshev')
					>>> timeArray, entArray = walk.entanglementHandle.getEntanglement()

					.. note::
						* The entanglement measure used in Von Neumann entropy, calculated \
							via :math:`S=-\\sum_{i}\\lambda_i\log_{2}\\lambda_i`, where :math:`\lambda_i` \
							are the eigenvalues of the reduced density matrix \
							:math:`\\rho_2 = \\text{Tr}_1(|\\psi(t)\\rangle\\langle\\psi(t)|)`
						* Nodes do not need to be specified, as entanglement is a global measurement.
						* As it is a global measurement, there is a large amount of node communication \
						which may increase overall program run time.

				To watch the probabilities,
					>>> walk.watch([0,1,4])
					>>> walk.propagate(2.,method='chebyshev')					
					>>>timeArray, probXArray, probYArray, probZArray = self.handle.getLocalNode(4,p=2)

					.. warning::
						Note that `walk.handle` attributes are **not** collective;
						if running on multiple nodes, only *local* values will be
						returned.
			"""      
		nodes = [i+self.N/2-1 for i in nodes]
		super(Line3P,self).watch(nodes,type=type)
		
	def plot(self,filename):
		""" Creates a plot of probability vs node.

			Args:
				filename (str): the absolute/relative path to the desired output file.

			Note:
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* if :func:`propagate` has not been called, the probability of the \
				initial state is plotted (:math:`t=0`). Otherwise, the last propagated \
				time (:attr:`Line3P.t`) is used.
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
			initstateLabels.append([sum(pair) for pair in zip(self.initState[i], [1-self.N/2,1-self.N/2,1-self.N/2,0])])

		plotStage = _PETSc.Log.Stage('Plotting'); plotStage.push()      
		_plot.plot3P(_np.arange(1-self.N/2,self.N/2+1),self.psiX,self.psiY,self.psiZ,filename,self.t,initstateLabels,
					self.defectNodes,self.defectAmp,self.N,self.rank)
		plotStage.pop()

	def plotNode(self,filename,node,t=None):
		""" Creates a plot of the marginal probablities on a *specified node* over time.

			Args:
				filename (str): the absolute/relative path to the desired output file.
				node (int): the node to plot.

			.. important::
				:func:`watch` **must** be called prior to propagation in order for\
				the probabilities at the specified nodes to be stored.

			Note:
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* Timesteps are recorded when propagations occur - to increase the number\
				of timesteps, increase the number of propagations.
				* See :func:`plotNodes` to plot multiple nodes.

			"""             
		node = node+self.N/2-1
		super(Line3P,self).plotNode(filename,node)

	def plotNodes(self,filename,p=1,t=None):
		""" Creates a plot of the node probablities over time.

			Args:
				filename (str): the absolute/relative path to the desired output file.
				p (int): (1|2|3) - choose whether to plot the marginal probability of particle 1, 2 or 3.

			.. important::
				:func:`watch` **must** be called prior to propagation in order for\
				the probabilities at the specified nodes to be stored.

			Note:
				* if multiple nodes are watched, they will **all** be plotted.
				* if you wish to plot marginal probabilities for all particles, see :func:`plotNode`.
				* ensure a file extension is present so that filetype is correctly set\
				(choose one of png, pdf, ps, eps or svg).
				* Timesteps are recorded when propagations occur - to increase the number\
				of timesteps, increase the number of propagations.

			""" 

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
	""" A graph isomorphism solver, containing functions for 
		creating graph certificates and checking isomorphism
		of adjacency matrices.

		>>> gi = pyCTQW.MPI.GraphISO(p=2,propagator='krylov')

		Args:
			p (int):   number of particles, 2 *(default)* or 3,
					to use in constructing the graph certificate.

			freqTol (float): the tolerance to use when constructing the
						frequency table (*default* ``1.e-2``). See also
						:py:func:`GIcert`.

			compareTol (float):    the tolerance used when comparing two Graph
							certificates (*default* ``1.e-10``). See also
							:py:func:`isomorphicQ`.

			propagator (str):  the CTQW propagator algorithm to use
							when calculating the graph certificate
							(``'chebyshev'`` *(default)* or ``'krylov'``).

		Note:
			* For ``freqTol=1.e-2``, all decimal places \
				below 0.01 are discarded from the probability \
				distribution.
			* Two isomorphic certificates satisfy \
				:math:`\max(|cert_1 - cert_2|) < compareTol`.

		"""

	def __init__(self,**kwargs):        
		self.__default = {
						'p'             : 2,
						'freqTol'       : 1.e-2,
						'compareTol'    : 1.e-10,
						'propagator'    : 'chebyshev'
						}

		self.__eigDefault = {
						'esolver'       : 'krylovschur',
						'emax_estimate' : 0.,
						'workType'      : 'null',
						'workSize'      : '35',
						'tol'           : 0.,
						'maxIt'         : 0,
						'verbose'       : False
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

			Keyword Args:
				esolver (str): the eigensolver algorithm to use. 

							* ``'krylovschur'`` *(default)* - Krylov-Schur
							* ``'arnoldi'`` - Arnoldi Method
							* ``'lanczos'`` - Lanczos Method
							* ``'power'`` - Power Iteration/Rayleigh Quotient Iteration
							* ``'gd'`` - Generalized Davidson
							* ``'jd'`` - Jacobi-Davidson,
							* ``'lapack'`` - Uses LAPACK eigenvalue solver routines
							* ``'arpack'`` - *only available if SLEPc is\
												compiled with ARPACK linking*

				workType (str):    can be used to set the eigensolver worktype
								(either ``'ncv'`` or ``'mpd'``). The default
								is to let SLEPc decide.

				workSize (int):    sets the work size **if** ``workType`` is set.

				tolIn (float): tolerance of the eigenvalue solver
								(*default* ``0.`` (SLEPc decides)).
				
				maxIt (int):   maximum number of iterations of the eigenvalue solver
								(*default* ``0`` (SLEPc decides)).
				
				verbose (bool): if ``True``, writes eigensolver information to the console

				emax_estimate (float): used to override the calculation
										of the graphs maximum eigenvalue.

			.. caution::
				* If supplied, the value of :attr:`emax_estimate`:math:`\hat{\lambda}_{\max}` \
				**must** satisfy :math:`\hat{\lambda}_{\max}\geq\lambda_{\max}`, \
				where :math:`\lambda_{\max}` is the actual maximum eigenvalue of the graph.
				* The greater the value of :math:`\hat{\lambda}_{\max} \
				-\lambda_{\max}`, the longer the convergence \
				time of the ``chebyshev`` propagator.

			Note:
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

			Args:
				adj (array or numpy.array): symmetric adjacency matrix in the form
						of a dense array of numpy array.

			Returns:
				graph certificate
			:rtype: :func:`numpy.array`

			.. admonition:: Fortran interface
			
				This function calls the Fortran function :f:func:`GraphISCert`.
			"""

		GIcertS = _PETSc.Log.Stage('GICert')
		GIcertS.push()
		cert, certSize = _ctqwmpi.GraphISCert(adj,self.p,self.freqTol,self.propagator,
			self.esolver,self.emax_estimate,self.workType,self.workSize,self.tol,self.maxIt,self.verbose)
		GIcertS.pop()

		return _np.array(cert).T[_np.lexsort(_np.array(cert)[:,0:certSize])[::-1]]

	def isomorphicQ(self,adj1,adj2):
		"""Returns ``True`` if two graphs are isomorphic.

			Args:
				adj(1|2) (array or numpy.array): \
					symmetric adjacency matrices in the form
					of a dense array of numpy array.
			
			:rtype: *bool*
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

			Args:
				folder (str):  path to a folder containing a collection of adjacency
								matrices in dense text file format.

				graphRange (array):    an array containing graph numbers to test.
										By default, all graphs in a folder are tested.

				info (bool):   if ``True``, information on each :math:`R_{ij}` comparison
								is printed to the console.

				checkSelf (bool):   if ``True``, each graph is also tested against itself.

			:rtype: *array of ints*

			Note:
				The text files must have filenames of the form \*X.txt
				where X represents a number (of any number of digits).
				These are used to order the graphs.
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
					if (self.__rank==0 and info):   print 'Adding graph ' + graph
			else:
				adj.append(_np.genfromtxt(graph))
				if (self.__rank==0 and info):   print 'Adding graph ' + graph

		NG = len(adj)
		comparisonTable = _np.zeros([NG,NG])

		for i in range(NG):
			for j in range(i if checkSelf else i+1,NG):
				if (self.__rank==0 and info):   print 'Testing graphs ' + str(i) + ',' + str(j)
				comparisonTable[i,j] = 1 if self.isomorphicQ(adj[i],adj[j]) else 0
				if (self.__rank==0 and info):   print '\tIsomorphic' if comparisonTable[i,j] == 1 else '\tNon-isomorphic'

		if checkSelf:
			return comparisonTable + comparisonTable.T - _np.diag(comparisonTable.diagonal())
		else:
			return comparisonTable + comparisonTable.T + _np.diag(_np.ones(NG))