#!/usr/bin/python
import os, errno, sys
from petsc4py import PETSc
from matplotlib import pyplot as plt
import numpy as np
import pylab as pl
import matplotlib
import fileinput


def vecToArray(obj):
	# scatter vector 'obj' to all processes
	comm = obj.getComm()
	rank = comm.getRank()
	scatter, obj0 = PETSc.Scatter.toAll(obj)
	scatter.scatter(obj, obj0, False, PETSc.Scatter.Mode.FORWARD)

	return np.asarray(obj0)

	# deallocate
	comm.barrier()
	scatter.destroy()
	obj0.destroy()
	
def vecToArray0(obj):
	# scatter vector 'obj' to process 0
	comm = obj.getComm()
	rank = comm.getRank()
	scatter, obj0 = PETSc.Scatter.toZero(obj)
	scatter.scatter(obj, obj0, False, PETSc.Scatter.Mode.FORWARD)

	if rank == 0:	return np.asarray(obj0)

	# deallocate
	comm.barrier()
	scatter.destroy()
	obj0.destroy()
	
def arrayToVec(vecArray):
	vec = PETSc.Vec().create(comm=PETSc.COMM_WORLD)
	vec.setSizes(len(vecArray))
	vec.setUp()
	(Istart,Iend) = vec.getOwnershipRange()
	return vec.createWithArray(vecArray[Istart:Iend],
			comm=PETSc.COMM_WORLD)
	vec.destroy()
	
def arrayToMat(matArray):
	try:
		import scipy.sparse as sparse
	except:
		print '\nERROR: loading matrices from txt files requires Scipy!'
		sys.exit()
						
	matSparse = sparse.csr_matrix(matArray)
				
	mat = PETSc.Mat().createAIJ(size=matSparse.shape,comm=PETSc.COMM_WORLD)
	(Istart,Iend) = mat.getOwnershipRange()
	
	ai = matSparse.indptr[Istart:Iend+1] - matSparse.indptr[Istart]
	aj = matSparse.indices[matSparse.indptr[Istart]:matSparse.indptr[Iend]]
	av = matSparse.data[matSparse.indptr[Istart]:matSparse.indptr[Iend]]
	
	mat.setValuesCSR(ai,aj,av)
	mat.assemble()
	
	return mat
	mat.destroy()

def adjToH(adj,d=[0],amp=[0.]):
	(Istart,Iend) = adj.getOwnershipRange()
	diagSum = []
	for i in range(Istart,Iend):
		diagSum.append(np.sum(adj.getRow(i)[-1]))
		for j,val in enumerate(d):
			if i==val:	diagSum[i-Istart] += amp[j]
	
	mat = PETSc.Mat().create(comm=PETSc.COMM_WORLD)
	mat.setSizes(adj.getSize())
	mat.setUp()
	
	for i in range(Istart,Iend):
		mat.setValue(i,i,diagSum[i-Istart])
	
	mat.assemble()
	mat.axpy(-1,adj)
	
	return mat
	mat.destroy()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---------------------- Vec I/O functions ---------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def exportVec(vec,filename,filetype):
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
	
	if filetype == 'txt':
		# scatter prob to process 0
		comm = vec.getComm()
		rank = comm.getRank()
		scatter, vec0 = PETSc.Scatter.toZero(vec)
		scatter.scatter(vec, vec0, False, PETSc.Scatter.Mode.FORWARD)

		# use process 0 to write to text file
		if rank == 0:
			array0 = np.asarray(vec0)
			with open(filename,'w') as f:
				for i in range(len(array0)):
					f.write('{0: .12e}\n'.format(array0[i]))

		# deallocate	
		comm.barrier()
		scatter.destroy()
		vec0.destroy()
		
	elif filetype == 'bin':
		binSave = PETSc.Viewer().createBinary(filename, 'w')
		binSave(vec)
		binSave.destroy()
	
	vec.comm.barrier()

def loadVec(filename,filetype):
	if filetype == 'txt':
		try:
			vecArray = np.loadtxt(filename,dtype=PETSc.ScalarType)
			return arrayToVec(vecArray)
		except:
			print "\nERROR: input state space file " + filename\
				+ " does not exist or is in an incorrect format"
			sys.exit()
		
	elif filetype == 'bin':
		binLoad = PETSc.Viewer().createBinary(filename, 'r')
		try:
			return PETSc.Vec().load(binLoad)
		except:
			print "\nERROR: input state space file " + filename\
				+ " does not exist or is in an incorrect format"
			sys.exit()			
		binLoad.destroy()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---------------------- Mat I/O functions ---------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def exportMat(mat,filename,filetype):
	rank = PETSc.Comm.Get_rank(PETSc.COMM_WORLD)

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
	
	if filetype == 'txt':
		txtSave = PETSc.Viewer().createASCII(filename, 'w',
			format=PETSc.Viewer.Format.ASCII_DENSE, comm=PETSc.COMM_WORLD)
		txtSave(mat)
		txtSave.destroy()
		
		if rank == 0:
			for line in fileinput.FileInput(filename,inplace=1):
				if line[2] != 't':
					line = line.replace(" i","j")
					line = line.replace(" -","-")
					line = line.replace("+-","-")
					print line,
		
	elif filetype == 'bin':
		binSave = PETSc.Viewer().createBinary(filename, 'w', comm=PETSc.COMM_WORLD)
		binSave(mat)
		binSave.destroy()
		
	mat.comm.barrier()

def loadMat(filename,filetype):
	if filetype == 'txt':
		try:
			try:
				matArray = np.loadtxt(filename,dtype=PETSc.ScalarType)
			except:
				filefix = []
				for line in fileinput.FileInput(filename,inplace=0):
					if line[2] != 't':
						line = line.replace(" i","j")
						line = line.replace(" -","-")
						line = line.replace("+-","-")
						filefix.append(line)
						
				matArray = np.loadtxt(filefix,dtype=PETSc.ScalarType)
				
			return arrayToMat(matArray)
		except:
			print "\nERROR: input state space file " + filename\
				+ " does not exist or is in an incorrect format"
			sys.exit()		
		
	elif filetype == 'bin':
		binLoad = PETSc.Viewer().createBinary(filename, 'r')
		try:
			return PETSc.Mat().load(binLoad)
		except:
			print "\nERROR: input state space file " + filename\
				+ " does not exist or is in an incorrect format"
			sys.exit()			
		binLoad.destroy()

def exportVecToMat(vec,filename,filetype):
	rank = PETSc.Comm.Get_rank(PETSc.COMM_WORLD)
	
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
	
	vecArray = vecToArray(vec)
	matArray = vecArray.reshape([np.sqrt(vecArray.size),np.sqrt(vecArray.size)])
	
	if filetype == 'txt':
		#if rank == 0:	np.savetxt(filename,matArray)
		txtSave = PETSc.Viewer().createASCII(filename, 'w',
			format=PETSc.Viewer.Format.ASCII_DENSE, comm=PETSc.COMM_WORLD)
		txtSave(arrayToMat(matArray))
		txtSave.destroy()
	
		if rank == 0:
			for line in fileinput.FileInput(filename,inplace=1):
				if line[2] != 't':
					line = line.replace(" i","j")
					line = line.replace(" -","-")
					line = line.replace("+-","-")
					print line,
		
	elif filetype == 'bin':
		binSave = PETSc.Viewer().createBinary(filename, 'w', comm=PETSc.COMM_WORLD)
		binSave(arrayToMat(matArray))
		binSave.destroy()
	vec.comm.barrier()
	
def loadMatToVec(filename,filetype):
	if filetype == 'txt':
		try:
			try:
				matArray = np.loadtxt(filename,dtype=PETSc.ScalarType)
			except:
				filefix = []
				for line in fileinput.FileInput(filename,inplace=0):
					if line[2] != 't':
						line = line.replace(" i","j")
						line = line.replace(" -","-")
						line = line.replace("+-","-")
						filefix.append(line)
						
				matArray = np.loadtxt(filefix,dtype=PETSc.ScalarType)
						
			vecArray = matArray.reshape(matArray.shape[0]**2)
			return arrayToVec(vecArray)
		except:
			print "\nERROR: input state space file " + filename\
				+ " does not exist or is in an incorrect format"
			sys.exit()		
		
	elif filetype == 'bin':
		print '\nERROR: only works for txt storage!'
		sys.exit()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------- Plotting functions -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plotGraph(prob,savefile,t,init_state,d,amp,N,rank):

	def prob_plot_p1(probVec,savefile,t,initstate,d,a,N):
		# convert vectors to arrays
		prob = np.real(np.asarray(probVec))
		# determine the plot range
		x = np.arange(0,N)
	
		# create plot
		fig = plt.figure()
		plt.plot(x, prob)
		plt.ylabel("$|\langle j|\psi\\rangle|^2$",rotation=90)
		plt.grid(b=None, which='major', axis='both', linestyle='-', alpha=0.3)
		plt.xlabel("$j$")
	
		# Plot titles
		if initstate[0]=='f':
			disp = "\nInitial state: {0}"
			IS_disp = [initstate]
		else:
			IS_disp = ()
			disp = "\nInitial state: $|\psi(0)\\rangle="

			for i in range(len(initstate)):
				IS_disp = IS_disp + ("({1: .3f})|{0}\\rangle".format(*(initstate[i])),)
	
				if i == 0:
					disp = disp + "{" + str(i) + "} "
				else:
					disp = disp + "+ {" + str(i) + "} "	
		
		plt.suptitle("CTQW probability distribution at time $t={}$".format(t))
	
		if (len(list(set(a))) == 1) and (list(set(a))[0] == 0.0):
			plt.title(disp.format(*IS_disp) + "$",
				horizontalalignment='right',multialignment='left', fontsize=11)
		else:
			def_disp = "\nDefects: $"
			for i in range(len(d)):
				def_disp += "{1: .3f}|{0}\\rangle +".format(i+1,d[i],a[i])
		
			plt.title(disp.format(*IS_disp) + "$" + def_disp[:-2] + "$",
				horizontalalignment='right',multialignment='left', fontsize=11)

		# save plot
		plt.subplots_adjust(top=0.85)
		pl.savefig(savefile)
		
	# scatter prob to process 0
	commX = prob.getComm()
	scatterX, prob0 = PETSc.Scatter.toZero(prob)
	scatterX.scatter(prob, prob0, False, PETSc.Scatter.Mode.FORWARD)
	
	# use process 0 to create the plot
	if rank==0:
		prob_plot_p1(prob0,savefile,t,init_state,d,amp,N)
	
	# deallocate	
	commX.barrier()
	scatterX.destroy()
	prob0.destroy()

def plotLine(prob,savefile,t,init_state,d,amp,N,rank):

	def prob_plot_p1(probVec,savefile,t,initstate,d,a,N):
		# convert vectors to arrays
		prob = np.real(np.asarray(probVec))
		# determine the plot range
		x = np.arange(1-N/2,N/2+1)
	
		# create plot
		fig = plt.figure()
		plt.plot(x, prob)
		plt.ylabel("$|\langle j|\psi\\rangle|^2$",rotation=90)
		plt.grid(b=None, which='major', axis='both', linestyle='-', alpha=0.3)
		plt.xlabel("$j$")
	
		# Plot titles
		if initstate[0]=='f':
			disp = "\nInitial state: {0}"
			IS_disp = [initstate]
		else:
			IS_disp = ()
			disp = "\nInitial state: $|\psi(0)\\rangle="

			for i in range(len(initstate)):
				IS_disp = IS_disp + ("({1: .3f})|{0}\\rangle".format(*(initstate[i])),)
	
				if i == 0:
					disp = disp + "{" + str(i) + "} "
				else:
					disp = disp + "+ {" + str(i) + "} "	
		
		plt.suptitle("CTQW probability distribution at time $t={}$".format(t))
	
		if (len(list(set(a))) == 1) and (list(set(a))[0] == 0.0):
			plt.title(disp.format(*IS_disp) + "$",
				horizontalalignment='right',multialignment='left', fontsize=11)
		else:
			def_disp = "\nDefects: $"
			for i in range(len(d)):
				def_disp += "{1: .3f}|{0}\\rangle +".format(i+1,d[i],a[i])
		
			plt.title(disp.format(*IS_disp) + "$" + def_disp[:-2] + "$",
				horizontalalignment='right',multialignment='left', fontsize=11)

		# save plot
		plt.subplots_adjust(top=0.85)
		pl.savefig(savefile)
		
	# scatter prob to process 0
	commX = prob.getComm()
	scatterX, prob0 = PETSc.Scatter.toZero(prob)
	scatterX.scatter(prob, prob0, False, PETSc.Scatter.Mode.FORWARD)
	
	# use process 0 to create the plot
	if rank==0:
		prob_plot_p1(prob0,savefile,t,init_state,d,amp,N)
	
	# deallocate	
	commX.barrier()
	scatterX.destroy()
	prob0.destroy()

def plotLine2P(psiX,psiY,savefile,t,init_state,d,amp,N,rank):
	
	def prob_plot_p2(psiX,psiY,savefile,t,initstate,d,a,N):
	
		# convert vectors to arrays
		probX = np.real(np.asarray(psiX))
		probY = np.real(np.asarray(psiY))
		# determine the plot range
		x = np.arange(1-N/2,N/2+1)

		# create plot
		lbl = [0,1];
	
		fig = plt.figure()
		lbl[0], = plt.plot(x, probX, 'b-')
		lbl[1], = plt.plot(x, probY, 'r--')
	
		plt.ylabel("$|\langle j|\psi\\rangle|^2$",rotation=90)
		plt.grid(b=None, which='major', axis='both', linestyle='-', alpha=0.3)
		plt.xlabel("$j$")
	
		# legend
		leg = plt.legend(lbl, ['P1','P2'], loc='center right', bbox_to_anchor=(1.015, 0.93))
		plt.setp(leg.get_texts(), fontsize='small') 
	
		# Plot titles
		if initstate[0]=='f':
			disp = "\nInitial state: {0}"
			IS_disp = [initstate]
		else:
			IS_disp = ()
			disp = "\nInitial state: $|\psi(0)\\rangle="

			for i in range(len(initstate)):
				IS_disp = IS_disp + ("({2: .3f})|{0},{1}\\rangle".format(*(initstate[i])),)
	
				if i == 0:
					disp = disp + "{" + str(i) + "} "
				else:
					disp = disp + "+ {" + str(i) + "} "	
	
		plt.suptitle("CTQW probability distribution at time $t={}$".format(t))
	
		if (len(list(set(a))) == 1) and (list(set(a))[0] == 0.0):
			plt.title(disp.format(*IS_disp) + "$",
				horizontalalignment='right',multialignment='left', fontsize=11)
		else:
			def_disp = "\nDefects: $"
			for i in range(len(d)):
				def_disp += "{1: .3f}|{0}\\rangle +".format(i+1,d[i],a[i])	
		
			plt.title(disp.format(*IS_disp) + "$" +  def_disp[:-2] + "$",
				horizontalalignment='right',multialignment='left', fontsize=11)

		# save plot
		plt.subplots_adjust(top=0.85)
		pl.savefig(savefile)
	
	# scatter psiX to process 0
	commX = psiX.getComm()
	scatterX, psiX0 = PETSc.Scatter.toZero(psiX)
	scatterX.scatter(psiX, psiX0, False, PETSc.Scatter.Mode.FORWARD)
	
	# scatter psiY to process 0
	commY = psiY.getComm()
	scatterY, psiY0 = PETSc.Scatter.toZero(psiY)
	scatterY.scatter(psiY, psiY0, False, PETSc.Scatter.Mode.FORWARD)
	
	# use process 0 to create the plot
	if rank==0:
		prob_plot_p2(psiX0,psiY0,savefile,t,init_state,d,amp,N)
	
	# deallocate	
	commX.barrier()
	scatterX.destroy()
	psiX0.destroy()
	
	commY.barrier()
	scatterY.destroy()
	psiY0.destroy()
