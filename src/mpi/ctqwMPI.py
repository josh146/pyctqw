#!/usr/bin/python
import sys, os, errno, petsc4py
petsc4py.init(sys.argv)

from petsc4py import PETSc
from libctqwMPI import ctqwmpi

from matplotlib import pyplot as plt
import numpy as np
import pylab as pl
import matplotlib

def plot_marginal(psiX,psiY,savefile,t,init_state,d,amp,N):
	commX = psiX.getComm()
	scatterX, psiX0 = PETSc.Scatter.toZero(psiX)
	scatterX.scatter(psiX, psiX0, False, PETSc.Scatter.Mode.FORWARD)
	
	commY = psiY.getComm()
	scatterY, psiY0 = PETSc.Scatter.toZero(psiY)
	scatterY.scatter(psiY, psiY0, False, PETSc.Scatter.Mode.FORWARD)
	
	if rank==0:
		prob_plot_p2(psiX0,psiY0,savefile,t,init_state,d,amp,N)
		
	commX.barrier()
	scatterX.destroy()
	psiX0.destroy()
	
	commY.barrier()
	scatterY.destroy()
	psiY0.destroy()

def prob_plot_p2(psiX,psiY,savefile,t,initstate,d,a,N):
	
	probX = np.real(np.asarray(psiX))
	probY = np.real(np.asarray(psiY))
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
	IS_disp = ()
	disp = "\nInitial state: $|\psi(0)\\rangle="

	for i in range(len(initstate)):
		IS_disp = IS_disp + ("({2: .3f})|{0},{1}\\rangle".format(*(initstate[i])),)
	
		if i == 0:
			disp = disp + "{" + str(i) + "} "
		else:
			disp = disp + "+ {" + str(i) + "} "
	
	def_disp = "\nDefects: $"
	for i in range(len(d)):
		def_disp += "{1: .3f}|{0}\\rangle +".format(i+1,d[i],a[i])	
	
	plt.suptitle("CTQW probability distribution at time $t={}$".format(t))
		
	plt.title(disp.format(*IS_disp) + "$" + def_disp[:-2] + "$",
		horizontalalignment='right',multialignment='left', fontsize=11)

	# save plot
	plt.subplots_adjust(top=0.85)
	pl.savefig(savefile)

OptDB = PETSc.Options()
N = OptDB.getInt('N', 10)
t = OptDB.getInt('t', 1)
d = [3,4]
amp = [2.0,1.5]
init_state = [[0.,1.,1.0/np.sqrt(2.0)], [1.,1.,1.0/np.sqrt(2.0)]]

rank = PETSc.Comm.Get_rank(PETSc.COMM_WORLD)

H = PETSc.Mat()
H.create(PETSc.COMM_WORLD)

psi0 = PETSc.Vec()
psi0.create(PETSc.COMM_WORLD)
psi0.setSizes(N**2)
psi0.setUp()

ctqwmpi.p2_init(psi0.fortran,init_state,N)
ctqwmpi.hamiltonian_2p_line(H.fortran,d,amp,N)

Emax,Emax_error,ierr = ctqwmpi.min_max_eigs(H.fortran,rank,'max','null','null',35,0.0,0,False)
Emin,Emin_error,ierr = ctqwmpi.min_max_eigs(H.fortran,rank,'min','null','null',35,0.0,0,False)

psi = psi0.duplicate()
ctqwmpi.qw_cheby(psi0.fortran,psi.fortran,t,H.fortran,Emin,Emax,rank,N)

psiX = PETSc.Vec()
psiX.create(PETSc.COMM_WORLD)
psiX.setSizes(N)
psiX.setUp()
psiY = psiX.duplicate()
ctqwmpi.marginal(psi.fortran,psiX.fortran,'x',N)
ctqwmpi.marginal(psi.fortran,psiY.fortran,'y',N)

# create output directory if it doesn't exist
try:
	os.mkdir('./out')
except OSError as exception:
	if exception.errno != errno.EEXIST:
		raise
		
plot_marginal(psiX,psiY,'out/plot.png',t,init_state,d,amp,N)

H.destroy()
psi.destroy()
psi0.destroy()
psiX.destroy()
psiY.destroy()
