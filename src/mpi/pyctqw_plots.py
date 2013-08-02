#!/usr/bin/python
from petsc4py import PETSc
from matplotlib import pyplot as plt
import numpy as np
import pylab as pl
import matplotlib

def plot_marginal(psiX,psiY,savefile,t,init_state,d,amp,N,rank):
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
