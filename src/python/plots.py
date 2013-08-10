#!/usr/bin/env python2.7
import argparse
import os
from matplotlib import pyplot as plt
import numpy as np
import pylab as pl
import matplotlib

def plotCoeff(coeff,savefile,alpha):	
	n = np.arange(len(coeff))

	# create plot
	fig = plt.figure()
	plt.plot(n, coeff,'r-',marker='.')
	plt.ylabel("$|J_n(\\alpha)|$",rotation=90)
	plt.grid(b=None, which='major', axis='both', linestyle='-', alpha=0.3)
	plt.xlabel("$n$")
	plt.yscale('log')
		
	plt.title('$\\alpha={}$'.format(alpha))

	# save plot
	pl.savefig(savefile)

def prob_plot_p1(data,savefile,t,initstate,d,a):
	# import data
	x,prob = np.genfromtxt(data, unpack=True)
	
	# create plot
	fig = plt.figure()
	plt.plot(x, prob)
	plt.ylabel("$|\langle j|\psi\\rangle|^2$",rotation=90)
	plt.grid(b=None, which='major', axis='both', linestyle='-', alpha=0.3)
	plt.xlabel("$j$")
	
	# Plot titles
	IS_disp = ()
	disp = "\nInitial state: $|\psi(0)\\rangle="

	for i in range(len(initstate)):
		IS_disp = IS_disp + ("({1: .3f})|{0}\\rangle".format(*(initstate[i])),)
	
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

def prob_plot_p2(dataX,dataY,savefile,t,initstate,d,a):
	# import data
	x,probX = np.genfromtxt(dataX, unpack=True)
	y,probY = np.genfromtxt(dataY, unpack=True)
	
	# create plot
	lbl = [0,1];
	
	fig = plt.figure()
	lbl[0], = plt.plot(x, probX, 'b-')
	lbl[1], = plt.plot(y, probY, 'r--')
	
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
