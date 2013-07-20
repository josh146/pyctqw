#!/usr/bin/env python2.7
import sys, time
import numpy as np

class ref:
	def __init__(self, obj):
    		self.obj = obj
	def get(self):
    		return self.obj
	def set(self, obj):
    		self.obj = obj

#~~~~~~~~~~~~~~~~~~~ ctqw functions ~~~~~~~~~~~~~~~~~~~
def qw_cheby(psi,dt,H,Emin,Emax):
	from scipy.special import j0, j1, jn
	
	alpha = (Emax-Emin)*dt/2.0
	
	phi0 = psi
	phi1 = -(2.0*H*psi - (Emax+Emin)*psi)/(Emax-Emin)
	U = j0(alpha)*phi0 + 2.0j*j1(alpha)*phi1
	
	terms = 0
	while abs(2.0*jn(terms,alpha)) > 1e-18:
		terms += 1
		
	for m in range(2,terms+1):
		#ctqw.progressbar(m,terms+1-2)
		#progressbar(m,terms+1-2)
		
		phi2 = -2.0*(2.0*H*phi1 - (Emax+Emin)*phi1)/(Emax-Emin) - phi0
		U = U + 2.0*(1j**m)*jn(m,alpha)*phi2
		
		phi0 = phi1
		phi1 = phi2
	
	return np.exp(-1j*(Emax+Emin)*dt/2.0)*U

def pymarginal(psi,p,N):
	marginal = np.ndarray((N))
	
	if p=='x':
		for i in range(N):
			marginal[i] = np.sum(np.abs( psi[i*N:(i+1)*N-1] )**2.0)
	elif p=='y':
		for i in range(N):
			marginal[i] = np.sum(np.abs( [psi[i+j*N] for j in range(N)] )**2.0)
			
	return marginal
	
def coord(x,y,N):
	return int(N*(x+N/2-1)+y+N/2)
    		
#~~~~~~~~~~~~~~~~~~~ output functions ~~~~~~~~~~~~~~~~~
def write_psi(psi, filename):
	N = len(psi)
	with open(filename,'w') as f:
		for i in range(N):
			f.write('{0}\t{1: .12e}\n'.format(i-np.floor(N/2)+1,psi[i].real))
			
def write_statespace(psi, filename, num_p):
	N = len(psi)
	
	if num_p==1:
		with open(filename,'w') as f:
			for i in range(N):
				f.write('{0: .12e}\n'.format(psi[i]))
	elif num_p==2:
		nn = int(np.sqrt(N))
		state_space = psi.reshape((nn,nn))		
		with open(filename,'w') as f:
			for i in range(nn):
				ss_disp = ""
				for j in range(nn):
					ss_disp = ss_disp + '{0: .12e}\t'.format(state_space[i,j])
				f.write(ss_disp+'\n')
