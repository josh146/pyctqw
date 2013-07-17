#!/usr/bin/env python2.7
import numpy as np
from libpyctqw import ctqw
import time

# set the variables
N = 50
t = 1.0;

# defect locations
a1 = 1.0;	d1 = 3
a2 = 1.0;	d2 = 4

# initial state
psi0 = np.array([0. for i in range(N**2)])
psi0[ctqw.coord(0,0,N)-1] = 1/np.sqrt(2)
psi0[ctqw.coord(1,1,N)-1] = 1/np.sqrt(2)

# Hamiltonian
start = time.time()
H1 = ctqw.hamiltonian_1p(d1,a1,d2,a2,N)
H2 = ctqw.hamiltonian_2p_noint(H1)
end = time.time()
print 'Hamiltonian creation time: {}'.format(end-start)

start = time.time()
evals = np.linalg.eigvalsh(H2)
Emax = max(evals.real)
Emin = min(evals.real)
end = time.time()
print 'Eigenvalue calc time: {}'.format(end-start)

# quantum walk
start = time.time()
#psi = ctqw.qw_burkadt(psi0,t,H2)
psi = ctqw.qw_cheby(psi0,t,H2,Emin,Emax,False)
end = time.time()
print 'QW propagation time: {}'.format(end-start)

# marginal probabilities
psiX = ctqw.pymarginal(psi,'x',N)
for i in range(len(psiX)):
	print '{0}\t{1: .12e}'.format(i-np.floor(N/2)+1,psiX[i].real)

