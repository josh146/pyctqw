#!/usr/bin/env python2.7
import sys
import os
import numpy as np
from libpyctqw import ctqw
import time
import options

args = options.parse_args()

# set the variables
N = args.grid_length
t = args.time

# defect locations
if len(args.defect_nodes) == len(args.defect_amp):
	d=args.defect_nodes
	a=args.defect_amp
else:
	print "ERROR: each specified defect location MUST have \
		a respective defect amplitude"
	sys.exit()

if args.particles == 1:
	# 1P initial state
	initialState = args.p1_initial_state
elif args.particles == 2:
	# 2P initial state
	initialState = args.p2_initial_state
	
#~~~~~~~~~~~~~~~~~~~~~~~ Create the initial statespace ~~~~~~~~~~~~~~~~~~~~~~~
if args.particles == 1:
	psi0 = [0. for i in range(N)]
	
	for i in range(len(initialState)):
		psi0[int(initialState[i][0])+N/2-1] = initialState[i][1]

elif args.particles == 2:
	psi0 = [0. for i in range(N**2)]
	
	for i in range(len(initialState)):
		psi0[ctqw.coord(*initialState[i][:2]+(N,))-1] = initialState[i][-1]
else:
	print 'ERROR: only 1 or 2 particle quantum walks supported'
	sys.exit()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Hamiltonian ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print 'Creating the Hamiltonian....'

start = time.time()

H1 = ctqw.hamiltonian_1p(d,a,N)

if args.particles == 1:
	H = H1
elif args.particles == 2:
	H = ctqw.hamiltonian_2p_noint(H1)

end = time.time()
print '\t\t\t\t\ttime: {: .12f}\n'.format(end-start)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ eigenvalues ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
start = time.time()

if args.evallib == 'lapack':
	if args.evalgen:
		print 'Finding eigenvalues via Lapack....'
		(Emin,Emax) = ctqw.extremeev(H)
	else:
		print 'Finding eigenvalues via Lapack and band storage....'
		(Emin,Emax) = ctqw.sband_extremeev(H)
elif args.evallib == 'numpy':
	if args.evalgen:
		print 'Finding eigenvalues via NumPy....'
		evals = np.linalg.eigvals(H)
	else:
		print 'Finding eigenvalues via NumPy and band storage....'
		evals = np.linalg.eigvalsh(H)
		
	Emax = max(evals.real)
	Emin = min(evals.real)
else:
	print 'ERROR: Unknown linear algebra library'
	sys.exit()

end = time.time()
print '\t\t\t\t\ttime: {: .12f}\n'.format(end-start)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Quantum Walk ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
start = time.time()

if args.expm == 'chebyshev':
	print 'Calculating exp(-iHt) via the Chebyshev method....'
	psi = ctqw.qw_cheby(psi0,t,H,Emin,Emax,False)
elif args.expm == 'burkadt':
	print 'Calculating exp(-iHt) via the Burkadt method....'
	psi = ctqw.qw_burkadt(psi0,t,H)
else:
	print 'ERROR: Unknown matrix exponential method'
	sys.exit()

end = time.time()
print '\ttime: {: .12f}\n'.format(end-start)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Post processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# marginal probabilities
if args.particles == 1:
	psiX = np.abs(psi)**2
	for i in range(len(psiX)):
		print '{0}\t{1: .12e}'.format(i-np.floor(N/2)+1,psiX[i].real)
		
elif args.particles == 2:
	psiX = ctqw.pymarginal(psi,'x',N)
	psiY = ctqw.pymarginal(psi,'y',N)
	for i in range(len(psiX)):
		print '{0}\t{1: .12e}'.format(i-np.floor(N/2)+1,psiX[i].real)















