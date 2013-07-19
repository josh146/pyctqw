#!/usr/bin/env python2.7
import sys
import os
import numpy as np
import time
import shutil
import errno

import options
import plots
import func

args = options.parse_args()

if (args.expm[:5]=='python' and args.eig_lib[:5]=='python') or args.sparse:
	Fortran = False
	if args.sparse or args.eig_lib=='python-scipy' or arg.expm=='python-chebyshev':
		try:
			from scipy import sparse
			from scipy.sparse.linalg import eigsh, eigs
		except:
			print "\nERROR: libpyctqw_intel.so cannot be found"
			sys.exit()
			
else:
	Fortran = True
	if args.fortran=="gcc":
		print "Importing libpyctqw_gcc.so..."
		try:	from libpyctqw_gcc import ctqw
		except:
			print "\nWARNING: libpyctqw_gcc.so cannot be found. Using libpyctqw_intel.so"
			try:	from libpyctqw_intel import ctqw
			except:
				print "\nERROR: libpyctqw_intel.so cannot be found"
				sys.exit()
	else:
		print "Importing libpyctqw_intel.so..."
		try:	from libpyctqw_intel import ctqw
		except:
			print "\nWARNING: libpyctqw_intel.so cannot be found. Using libpyctqw_gcc.so"
			try:	from libpyctqw_gcc import ctqw
			except:
				print "\nERROR: libpyctqw_gcc.so cannot be found"
				sys.exit()

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
else:
	print '\nERROR: only 1 or 2 particle quantum walks supported'
	sys.exit()
	
#~~~~~~~~~~~~~~~~~~~~~~~~~ Print out simulation info ~~~~~~~~~~~~~~~~~~~~~~~~
print "\nPerforming a {} particle CTQW:\n".format(args.particles)
print "\t N={0} \t t={1}".format(N,t)

IS_disp = ()
disp = "\t Initial state: "

for i in range(len(initialState)):
	if args.particles == 1:
		IS_disp = IS_disp + ("{1}|{0}>".format(*(initialState[i])),)
	else:
		IS_disp = IS_disp + ("{2}|{0},{1}>".format(*(initialState[i])),)
	
	if i == 0:
		disp = disp + "{" + str(i) + "} "
	else:
		disp = disp + "+ {" + str(i) + "} "
		
if args.input_state=="":
	print disp.format(*IS_disp)
else:
	print "\t Loading initial state from file " + args.input_state

for i in range(len(d)):
	print "\t Defect {0}: node {1}, amplitude {2}".format(i+1,d[i],a[i])
	
#~~~~~~~~~~~~~~~~~~~~~~~ Create the initial statespace ~~~~~~~~~~~~~~~~~~~~~~~
if args.particles == 1:
	
	if args.input_state=="":
		if args.sparse:
			psi0 = sparse.lil_matrix((N,1),dtype=complex)
			for i in range(len(initialState)):
				psi0[int(initialState[i][0])+N/2-1,0] = initialState[i][1]
		else:
			psi0 = np.array([0. for i in range(N)],dtype=complex)
			for i in range(len(initialState)):
				psi0[int(initialState[i][0])+N/2-1] = initialState[i][1]
	else:
		try:	psi0 = np.loadtxt(args.input_state,dtype=complex).reshape(N)
		except:
			print "\nERROR: input state space file " + args.input_state\
				+ " does not exist or is in an incorrect format"
			sys.exit()

elif args.particles == 2:
	if args.input_state=="":
		if args.sparse:
			psi0 = sparse.lil_matrix((N**2,1),dtype=complex)
			for i in range(len(initialState)):
				if Fortran:	psi0[ctqw.coord(*initialState[i][:2]+(N,))-1,0] = initialState[i][-1]
				else:		psi0[func.coord(*initialState[i][:2]+(N,))-1,0] = initialState[i][-1]
		else:	
			psi0 = np.array([0. for i in range(N**2)],dtype=complex)
			for i in range(len(initialState)):
				if Fortran:	psi0[ctqw.coord(*initialState[i][:2]+(N,))-1] = initialState[i][-1]
				else:		psi0[func.coord(*initialState[i][:2]+(N,))-1] = initialState[i][-1]

			
	else:
		try:	psi0 = np.loadtxt(args.input_state,dtype=complex).reshape(N**2)
		except:	
			print "\nERROR: input state space file " + args.input_state\
				+ " does not exist or is in an incorrect format"
			sys.exit()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Hamiltonian ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
start = time.time()
if args.sparse:
	print '\nCreating the Hamiltonian sparsely using SciPy....'
	
	H1 = sparse.lil_matrix((N, N))
	
	H1.setdiag([2.0 for x in range(N)], k=0)
	H1.setdiag([-1.0 for x in range(N)], k=1)
	H1.setdiag([-1.0 for x in range(N)], k=-1)
	
	for i in range(len(d)):
		H1[d[i]+N/2-1,d[i]+N/2-1] = 2.0 + a[i]
		
	if args.particles == 1:
		H = H1.tocsc()
	elif args.particles == 2:
		H = sparse.kronsum(H1, H1, format='csc')
	
else:
	print '\nCreating the Hamiltonian using Fortran....'

	H1 = ctqw.hamiltonian_1p(d,a,N)

	if args.particles == 1:
		H = H1
	elif args.particles == 2:
		H = ctqw.hamiltonian_2p_noint(H1)
		
end = time.time()
print '\t\t\t\t\ttime: {: .12f}\n'.format(end-start)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ eigenvalues ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if args.expm!='burkadt':
	start = time.time()
	
	if args.eig_lib == 'python-scipy' or args.sparse:
		if args.sparse:	H_sparse = H
		else:		H_sparse = sparse.csc_matrix(H)		
		
		Emax = 'none'
		Emin = 'none'
		count = 1
		
		if args.eig_general:
			sys.stdout.write('Finding eigenvalues via SciPy\'s ARPACK bindings....'+'\b')
			while (Emin=='none' or Emax=='none') and count<11:
				sys.stdout.write('Attempt {0:02d}'.format(count)+"\b"*10)
				sys.stdout.flush()
				try:
					start = time.time()
					Emax = eigs(H_sparse,count,which='LA')[0][-1].real
					Emin = eigs(H_sparse,count,which='LA',sigma=0)[0][-1].real
					break
				except:
					count += 1
			else:
				print "\nERROR: SciPy band storage eigenvalue solver could not converge"
				sys.exit()
		else:
			sys.stdout.write('Finding eigenvalues via SciPy\'s ARPACK bindings and band storage....'+'\b')
			while (Emin=='none' or Emax=='none') and count<11:
				sys.stdout.write('Attempt {0:02d}'.format(count)+"\b"*10)
				sys.stdout.flush()
				try:
					start = time.time()
					Emax = eigsh(H_sparse,count,which='LA')[0][-1].real
					Emin = eigsh(H_sparse,count,which='LA',sigma=0)[0][-1].real
					sys.stdout.write("\n")
					break
				except:
					count += 1
			else:
				print "\nERROR: SciPy band storage eigenvalue solver could not converge"
				sys.exit()

	elif args.eig_lib == 'lapack':
		if args.eig_general:
			print 'Finding eigenvalues via Lapack....'
			(Emin,Emax) = ctqw.extremeev(H)
		else:
			print 'Finding eigenvalues via Lapack and band storage....'
			(Emin,Emax) = ctqw.sband_extremeev(H)
		
	elif args.eig_lib == 'python-numpy':
		if args.eig_general:
			print 'Finding eigenvalues via NumPy....'
			evals = np.linalg.eigvals(H)
		else:
			print 'Finding eigenvalues via NumPy and band storage....'
			evals = np.linalg.eigvalsh(H)
		
		Emax = max(evals.real)
		Emin = min(evals.real)
	
	else:
		print '\nERROR: Unknown linear algebra library'
		sys.exit()

	print 'Min eigenvalue: {0: .4f}  Max eigenvalue: {1: .4f}'.format(Emin, Emax)

	end = time.time()
	print '\t\t\t\t\ttime: {: .12f}\n'.format(end-start)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Quantum Walk ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
start = time.time()

if args.expm == 'python-chebyshev' or args.sparse:
	print 'Calculating exp(-iHt) via the Chebyshev method in Python....'
	if args.sparse:
		psi0 = psi0.tocsc()	
	else:
		psi0 = np.matrix(psi0).T
		H = np.matrix(H)
		
	psi = func.qw_cheby(psi0,t,H,Emin,Emax)
	
	if not args.sparse:
		psi = np.array(psi.T)[0]
	
elif args.expm == 'chebyshev':
	print 'Calculating exp(-iHt) via the Chebyshev method....'
	psi = ctqw.qw_cheby(psi0,t,H,Emin,Emax)
	
elif args.expm == 'burkadt':
	print 'Calculating exp(-iHt) via the Burkadt method....'
	psi = ctqw.qw_burkadt(psi0,t,H)
	
else:
	print '\nERROR: Unknown matrix exponential method'
	sys.exit()

end = time.time()
if Fortran:	print '\ttime: {: .12f}\n'.format(end-start)
else:		print '\t\t\t\t\ttime: {: .12f}\n'.format(end-start)

if os.path.exists('coeff.txt'):
	coeff = np.loadtxt('coeff.txt')
	os.remove('coeff.txt')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
# create output directory if it doesn't exist
try:
	os.mkdir(args.output)
except OSError as exception:
	if exception.errno != errno.EEXIST:
		raise

# marginal probabilities and statespace output
if args.particles == 1:
	if args.sparse: psi = np.array(psi.T.todense())[0]
	
	psiX = np.abs(psi)**2.0
	
	print "Creating marginal probability:\t"+args.output+"/output_psi_t"+str(t)+".txt"
	func.write_psi(psiX, args.output + "/output_psi_t"+str(t)+".txt")
	
	if args.statespace:
		print "Outputting final state space:\t"+args.output+"/output_statespace_t"+str(t)+".txt"
		func.write_statespace(psi, args.output+"/output_statespace_t"+str(t)+".txt", 1)
		
elif args.particles == 2:
	if args.sparse:	psi = np.array(psi.T.todense())[0]
	
	if Fortran:
		psiX = ctqw.pymarginal(psi,'x',N)
		psiY = ctqw.pymarginal(psi,'y',N)
	else:	
		psiX = func.pymarginal(psi,'x',N)
		psiY = func.pymarginal(psi,'y',N)
	
	print "Creating marginal probability for particle 1:\t"\
		+args.output+"/output_psiX_t"+str(t)+".txt"
	func.write_psi(psiX, args.output + "/" + "output_psiX_t"+str(t)+".txt")
	
	print "Creating marginal probability for particle 2:\t"\
		+args.output+"/output_psiX_t"+str(t)+".txt"
	func.write_psi(psiY, args.output + "/" + "output_psiY_t"+str(t)+".txt")
	
	if args.statespace:
		print "Outputting final state space:\t\t\t"+args.output+"/output_statespace_t"+str(t)+".txt"
		func.write_statespace(psi, args.output+"/output_statespace_t"+str(t)+".txt", 2)
				
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if args.particles == 1:
	print "Creating probability plot:\t"+args.output+"/output_prob_t"+str(t)+".png"
	plots.prob_plot_p1(args.output+"/output_psi_t"+str(t)+".txt",
			args.output+"/output_prob_t"+str(t)+".png",
			t,initialState,d,a)
else:
	print "Creating marginal probability plot:\t\t"+args.output+"/output_prob_t"+str(t)+".png"
	plots.prob_plot_p2(args.output+"/output_psiX_t"+str(t)+".txt",
			args.output+"/output_psiY_t"+str(t)+".txt",
			args.output+"/output_prob_t"+str(t)+".png",
			t,initialState,d,a)













