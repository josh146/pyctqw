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

#~~~~~~~~~~~ determine what modules must be imported ~~~~~~~~~~~~~~~~~~~~~~~~
sparseMat = False
calc_eig = False
LoadFortran = False
LoadSciPy = False
LoadQutip = False

if args.matrix_form[-6:]=='sparse': sparseMat = True

if args.propagator[-9:]=='chebyshev': calc_eig = True

if args.propagator.lower() == 'qutip-expm':	LoadQutip = True

if sparseMat:
	pass
elif args.matrix_form[:4]=='fort'\
    or (args.eig_lib=='lapack' and calc_eig)\
    or args.propagator[:4]=='fort':
	LoadFortran = True

if sparseMat or (args.eig_lib[-5:]=='scipy' and calc_eig)\
    or args.propagator[:6]=='python':
    	LoadSciPy = True

if LoadSciPy:
	try:
		import scipy
		from scipy import sparse
		from scipy.linalg import expm
		scipyVersion = int(scipy.__version__.split('.')[1])
	
		if scipyVersion < 12 and (sparseMat and args.propagator=='python-expm'):
			print "\nERROR: SciPy module is version " + scipy.__version__ \
				+ ", but python-expm for\nsparse matrices requires version >0.12.0"
			print "Please use the python-chebyshev propagator method instead"
			sys.exit(1)
		elif sparseMat and args.propagator=='python-expm':
			from scipy.sparse.linalg import expm
	except SystemExit:
		sys.exit(1)
	except:
		print "\nERROR: SciPy python module cannot be found"
		sys.exit(1)

if LoadQutip:
	try:
		import qutip as qt
	except:
		print "\nERROR: QuTiP python module cannot be found"
		sys.exit(1)

if LoadFortran:
	if args.fortran=="gcc":
		print "Importing libpyctqw_gcc.so..."
		try:	from libpyctqw_gcc import ctqw as fctqw
		except:
			print "\nWARNING: libpyctqw_gcc.so cannot be found. Using libpyctqw_intel.so"
			try:	from libpyctqw_intel import ctqw as fctqw
			except:
				print "\nERROR: libpyctqw_intel.so cannot be found"
				sys.exit()
	else:
		print "Importing libpyctqw_intel.so..."
		try:	from libpyctqw_intel import ctqw as fctqw
		except:
			print "\nWARNING: libpyctqw_intel.so cannot be found. Using libpyctqw_gcc.so"
			try:	from libpyctqw_gcc import ctqw as fctqw
			except:
				print "\nERROR: libpyctqw_gcc.so cannot be found"
				sys.exit()	
else:
	fctqw = func

#~~~~~~~~~~~~~~~~~~~~~~~~~~~ set the variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
		if sparseMat:
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
		if sparseMat:
			psi0 = sparse.lil_matrix((N**2,1),dtype=complex)
			for i in range(len(initialState)):
				psi0[fctqw.coord(*initialState[i][:2]+(N,))-1,0] = initialState[i][-1]
		else:	
			psi0 = np.array([0. for i in range(N**2)],dtype=complex)
			for i in range(len(initialState)):
				psi0[fctqw.coord(*initialState[i][:2]+(N,))-1] = initialState[i][-1]

	else:
		try:	psi0 = np.loadtxt(args.input_state,dtype=complex).reshape(N**2)
		except:	
			print "\nERROR: input state space file " + args.input_state\
				+ " does not exist or is in an incorrect format"
			sys.exit()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Hamiltonian ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
start = time.time()
if sparseMat:
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

elif args.matrix_form[:4]=='fort':
	print '\nCreating the Hamiltonian using Fortran....'

	H1 = fctqw.hamiltonian_1p(d,a,N)

	if args.particles == 1:
		H = H1
	elif args.particles == 2:
		H = fctqw.hamiltonian_2p_noint(H1)

else:
	print '\nCreating the Hamiltonian using NumPy....'

	H1 = np.zeros((N,N))
	np.fill_diagonal(H1,2.0)
	H1.reshape(H1.size)[1:H1.shape[1]*(H1.shape[0]-1):H1.shape[1]+1] = -np.ones((N-1,))
	H1.reshape(H1.size)[H1.shape[1]::H1.shape[1]+1] = -np.ones((N-1,))

	for i in range(len(d)):
		H1[d[i]+N/2-1,d[i]+N/2-1] = 2.0 + a[i]
	
	if args.particles == 1:
		H = H1
	elif args.particles == 2:
		H = np.kron(np.identity(N),H1)+np.kron(H1,np.identity(N))
	
end = time.time()
print '\t\t\t\t\ttime: {: .12f}\n'.format(end-start)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ eigenvalues ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if calc_eig:
	start = time.time()

	if args.eig_lib == 'estimate':
		print 'Estimating the eigenvalues based off the Hamiltonian...'
		if sparseMat:	(Emin,Emax) = (0,1.5*np.sort(H.data)[-1])
		else:		(Emin,Emax) = (0,1.5*np.max(H))

	elif args.eig_lib == 'python-scipy' or sparseMat:		
		if args.eig_solver == 'general':
			(Emin,Emax,start) = func.arpack_extremeev(H,method='gen')
		else:
			(Emin,Emax,start) = func.arpack_extremeev(H,method='SH')
		
	elif args.eig_lib == 'lapack':
		if args.eig_solver == 'general':
			print 'Finding eigenvalues via Lapack....'
			(Emin,Emax) = fctqw.extremeev(H)
		else:
			print 'Finding eigenvalues via Lapack and SH algorithms....'
			(Emin,Emax) = fctqw.sband_extremeev(H)
	
	elif args.eig_lib == 'python-numpy':
		if args.eig_solver == 'general':
			print 'Finding eigenvalues via NumPy....'
			(Emin,Emax) = func.np_extremeev(H,method='SH')
		else:
			print 'Finding eigenvalues via NumPy and SH algorithms....'
			(Emin,Emax) = func.np_extremeev(H,method='gen')
	else:
	
		print '\nWARNING: Unknown linear algebra library'
		print 'Estimating the eigenvalues based off the Hamiltonian...'
		if sparseMat:	(Emin,Emax) = (0,1.5*np.sort(H.data)[-1])
		else:		(Emin,Emax) = (0,1.5*np.max(H))

	print 'Min eigenvalue: {0: .4f}  Max eigenvalue: {1: .4f}'.format(Emin, Emax)
	#print 'Min H element: {0: .4f}  Max H element: {1: .4f}'.format(*np.sort(H.data)[[0,-1]])

	end = time.time()
	print '\t\t\t\t\ttime: {: .12f}\n'.format(end-start)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Quantum Walk ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
start = time.time()
if sparseMat:
	if args.propagator == 'python-expm':
		print 'Calculating exp(-iHt) via SciPy....'
		psi0 = psi0.tocsc()
		psi = sparse.linalg.expm(-1j*H*t)*psi0
		
	elif args.propagator == 'qutip-expm':
		print 'Calculating exp(-iHt) via QuTiP\'s Padre approximation....'
		psi0 = qt.Qobj(psi0.tocsc())
		H = qt.Qobj(H)
		psi = H.expm()*psi0
		psi = psi.full().T[0]
	
	else:
		print 'Calculating exp(-iHt) via the Chebyshev method in Python....'
		psi0 = psi0.tocsc()
		psi = func.qw_cheby(psi0,t,H,Emin,Emax)
	
else:
	if args.propagator == 'python-chebyshev':
		print 'Calculating exp(-iHt) via the Chebyshev method in Python....'
		psi0 = np.matrix(psi0).T
		H = np.matrix(H)
	
		psi = func.qw_cheby(psi0,t,H,Emin,Emax)
		psi = np.array(psi.T)[0]
	
	elif args.propagator == 'python-expm':
		print 'Calculating exp(-iHt) via NumPy....'
		psi0 = np.matrix(psi0).T
		psi = np.matrix(expm(-1j*H*t))*psi0	
		psi = np.array(psi.T)[0]

	elif args.propagator == 'fortran-chebyshev':
		print 'Calculating exp(-iHt) via the Chebyshev method in Fortran....'
		psi = fctqw.qw_cheby(psi0,t,H,Emin,Emax)

	elif args.propagator == 'fortran-burkadt':
		print 'Calculating exp(-iHt) via the Burkadt method in Fortran....'
		psi = fctqw.qw_burkadt(psi0,t,H)

	else:
		print '\nERROR: Unknown CTQW propagator method'
		sys.exit()

end = time.time()
if LoadFortran:	print '\ttime: {: .12f}\n'.format(end-start)
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
	if sparseMat: psi = np.array(psi.T.todense())[0]

	psiX = np.abs(psi)**2.0

	print "Creating marginal probability:\t"+args.output+"/output_psi_t"+str(t)+".txt"
	func.write_psi(psiX, args.output + "/output_psi_t"+str(t)+".txt")

	if args.statespace:
		print "Outputting final state space:\t"+args.output+"/output_statespace_t"+str(t)+".txt"
		func.write_statespace(psi, args.output+"/output_statespace_t"+str(t)+".txt", 1)
	
elif args.particles == 2:
	if sparseMat:	psi = np.array(psi.T.todense())[0]

	if LoadFortran:
		psiX = fctqw.pymarginal(psi,'x',N)
		psiY = fctqw.pymarginal(psi,'y',N)
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
				
				
				
