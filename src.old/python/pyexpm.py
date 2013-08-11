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
from scipy import sparse

args = options.parse_args()

# import matrix
if args.input_matrix == "":
	if type(args.matrix[0]) is str:
			
		if "".join(args.matrix) == 'random_real':
			N = args.leading_dim
			print "Creating a random {0}x{0} matrix...".format(N)
			A = np.random.random((N,N))
			if args.eig_solver == 'SH':
				print "Creating a symmetric matrix..."
				A = A + A.T.conjugate()
			sparseMat = False
		
			print "Outputting A:\t\t\t"+args.output+"/A_input.txt"
			func.write_matrix(A,args.output + "/A_input.txt")
		
		elif "".join(args.matrix)[:12] == 'hamiltonian':
			N = args.leading_dim
			A = sparse.lil_matrix((N, N))
			print "Creating a sparse {0}x{0} CTQW Hamiltonian matrix...".format(N)
			A.setdiag([2.0 for x in range(N)], k=0)
			A.setdiag([-1.0 for x in range(N)], k=1)
			A.setdiag([-1.0 for x in range(N)], k=-1)
			
			if args.matrix[-2:] == "2D":
				print "Creating a sparse {0}x{0} 2D CTQW Hamiltonian matrix...".format(N**2)
				A = sparse.kronsum(A, A, format='csc')
		
			if not args.sparse_alg:
				print "Converting the Hamiltonian to a dense matrix"
				A = A.toarray()
				sparseMat = False
			else:
				A = A.tocsc()
		
			sparseMat = True
		
		
		else:
			N = args.leading_dim
			print "Creating a complex random {0}x{0} matrix...".format(N)
			A = np.random.random((N,N)) + 1j*np.random.random((N,N))
			if args.eig_solver == 'SH':
				print "Creating a Hermitian matrix..."
				A = A + A.T.conjugate()
			sparseMat = False
		
			print "Outputting A:\t\t\t"+args.output+"/A_input.txt"
			func.write_matrix(A,args.output + "/A_input.txt")
		
	elif type(args.matrix[0]) is tuple:
		N = args.leading_dim
		A = sparse.lil_matrix((N,N),dtype=complex)
		for i in range(len(args.matrix)):
			A[args.matrix[i][0],args.matrix[i][1]] = args.matrix[i][2]
		A = A.tocsc()
		sparseMat = True
		
		print "Imported a sparse matrix"
		
	else:	
		A = np.array(args.matrix)
		N = A.shape[0]
		sparseMat = False
		
		print "Imported a dense matrix"

else:
	try:
		A = np.loadtxt(args.input_matrix,dtype=complex)
		if args.sparse_alg:
			A = sparse.csc_matrix(A)
	except:
		print "\nERROR: input state space file " + args.input_matrix\
				+ " does not exist or is in an incorrect format"
		sys.exit()

#~~~~~~~~~~~ determine what modules must be imported ~~~~~~~~~~~~~~~~~~~~~~~~
if args.eig_lib=='lapack':
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

# calc the eigenvalues
start = time.time()	
if args.eig_lib == 'python-scipy' or sparseMat:		
	if args.eig_solver == 'general':
		try:	
			print '\nFinding eigenvalues via SciPy....'
			(Emin,Emax,start) = func.arpack_extremeev(A,method='gen')
		except:
			print '\nFinding eigenvalues via NumPy....'
			(Emin,Emax) = func.np_extremeev(A.toarray(),method='gen')			
	else:
		try:	
			print '\nFinding eigenvalues via SciPy and SH algorithms....'
			(Emin,Emax,start) = func.arpack_extremeev(A,method='SH')
		except:
			print '\nFinding eigenvalues via NumPy and SH algorithms....'
			(Emin,Emax) = func.np_extremeev(A.toarray(),method='SH')	
	
elif args.eig_lib == 'lapack':
	if args.eig_solver == 'general':
		print '\nFinding eigenvalues via Lapack....'
		(Emin,Emax) = fctqw.extremeev(A)
	else:
		print '\nFinding eigenvalues via Lapack and SH algorithms...'
		(Emin,Emax) = fctqw.sband_extremeev(A)

elif args.eig_lib == 'python-numpy':
	if args.eig_solver == 'general':
		print '\nFinding eigenvalues via NumPy....'
		(Emin,Emax) = func.np_extremeev(A,method='gen')
	else:
		print '\nFinding eigenvalues via NumPy and SH algorithms....'
		(Emin,Emax) = func.np_extremeev(A,method='SH')

else:
	print '\nERROR: Unknown linear algebra library'
	sys.exit()

print 'Min eigenvalue: {0: .4f}  Max eigenvalue: {1: .4f}'.format(Emin, Emax)

end = time.time()
print '\t\t\t\t\ttime: {: .12f}\n'.format(end-start)

# calculate exp(A)
start = time.time()
print 'Finding exp(A)....'
expA = func.expm_cheby(A,Emin,Emax,sparseMode=args.sparse_alg)
end = time.time()
print '\t\t\t\t\ttime: {: .12f}\n'.format(end-start)

if N < 4:
	print expA

# write matrix to a file
print "Outputting exp(A):\t\t\t"+args.output+"/expA_output.txt"
func.write_matrix(expA,args.output + "/expA_output.txt")





