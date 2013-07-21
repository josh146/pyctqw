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
    		
#~~~~~~~~~~~~~~~~~~~ eigenvalue solvers ~~~~~~~~~~~~~~~
def arpack_extremeev(H,method='gen',numpy_fallback=False):
	from scipy.sparse.linalg import eigs, eigsh
	
	if method == 'SH':
		eigFunc = eigsh
	elif method == 'gen':
		eigFunc = eigs
	else:
		print "\nERROR: unknown eigenvalue solver method"
		sys.exit(1)

	Emax = 'none'
	Emin = 'none'
	count = 1

	if method == 'SH':	sys.stdout.write('Finding eigenvalues via SciPy\'s ARPACK bindings and band storage....'+'\b')	
	elif method == 'gen':	sys.stdout.write('Finding eigenvalues via SciPy\'s ARPACK bindings....'+'\b')
	
	while (Emin=='none' or Emax=='none') and count<11:
		sys.stdout.write('Attempt {0:02d}'.format(count)+"\b"*10)
		sys.stdout.flush()
		try:
			start = time.time()
			Emax = eigFunc(H,count,which='LA')[0][-1].real
			Emin = eigFunc(H,count,which='LA',sigma=0)[0][-1].real
			sys.stdout.write("\n")
			break
		except:
			count += 1
	else:
		print "\nERROR: SciPy eigenvalue solver could not converge"
		if numpy_fallback:
			print "WARNING: NumPy eigenvalue solver will be used as a fallback; this may hang for very large matrices"
			try:
				(Emin,Emax) = np_extremeev(H,method=method)
			except:
				print "ERROR: NumPy fallback will not work with sparse matrices"
				sys.exit(1)
		else:
			sys.exit()
	
	return (Emin,Emax,start)


def np_extremeev(H,method='gen'):
	if method == 'SH':
		eigFunc = np.linalg.eigvalsh
	elif method == 'gen':
		eigFunc = np.linalg.eigvals
	else:
		print "\nERROR: unknown eigenvalue solver method"
		sys.exit(1)

	return np.sort(np.linalg.eigvals(H))[[0,-1]].tolist()

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
	
#~~~~~~~~~~~~ chebyshev matrix exponential ~~~~~~~~~~~~
def expm_cheby(H,Emin,Emax,sparseMode=False):
	from scipy.special import jn
	from scipy import sparse

	if sparseMode or sparse.issparse(H):
		ident = sparse.identity(H.shape[0]).tocsc()
		if sparse.issparse(H):
			inputmatrix = H.format
			if not sparse.isspmatrix_csc(H):
				H = H.tocsc()
		else:
			inputmatrix = 'dense'
			H = sparse.csc_matrix(H)
	else:
		inputmatrix = ''
		ident = np.matrix(np.identity(H.shape[0]))
		H = np.matrix(H)
	
	alpha = (Emax-Emin)/2.0
	
	phi0 = ident
	phi1 = (2.0*H - (Emax+Emin)*phi0)/(Emax-Emin)
	U = jn(0,alpha)*phi0 + 2.0*jn(1,alpha)*phi1
	
	terms = 0
	while abs(2.0*jn(terms,alpha)) > 1e-18:
		terms += 1
		
	for m in range(2,terms+1):
		phi2 = 2.0*(2.0*H*phi1 - (Emax+Emin)*phi1)/(Emax-Emin) + phi0
		U = U + 2.0*jn(m,alpha)*phi2
		
		phi0 = phi1
		phi1 = phi2
		
	expm = np.exp((Emax+Emin)/2.0)*U
		
	if sparseMode or inputmatrix:
		if inputmatrix == 'dense':
			return expm.toarray()
		else:
			return getattr(expm,'to'+inputmatrix)()
	else:
		return np.array(expm)

    		
#~~~~~~~~~~~~~~~~~~~ output functions ~~~~~~~~~~~~~~~~~
def write_matrix(A, filename):
	N = A.shape[0]
	with open(filename,'w') as f:
		for i in range(N):
			row_tuple = ()
			row_disp = ""
			for j in range(N):
				row_tuple += ("{0: .12e}".format(A[i,j]),)
				if j == N-1:	row_disp += "{" + str(j) + "}\n"
				else:		row_disp += "{" + str(j) + "} "
			f.write(row_disp.format(*row_tuple))

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
