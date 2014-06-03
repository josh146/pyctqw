#!/usr/bin/env python2.7
from __future__ import print_function
from pprint import pprint

# initialize PETSc
import sys
import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np

# import pyCTQW as qw
import pyCTQW.MPI as qw

# get the MPI rank
rank = PETSc.Comm.Get_rank(PETSc.COMM_WORLD)

# create a graph isomorphism object
gi = qw.GraphISO(verbose=True,freqTol=0.01)

#=========================
#Two Non-Isomorphic Graphs
#=========================
adj1 = np.genfromtxt('../graphs/furer/grid4_furer1.txt')
adj2 = np.genfromtxt('../graphs/furer/grid4_furer2.txt')

# GIcert1 = gi.GIcert(adj1)

# if rank==0:
#     pprint(GIcert1)
#     np.savetxt('out/g4f1cert.txt',GIcert1)

# verify using the inbuilt checking method
if rank == 0:
    print('Testing isomorphism of g4f1 and g4f2:')

cert1, cert2, checkISO = gi.isomorphicQ(adj1, adj2,
                            saveCert1='out/g4f1cert-new.txt',
                            saveCert2='out/g4f2cert-new.txt')

if rank == 0:
    print(checkISO)

#=========================
#Two Isomorphic Graphs
#=========================
adj1p = np.genfromtxt('../graphs/furer/grid4_furer1_perm1.txt')

# verify using the inbuilt checking method
if rank == 0:
    print('Testing isomorphism of g4f1 and g4f1p1:')

cert1, cert1p, checkISO = gi.isomorphicQ(cert1, adj1p,
                            saveCert2='out/g4f1p1cert-new.txt')
if rank == 0:
    print(checkISO)
