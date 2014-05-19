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
gi = qw.GraphISO(verbose=True)

#=========================
#Two Non-Isomorphic Graphs
#=========================
adj1 = np.genfromtxt('../graphs/furer/grid3_furer1.txt')
adj2 = np.genfromtxt('../graphs/furer/grid3_furer2.txt')

GIcert1 = gi.GIcert(adj1)

print(rank, GIcert1)

# # verify using the inbuilt checking method
# if rank == 0:S
#     print('Testing isomorphism of g3f1 and g3f2:')

# cert1, cert2, checkISO = gi.isomorphicQ(adj1, adj2,
#                             saveCert1='out/g3f1cert.txt',
#                             saveCert2='out/g3f2cert.txt')

# if rank == 0:
#     print(checkISO)

# #=========================
# #Two Isomorphic Graphs
# #=========================
# adj1p = np.genfromtxt('../graphs/furer/grid3_furer1_perm1.txt')

# # verify using the inbuilt checking method
# if rank == 0:
#     print('Testing isomorphism of g3f1 and g3f1p1:')

# cert1, cert1p, checkISO = gi.isomorphicQ(cert1, adj1p,
#                             saveCert2='out/g3f1p1cert.txt')
# if rank == 0:
#     print(checkISO)
