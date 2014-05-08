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
gi = qw.GraphISO()

#=========================
#Two Non-Isomorphic Graphs
#=========================
adj1 = np.genfromtxt('../graphs/k-equiv/K8.txt')
adj2 = np.genfromtxt('../graphs/k-equiv/K8_twisted01.txt')

# generate certificates for these two graphs
# cert1 = gi.GIcert(adj1)
# cert2 = gi.GIcert(adj2)

# print the certificates
# if rank == 0:
#     print('The GI certificate of K1.txt:')
#     pprint(cert1)

#     print('\nThe GI certificate of K1-twisted.txt:')
#     pprint(cert2)

# verify using the inbuilt checking method
if rank == 0:
    print('\nTesting isomorphism of K8.txt and K8_twisted01.txt:')

checkISO = gi.isomorphicQ(adj1, adj2)
if rank == 0:
    print(checkISO)

#=========================
#    Isomorphic Graphs
#=========================
# if rank == 0:
#     print('\nTesting isomorphism of K1.txt and K1-perm1.txt:')

# adj1 = np.genfromtxt('../graphs/k-equiv/K1.txt')
# adj2 = np.genfromtxt('../graphs/k-equiv/K1-perm1.txt')

# checkISO = gi.isomorphicQ(adj1, adj2)
# if rank == 0:
#     print(checkISO)

# if rank == 0:
#     print('\nTesting isomorphism of K1-twisted.txt and K1-twisted-perm1.txt:')

# adj1 = np.genfromtxt('../graphs/k-equiv/K1-twisted.txt')
# adj2 = np.genfromtxt('../graphs/k-equiv/K1-twisted-perm1.txt')

# checkISO = gi.isomorphicQ(adj1, adj2)
# if rank == 0:
#     print(checkISO)
