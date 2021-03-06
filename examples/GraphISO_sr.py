#!/usr/bin/env python2.7
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
#  Two Isomorphic Graphs
#=========================
# import two adjacency matrices; these are permutations
# of the first graph in the SR-25-12-5-6, and thus are isomorphic
adj1 = np.genfromtxt('../graphs/strong-regular-25-12-5-6/1-permutations/SR1_perm_1.txt')
adj1p = np.genfromtxt('../graphs/strong-regular-25-12-5-6/1-permutations/SR1_perm_3.txt')

# generate certificates for these two graphs
cert1 = gi.GIcert(adj1)
cert1p = gi.GIcert(adj1p)

# print the certificates side by side;
# they should be identical
if rank == 0:
    print 'The GI certificates of two permutations of SR(25,12,5,6) graph 1:'
    for i in range(len(cert1)):
        print cert1[i], cert1p[i]

if rank == 0:
    print '\nTesting isomorphism of the two permutations of SR(25,12,5,6) graph 1:'
# verify using the inbuilt checking method
cert1, cert1p, checkISO = gi.isomorphicQ(cert1, cert1p)
if rank == 0:
    print checkISO

#=========================
#Two Non-Isomorphic Graphs
#=========================
if rank == 0:
    print '\nTesting isomorphism of graphs 1 and 11 of the SR(25,12,5,6) family:'

adj2 = np.genfromtxt('../graphs/strong-regular-25-12-5-6/11.txt')

cert1, cert2, checkISO  = gi.isomorphicQ(cert1, adj2)
if rank == 0:
    print checkISO
