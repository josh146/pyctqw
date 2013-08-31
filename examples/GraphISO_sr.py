#!/usr/bin/env python2.7
# initialize PETSc
import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np

# import pyCTQW as qw
import pyCTQW.MPI as qw

# get the MPI rank
rank =  PETSc.Comm.Get_rank(PETSc.COMM_WORLD)

# create a graph isomorphism object
gi = qw.GraphISO()

#=========================
#  Two Isomorphic Graphs
#=========================

# import two adjacency matrices; these are permutations
# of the first graph in the SR-25-12-5-6, and thus are isomorphic
adj1 = np.genfromtxt('../graphs/strong-regular-25-12-5-6/1-permutations/SR1_perm_1.txt')
adj2 = np.genfromtxt('../graphs/strong-regular-25-12-5-6/1-permutations/SR1_perm_3.txt')

# generate certificates for these two graphs
cert1 = gi.GIcert(adj1)
cert2 = gi.GIcert(adj2)

# print the certificates side by side;
# they should be identical
if rank == 0:
	for i in range(len(cert1)):
		print cert1[i], cert2[i]

# verify using the inbuilt checking method
checkISO = gi.isomorphicQ(adj1,adj2)
if rank == 0: print checkISO

#=========================
#Two Non-Isomorphic Graphs
#=========================

adj1 = np.genfromtxt('../graphs/strong-regular-25-12-5-6/1.txt')
adj2 = np.genfromtxt('../graphs/strong-regular-25-12-5-6/11.txt')

checkISO = gi.isomorphicQ(adj1,adj2)
if rank == 0: print checkISO