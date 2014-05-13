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

# verify using the inbuilt checking method
if rank == 0:
    print('Testing isomorphism of K8 and K8_twisted01:')

cert1, cert2, checkISO = gi.isomorphicQ(adj1, adj2,
                            saveCert1='out/k8cert.txt',
                            saveCert2='out/k8cert_twisted01.txt')

if rank == 0:
    print(checkISO)
