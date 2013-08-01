#!/usr/bin/python
import sys, petsc4py
petsc4py.init(sys.argv)

from petsc4py import PETSc
from libctqwMPI import ctqwmpi
OptDB = PETSc.Options()

N = OptDB.getInt('N', 6)
lambda_ = OptDB.getReal('lambda', 6.0)
do_plot = OptDB.getBool('plot', False)

H = PETSc.Mat()
H.create(PETSc.COMM_WORLD)
mat_H = H.fortran

psi0 = PETSc.Vec()
psi0.create(PETSc.COMM_WORLD)
psi0.setSizes(N)
psi0.setUp()

psi0.setValues([0,1,2,3,4,5],[0.01,0.02,1,2,5,0])
psi0.assemble()

psi = psi0.duplicate()

d = [-2,1]
amp = [5.0,1.0]
ctqwmpi.hamiltonian_1p_line(mat_H,d,amp,N)
print H.view()

rank = PETSc.Comm.Get_rank(PETSc.COMM_WORLD)

ctqwmpi.qw_cheby(psi0.fortran,psi.fortran,5.0,mat_H,0.496405,7.20008,rank,N)

psi.view()

H.destroy()
psi.destroy()
psi0.destroy()
