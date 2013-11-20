.. _example-MPI:

======================
:file:`exampleMPI.F90`
======================

Description 
-----------

This Fortran example propagates a 1 particle continuous-time quantum walk on an infinite line

Amongst the features used, it illustrates:
	*   calling the :f:mod:`libctqwMPI` library from Fortran
	*   recieving command line options using PETSc 
	*   creating PETSc log stages
	*   the use of the chebyshev algorithm
	*   adding a diagonal defects to various nodes
	*   viewing vectors via PetscViewer

Compiling and executing
------------------------

To compile this example, open a terminal in the top directory of the ``pyCTQW-X.Y`` folder, and run

.. code-block:: bash

	$ make fortran
	$ make examples

The example can then be run by changing into the :file:`examples/fortran` folder, and running

.. code-block:: bash

	$ mpirun -np X exampleMPI

where ``X`` is the number of MPI processes to run.

Output
------------

.. code-block:: bash

	$ mpirun -np 4 exampleMPI
	Vector Object: 4 MPI processes
	  type: mpi
	Process [0]
	1.30318e-15 + 9.41288e-16 i
	4.26061e-15 - 5.83926e-15 i
	-2.6771e-14 - 1.97231e-14 i
	-8.90705e-14 + 1.19677e-13 i
	5.22435e-13 + 3.92973e-13 i
	1.69261e-12 - 2.2254e-12 i
	-9.24365e-12 - 7.11278e-12 i
	-2.91412e-11 + 3.74127e-11 i
	1.47434e-10 + 1.16317e-10 i
	4.51971e-10 - 5.6522e-10 i
	-2.10616e-09 - 1.70822e-09 i
	-6.27416e-09 + 7.62086e-09 i
	2.67489e-08 + 2.23733e-08 i
	7.73784e-08 - 9.09727e-08 i
	-2.99429e-07 - 2.59261e-07 i
	-8.40536e-07 + 9.52525e-07 i
	2.92435e-06 + 2.63331e-06 i
	7.96061e-06 - 8.65084e-06 i
	-2.46145e-05 - 2.31843e-05 i
	-6.49345e-05 + 6.72304e-05 i
	0.000175879 + 0.000174556 i
	0.00044937 - 0.000439577 i
	-0.00104657 - 0.00110507 i
	-0.00258837 + 0.0023656 i
	0.00505632 + 0.00575514 i
	Process [1]
	0.012099 - 0.0101714 i
	-0.0191447 - 0.0239345 i
	-0.0442924 + 0.0334726 i
	0.053856 + 0.0761079 i
	0.12025 - 0.0787511 i
	-0.102846 - 0.172366 i
	-0.219755 + 0.116923 i
	0.111189 + 0.241323 i
	0.214786 - 0.0830223 i
	-0.0457964 - 0.132734 i
	-0.0208922 + 0.0294486 i
	0.0624386 - 0.0589731 i
	-0.0404472 - 0.137791 i
	-0.191971 - 0.0860889 i
	-0.229877 + 0.138719 i
	-0.0383627 + 0.243064 i
	0.0596756 + 0.217395 i
	0.21466 + 0.189964 i
	0.261424 + 0.00431375 i
	0.21683 - 0.0655021 i
	0.178224 - 0.166243 i
	0.0818483 - 0.179087 i
	0.0335323 - 0.183573 i
	-0.01719 - 0.145368 i
	-0.0340458 - 0.120739 i
	Process [2]
	-0.0344425 - 0.0737224 i
	-0.0328421 - 0.0407838 i
	0.0103582 - 0.0348961 i
	-0.0119036 + 0.0257429 i
	-0.00586558 - 0.0254129 i
	0.010251 - 0.000336933 i
	-0.0147693 - 0.0158695 i
	0.0286581 - 0.0184521 i
	-0.00508709 + 0.00712859 i
	0.0298955 - 0.0310717 i
	0.026666 + 0.0413678 i
	-0.011159 - 0.00816756 i
	0.0425462 + 0.0359954 i
	-0.0610777 + 0.0473539 i
	-0.018832 - 0.0447773 i
	-0.00248435 + 0.0247191 i
	-0.0610373 - 0.0547772 i
	0.0902541 - 0.0779196 i
	0.0756292 + 0.10105 i
	-0.0915928 + 0.061552 i
	-0.0438766 - 0.0715585 i
	0.049748 - 0.0280784 i
	0.0163906 + 0.0313776 i
	-0.0181928 + 0.00882673 i
	-0.0044223 - 0.00978976 i
	Process [3]
	0.00492519 - 0.00207491 i
	0.000916545 + 0.00233018 i
	-0.0010417 + 0.000382831 i
	-0.000151757 - 0.000441799 i
	0.000178365 - 5.72714e-05 i
	2.06318e-05 + 6.875e-05 i
	-2.53648e-05 + 7.11154e-06 i
	-2.35018e-06 - 8.97787e-06 i
	3.05473e-06 - 7.45989e-07 i
	2.27794e-07 + 1.00097e-06 i
	-3.16392e-07 + 6.70093e-08 i
	-1.90127e-08 - 9.66128e-08 i
	2.85389e-08 - 5.20876e-09 i
	1.37913e-09 + 8.16537e-09 i
	-2.26541e-09 + 3.53182e-10 i
	-8.75365e-11 - 6.10116e-10 i
	1.59661e-10 - 2.10083e-11 i
	4.88366e-12 + 4.06353e-11 i
	-1.00669e-11 + 1.0998e-12 i
	-2.39917e-13 - 2.42956e-12 i
	5.71638e-13 - 5.06798e-14 i
	1.03592e-14 + 1.31214e-13 i
	-2.94057e-14 + 2.04654e-15 i
	-3.89668e-16 - 6.42454e-15 i
	1.43612e-15 - 7.35772e-17 i



Source Code
--------------------------------------------------------
[:download:`Download source code </../examples/fortran/exampleMPI.F90>`]

.. literalinclude:: /../examples/fortran/exampleMPI.F90
    :linenos:

