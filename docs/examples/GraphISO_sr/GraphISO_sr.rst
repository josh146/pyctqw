======================
:file:`GraphISO_sr.py`
======================

Description 
-----------
This example script returns the GI certificates, and also checks isomorphism of various graphs in the strongly regular :math:`(25,12,5,6)` family.

1)  Firstly, two permutations of graph 1 in the SR family are formed by reordering the rows and columns in such a way as to preserve the symmetry of the graph. For example,

    .. math::
        :nowrap:

        \begin{align}
        &\text{i)}\text{ swap rows }2\text{ & }3: \left[\begin{matrix}
                a & b & c\\
                b & d & e\\
                c & e & f
            \end{matrix}\right] \Rightarrow
         \left[\begin{matrix}
                b & d & e\\
                a & b & c\\
                c & e & f
            \end{matrix}\right]\\[12pt]

         &\text{ii)} \text{ swap columns }2\text{ & }3:   \left[\begin{matrix}
                b & d & e\\
                a & b & c\\
                c & e & f
            \end{matrix}\right] \Rightarrow
            \left[\begin{matrix}
                d & b & e\\
                b & a & c\\
                e & c & f
            \end{matrix}\right]
        \end{align}


    -- for adjacency matrices, this is process is equivalent to 'swapping' the label of nodes 1 and 2, and thus is isomorphic to the original adjacency matrix.

    The GI certificates of both permutations of graph 1 are then calculated and displayed side by side. For each graph, these consist of a column of probabilities resulting from performing an interacting two particle quantum walk over *all* possible initial bosonic edge states; i.e.

    .. math::
        |\psi(0)\rangle = \frac{1}{\sqrt{2}}\left(|j\rangle\otimes|j\rangle + |k\rangle\otimes|k\rangle  \right)~~~~\forall j,k\in\{0,1,\dots,N-1\},

    coupled with the frequency with which they appear. The quantum walk is performed for time :math:`t=2N`, allowing the walk time to propagate to all nodes. Therefore, as the graphs are isomorphic, the GI certificates should be identical.

    Finally, :func:`pyCTQW.MPI.GraphISO.isomorphicQ` is used to verify that the two graphs are in fact isomorphic.

2)  The second section of the program uses :func:`pyCTQW.MPI.GraphISO.isomorphicQ` to check whether graphs 1 and 11 of the SR family are isomorphic. They are not, and this is correctly returned by the program.

Output
------------
.. code-block :: none
    :linenos:  

    $ mpirun -np 2 GraphISO_sr.py

    The GI certificates of two permutations of SR(25,12,5,6) graph 1:
    [  1.14200170e-01   1.88140000e+04] [  1.14200170e-01   1.88140000e+04]
    [  1.53741042e-01   1.54380000e+04] [  1.53741042e-01   1.54380000e+04]
    [  1.76558859e-01   1.02600000e+04] [  1.76558859e-01   1.02600000e+04]
    [  1.28264401e-01   9.89800000e+03] [  1.28264401e-01   9.89800000e+03]
    [  1.02738389e-01   9.75000000e+03] [  1.02738389e-01   9.75000000e+03]
    [  2.22889534e-01   5.45200000e+03] [  2.22889534e-01   5.45200000e+03]
    [  1.43460570e-01   5.32600000e+03] [  1.43460570e-01   5.32600000e+03]
    [  2.08976427e-01   4.86000000e+03] [  2.08976427e-01   4.86000000e+03]
    [  1.88334629e-01   4.25600000e+03] [  1.88334629e-01   4.25600000e+03]
    [  9.23774555e-02   2.46800000e+03] [  9.23774555e-02   2.46800000e+03]
    [  1.98534713e-01   1.96000000e+03] [  1.98534713e-01   1.96000000e+03]
    [  1.64488997e-01   1.79000000e+03] [  1.64488997e-01   1.79000000e+03]
    [  2.38430321e-01   1.76200000e+03] [  2.38430321e-01   1.76200000e+03]
    [  2.52977006e-01   5.34000000e+02] [  2.52977006e-01   5.34000000e+02]
    [  2.65853983e-01   4.04000000e+02] [  2.65853983e-01   4.04000000e+02]
    [  7.43487849e-02   3.01000000e+02] [  7.43487849e-02   3.01000000e+02]
    [   0.27720404  178.        ] [   0.27720404  178.        ]

    Testing isomorphism of the two permutations of SR(25,12,5,6) graph 1:
    True

    Testing isomorphism of graphs 1 and 11 of the SR(25,12,5,6) family:
    False

Required Files
-----------------

    * :download:`SR1_perm_1 adjacency matrix </../graphs/strong-regular-25-12-5-6/1-permutations/SR1_perm_1.txt>`
    * :download:`SR1_perm_3 adjacency matrix</../graphs/strong-regular-25-12-5-6/1-permutations/SR1_perm_3.txt>`
    * :download:`1.txt adjacency matrix </../graphs/strong-regular-25-12-5-6/1.txt>`
    * :download:`11.txt adjacency matrix </../graphs/strong-regular-25-12-5-6/11.txt>`

Source Code
-------------

[:download:`Download source code </../examples/GraphISO_sr.py>`]

.. literalinclude:: /../examples/GraphISO_sr.py   
    :linenos:  





