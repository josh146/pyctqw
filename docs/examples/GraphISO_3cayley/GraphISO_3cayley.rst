============================
:file:`GraphISO_3cayley.py`
============================

Description 
-----------
This example script returns the GI certificates, and also checks isomorphism of various graphs in the strongly regular :math:`(25,12,5,6)` family.

1)  Firstly, two permutations of the 3-Cayley tree are formed by reordering the rows and columns in such a way as to preserve the symmetry of the graph. For example,

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

    The script then uses :func:`pyCTQW.MPI.GraphISO.AllIsomorphicQ` to test the isomorphism of all combinations of these two matrices (i.e. 1 with 2, 1 with 1, 2 with 2, and 2 with 1), and constructs a matrix containing :math:`1` for isomorphic graphs, and :math:`0` otherwise.

    As these graph *are* isomorphic, the result should be a matrix where all elements are :math:`1`.

2)  The second section of the program uses :func:`pyCTQW.MPI.GraphISO.AllIsomorphicQ` to test the isomorphism of all combinations of a 3-Cayley graph and a variant of the 3-Cayley graph, where an additional edge is placed between nodes 1 and 9.

    As these graph are *not* isomorphic, the result should be a :math:`2\times 2`: identity matrix.

Output
------------
.. code-block :: none
    :linenos:  

    $ mpirun -np 2 GraphISO_3cayley.py
    1) Testing isomorphism of all permutations:
    [[ 1.  1.]
     [ 1.  1.]]
    
    2) Testing isomorphism of all permutations:
    [[ 1.  0.]
     [ 0.  1.]]



Required Files
-----------------

    * :download:`3-cayley permutation 1 adjacency matrix </../graphs/cayley/3-cayley-permutations/3-caley-v1.txt>`
    * :download:`3-cayley permutation 2 adjacency matrix </../graphs/cayley/3-cayley-permutations/3-caley-v2.txt>`
    * :download:`3-cayley variant adjacency matrix </../graphs/cayley/3-cayley-variant/2.txt>`

Source Code
-------------

[:download:`Download source code </../examples/GraphISO_3cayley.py>`]

.. literalinclude:: /../examples/GraphISO_3cayley.py   
    :linenos:  





