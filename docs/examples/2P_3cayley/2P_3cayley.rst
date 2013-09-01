=====================
:file:`2P_3cayley.py`
=====================

Description
------------

.. |3-cayley-graph| image:: 3cayley.png
    :width: 300pt

+-----------------------+-------------------------------------------------------------------------------------+
|                       |This example propagates a 2 particle continuous-time quantum walk on a 3-cayley tree.|
|   |3-cayley-graph|    |                                                                                     |
|                       |Amongst the features used, it illustrates:                                           |
|                       |    *   recieving command line options using PETSc                                   |
|                       |    *   the use of the chebyshev algorithm, with minimum eigenvalue pre-set          |
|                       |    *   adding a diagonal defect to various nodes                                    |
|                       |    *   same-node interactions between particles                                     |
|                       |    *   creating node handles to watch the probability at specified nodes            |
|                       |    *   creating entanglement handles to watch the entanglement                      |
|                       |    *   various plotting abilities:                                                  |
|                       |            - probability vs node plots                                              |
|                       |            - probability vs time plots                                              |
|                       |            - graph plots                                                            |
|                       |            - entanglement plots                                                     |
|                       |    * outputs the partial trace of the density matrix in binary form                 |
+-----------------------+-------------------------------------------------------------------------------------+

Output
------------

.. image:: 2p_3cayley_graph.png
    :width: 360pt

.. image:: 2p_3cayley_plot.png
    :width: 360pt

.. image:: 2p_3cayley_nodes_particle1.png
    :width: 360pt

.. image:: 2p_3cayley_nodes_particle2.png
    :width: 360pt

.. image:: 2p_3cayley_node1.png
    :width: 360pt

.. image:: 2p_3cayley_ent.png
    :width: 360pt

Required Files
-----------------
    * :download:`Adjacency matrix </../graphs/cayley/3-cayley.txt>`

Source Code
--------------------------------------------------------
[:download:`Download source code </../examples/2P_3cayley.py>`]

.. literalinclude:: /../examples/2P_3cayley.py
    :linenos:

