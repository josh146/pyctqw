.. pyCTQW documentation master file, created by
   sphinx-quickstart on Wed Aug 28 23:14:36 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

----------------------------------
pyCTQW Manual
----------------------------------

.. only:: html


    +----------------------------+----------------------------------------------------------------------------------------------------------+
    |  **Contents**              |                                               **Features**                                               |
    +============================+==========================================================================================================+
    | .. toctree::               |                                                                                                          |
    |    :glob:                  |  * Fortran and Python bindings (in the form of a library and module respectively)                        |
    |    :maxdepth: 2            |  * Supports MPI through the use of the PETSc and SLEPc high-performance sparse                           |
    |                            |    matrix libraries (CUDA support planned)                                                               |
    |    downloads               |                                                                                                          |
    |    installation            |  * Has in-built support for infinite-line Hamiltonians                                                   |
    |    quickstart              |                                                                                                          |
    |    pyCTQW.MPI              |  * Can import custom adjancency matrices                                                                 |
    |    fortran                 |  * Supports one, two and three continuous-time quantum walkers (with or without interactions)            |
    |    examples                |  * Import and export matrices/states in binary or text format                                            |
    |                            |  * Python module supports plotting and visualisation using matplotlib and networkx                       |
    |                            |  * Entanglement calculations                                                                             |
    |                            |  * Ability to place diagonal defects/barriers on graph vertices                                          |
    |                            |  * MPI graph isomorphism methods, for both 2 and 3 interacting particles                                 |
    |                            |                                                                                                          |
    |                            +----------------------------------------------------------------------------------------------------------+
    |                            | **Gallery**                                                                                              |
    |                            +----------------------------------------------------------------------------------------------------------+
    |                            | .. image:: examples/2P_3cayley/2p_3cayley_nodes_particle1.png                                            |
    |                            |   :width: 300pt                                                                                          |
    |                            |                                                                                                          |
    |                            | .. image:: examples/1P_3cayley/1p_3cayley_graph.png                                                      |
    |                            |   :width: 300pt                                                                                          |
    |                            |                                                                                                          |
    |                            | .. image:: examples/1P_sr/1p_sr_graph.png                                                                |
    |                            |   :width: 300pt                                                                                          |
    |                            |                                                                                                          |
    |                            | .. image:: examples/2P_line/2p_line_plot.png                                                             |
    |                            |   :width: 300pt                                                                                          |
    |                            |                                                                                                          |
    |                            | .. image: examples/2P_3cayley/2p_3cayley_ent.png                                                         |
    |                            |   :width: 250pt                                                                                          |
    +----------------------------+----------------------------------------------------------------------------------------------------------+
    