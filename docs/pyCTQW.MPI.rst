==========================
:mod:`pyCTQW.MPI` Package
==========================

.. module:: pyCTQW.MPI

Contains distributed memory methods for calculating:
	* the continuous-time quantum walk on arbitrary graphs and lines for 1, 2 and 3 particles
	* Non-interacting and interacting walks
	* Entanglement calculations
	* ability to place diagonal defects/barriers on graph vertices
	* contains MPI graph isomorphism methods, for both 2 and 3 interacting particles

These classes also contain plotting methods, to visualise the quantum walking dynamics

.. currentmodule:: pyCTQW.MPI

Arbitrary graph QW
--------------------------

.. autosummary::
	:toctree:

	Graph
	Graph2P
	Graph3P

.. seealso::
	Examples using these objects and methods:
		:doc:`examples/1P_3cayley/1P_3cayley`,
		:doc:`examples/2P_3cayley/2P_3cayley`,
		:doc:`examples/3P_3cayley/3P_3cayley`


Infinite line QW
--------------------------

.. autosummary::
	:toctree:

	Line
	Line2P
	Line3P

.. seealso::
	Examples using these objects and methods:
		:doc:`examples/1P_line/1P_line`,
		:doc:`examples/2P_line/2P_line`,
		:doc:`examples/3P_line/3P_line`


.. .. important::
.. 	*	only **even** number nodes `N` are permitted for the :class:`Line` classes.
.. 	*	unlike the :class:`Graph` classes, which label the graph nodes \
.. 		from :math:`j\in\{0,\dots,N-1\}`, the :class:`Line` classes \
.. 		label the nodes on the line from :math:`j\in\{1-N/2,\dots,N/2\}`.


Graph Isomorphism
--------------------------

.. autosummary::
	:toctree:

	GraphISO

.. seealso::
	Examples using this object:
		:doc:`examples/GraphISO_sr/GraphISO_sr`,
		:doc:`examples/GraphISO_3cayley/GraphISO_3cayley`


Submodules
--------------------------

.. caution::
	The methods and classes available in these submodules have been
	designed to work with the classes listed above; that is, they are
	being called implicitly by the above classes, invisible to the user.

	It is not recommended that they be used by themselves, however it is
	possible to do so.


.. autosummary::
	:toctree:

	pyCTQW.MPI.io
	pyCTQW.MPI.plot
	pyCTQW.MPI.ctqw