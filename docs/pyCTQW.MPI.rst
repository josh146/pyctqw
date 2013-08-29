==========================
:mod:`pyCTQW.MPI` Glossary
==========================

Contains distributed memory methods for calculating:
	* the continuous-time quantum walk on arbitrary graphs and lines for 1, 2 and 3 particles
	* Non-interacting and interacting walks
	* Entanglement calculations
	* ability to place diagonal defects/barriers on graph vertices
	* contains MPI graph isomorphism methods, for both 2 and 3 interacting particles

These classes also contain plotting methods, to visualise the quantum walking dynamics

.. currentmodule:: pyCTQW.MPI

.. rst-class:: html-toggle

.. rst-class:: emphasize-children

Arbitrary graph CTQW
--------------------------

.. autosummary::
	:toctree:

	Graph
	Graph2P
	Graph3P

.. rst-class:: html-toggle

.. rst-class:: emphasize-children

Infinite line CTQW
--------------------------

.. autosummary::
	:toctree:

	Line
	Line2P
	Line3P


.. important::
	*	only **even** number nodes `N` are permitted for the :class:`Line` classes.
	*	unlike the :class:`Graph` classes, which label the graph nodes \
		from :math:`j\in\{0,\dots,N-1\}`, the :class:`Line` classes \
		label the nodes on the line from :math:`j\in\{1-N/2,\dots,N/2\}`.


.. rst-class:: emphasize-children

.. rst-class:: html-toggle

Graph Isomorphism
--------------------------

.. autosummary::
	:toctree:

	GraphISO

.. rst-class:: emphasize-children

.. rst-class:: html-toggle

Submodules
--------------------------

.. caution::
	The methods and classes available in these submodules have been
	designed to work with the classes listed above; that is, they are
	being called implicitly by the above classes, invisible to the user.
	Therefore, it is not recommended that they be used by themselves, however it is
	possible to do so.

.. autosummary::
	:toctree:

	pyCTQW.MPI.io
	pyCTQW.MPI.plot
	pyCTQW.MPI.ctqw