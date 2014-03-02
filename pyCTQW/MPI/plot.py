#!/usr/bin/python
#
#  This file provides a variety of plotting functions.
#  ----------------------------------------------------------------------------
#  pyCTQW - Distributed memory CTQW Fortran library and Python module
#  Copyright (C) 2013-2014, Joshua Izaac
#
#  pyCTQW is free software: you can redistribute it and/or modify it under the
#  terms of version 3 of the GNU General Public License as published by
#  the Free Software Foundation.
#
#  pyCTQW is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
#  more details.
#
#  You should have received a copy of the GNU General Public License
#  along with pyCTQW. If not, see <http://www.gnu.org/licenses/>.
#  ----------------------------------------------------------------------------
from petsc4py import PETSc as _PETSc
from matplotlib import pyplot as _plt
import numpy as _np
import pylab as _pl

import io as _io

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------- Plotting functions -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def plot(x,prob,savefile,t,init_state,d,amp,N,rank):
	""" Creates a plot of probability vs graph node for a specified time.

		Args:
			N (int) : the number of nodes to be plotted
			x (array or numpy.array of ints) : array of length :math:`N` containing the nodes to be plotted
			prob (petsc4py.PETSc.Vec) : vector of length :math:`N` containing probabilities for each node in ``x``.
			savefile (str) : the absolute/relative path to the desired output file.
			t (float) : the time of the measurement
			init_state (array) : an :math:`n\\times 2` array, containing the initial \
						state of the quantum walker in the format ``[[j1,amp1],[j2,amp2],...]``.
			d (array of ints): an array containing *integers* indicating the nodes
						where diagonal defects are to be placed (e.g. ``d=[0,1,4]``).
			amp (array of floats):   an array containing *floats* indicating the diagonal defect
						amplitudes corresponding to each element in ``d`` (e.g. ``amp=[0.5,-1,4.2]``).
			rank (int): the rank of each MPI node calling this function.

		Warning:
			* The size of ``amp`` and ``d`` must be identical
			* ensure a file extension is present so that filetype is correctly set\
			(choose one of png, pdf, ps, eps or svg).
		"""

	def prob_plot_p1(probVec,savefile,t,initstate,d,a,N):
		# convert vectors to arrays
		prob = _np.real(_np.asarray(probVec))
	
		# create plot
		fig = _plt.figure()
		_plt.plot(x, prob)
		_plt.ylabel("$|\langle j|\psi\\rangle|^2$",rotation=90)
		_plt.grid(b=None, which='major', axis='both', linestyle='-', alpha=0.3)
		_plt.xlabel("$j$")
	
		# Plot titles
		if initstate[0]=='f':
			disp = "\nInitial state: {0}"
			IS_disp = [initstate]
		else:
			IS_disp = ()
			disp = "\nInitial state: $|\psi(0)\\rangle="

			for i in range(len(initstate)):
				IS_disp = IS_disp + ("({1: .3g})|{0.real: .0f}\\rangle".format(*(initstate[i])),)
	
				if i == 0:
					disp = disp + "{" + str(i) + "} "
				else:
					disp = disp + "+ {" + str(i) + "} "	
		
		_plt.suptitle("CTQW probability distribution at time $t={}$".format(t))
	
		if (len(list(set(a))) == 1) and (list(set(a))[0] == 0.0):
			_plt.title(disp.format(*IS_disp) + "$",
				multialignment='left', fontsize=11)
		else:
			def_disp = "\nDefects: $"
			for i in range(len(d)):
				def_disp += "{1: .3}|{0}\\rangle +".format(d[i],a[i])
		
			_plt.title(disp.format(*IS_disp) + "$" + def_disp[:-2] + "$",
				multialignment='left', fontsize=11)

		# save plot
		_plt.subplots_adjust(top=0.85)
		_pl.savefig(savefile)
		_plt.close()
		
	# scatter prob to process 0
	commX = prob.getComm()
	scatterX, prob0 = _PETSc.Scatter.toZero(prob)
	scatterX.scatter(prob, prob0, False, _PETSc.Scatter.Mode.FORWARD)
	
	# use process 0 to create the plot
	if rank==0:
		prob_plot_p1(prob0,savefile,t,init_state,d,amp,N)
	
	# deallocate	
	commX.barrier()
	scatterX.destroy()
	prob0.destroy()

def plot2P(x,psiX,psiY,savefile,t,init_state,d,amp,N,rank):
	""" Creates a plot of probability vs graph node for two probability distributions at a specified time.

		Args:
			N (int) : the number of nodes to be plotted
			x (array or numpy.array of ints) : array of length :math:`N` containing the nodes to be plotted
			psiX (petsc4py.PETSc.Vec) : vector of length :math:`N` containing probabilities for each node in ``x``.
			psiY (petsc4py.PETSc.Vec) : vector of length :math:`N` containing probabilities for each node in ``x``.
			savefile (str) : the absolute/relative path to the desired output file.
			t (float) : the time of the measurement
			init_state (array) : an :math:`n\\times 3` array, containing the initial \
						state of the quantum walker in the format ``[[x1,y1,amp1],[x2,y2,amp2],...]``.
			d (array of ints): an array containing *integers* indicating the nodes
						where diagonal defects are to be placed (e.g. ``d=[0,1,4]``).
			amp (array of floats):   an array containing *floats* indicating the diagonal defect
						amplitudes corresponding to each element in ``d`` (e.g. ``amp=[0.5,-1,4.2]``).
			rank (int): the rank of each MPI node calling this function.

		Warning:
			* The size of ``a`` and ``d`` must be identical
			* ensure a file extension is present so that filetype is correctly set\
			(choose one of png, pdf, ps, eps or svg).
		"""
	
	def prob_plot_p2(psiX,psiY,savefile,t,initstate,d,a,N):
	
		# convert vectors to arrays
		probX = _np.real(_np.asarray(psiX))
		probY = _np.real(_np.asarray(psiY))

		# create plot
		lbl = [0,1];
	
		fig = _plt.figure()
		lbl[0], = _plt.plot(x, probX, 'b-')
		lbl[1], = _plt.plot(x, probY, 'r--')
	
		_plt.ylabel("$|\langle j|\psi\\rangle|^2$",rotation=90)
		_plt.grid(b=None, which='major', axis='both', linestyle='-', alpha=0.3)
		_plt.xlabel("$j$")
	
		# legend
		leg = _plt.legend(lbl, ['P1','P2'], loc='center right', bbox_to_anchor=(1.015, 0.93))
		_plt.setp(leg.get_texts(), fontsize='small') 
	
		# Plot titles
		if initstate[0]=='f':
			disp = "\nInitial state: {0}"
			IS_disp = [initstate]
		else:
			IS_disp = ()
			disp = "\nInitial state: $|\psi(0)\\rangle="

			for i in range(len(initstate)):
				IS_disp = IS_disp + ("({2: .3g})|{0.real: .0f},{1.real: .0f}\\rangle".format(*(initstate[i])),)
	
				if i == 0:
					disp = disp + "{" + str(i) + "} "
				else:
					disp = disp + "+ {" + str(i) + "} "	
	
		_plt.suptitle("CTQW probability distribution at time $t={}$".format(t))
	
		if (len(list(set(a))) == 1) and (list(set(a))[0] == 0.0):
			_plt.title(disp.format(*IS_disp) + "$",
				multialignment='left', fontsize=11)
		else:
			def_disp = "\nDefects: $"
			for i in range(len(d)):
				def_disp += "{1: .3}|{0}\\rangle +".format(d[i],a[i])	
		
			_plt.title(disp.format(*IS_disp) + "$" +  def_disp[:-2] + "$",
				multialignment='left', fontsize=11)

		# save plot
		#_plt.ylim((0,0.3))
		_plt.subplots_adjust(top=0.85)
		_pl.savefig(savefile)
		_plt.close()
	
	# scatter psiX to process 0
	commX = psiX.getComm()
	scatterX, psiX0 = _PETSc.Scatter.toZero(psiX)
	scatterX.scatter(psiX, psiX0, False, _PETSc.Scatter.Mode.FORWARD)
	
	# scatter psiY to process 0
	commY = psiY.getComm()
	scatterY, psiY0 = _PETSc.Scatter.toZero(psiY)
	scatterY.scatter(psiY, psiY0, False, _PETSc.Scatter.Mode.FORWARD)
	
	# use process 0 to create the plot
	if rank==0:
		prob_plot_p2(psiX0,psiY0,savefile,t,init_state,d,amp,N)
	
	# deallocate	
	commX.barrier()
	scatterX.destroy()
	psiX0.destroy()
	
	commY.barrier()
	scatterY.destroy()
	psiY0.destroy()

def plot3P(x,psiX,psiY,psiZ,savefile,t,init_state,d,amp,N,rank):
	""" Creates a plot of probability vs graph node for three probability distributions at a specified time.

		Args:
			N (int) : the number of nodes to be plotted
			x (array or numpy.array of ints) : array of length :math:`N` containing the nodes to be plotted
			psiX (petsc4py.PETSc.Vec) : vector of length :math:`N` containing probabilities for each node in ``x``.
			psiY (petsc4py.PETSc.Vec) : vector of length :math:`N` containing probabilities for each node in ``x``.
			psiZ (petsc4py.PETSc.Vec) : vector of length :math:`N` containing probabilities for each node in ``x``.
			savefile (str) : the absolute/relative path to the desired output file.
			t (float) : the time of the measurement
			init_state (array) : an :math:`n\\times 3` array, containing the initial \
						state of the quantum walker in the format ``[[x1,y1,amp1],[x2,y2,amp2],...]``.
			d (array of ints): an array containing *integers* indicating the nodes
						where diagonal defects are to be placed (e.g. ``d=[0,1,4]``).
			amp (array of floats):   an array containing *floats* indicating the diagonal defect
						amplitudes corresponding to each element in ``d`` (e.g. ``amp=[0.5,-1,4.2]``).
			rank (int): the rank of each MPI node calling this function.

		Warning:
			* The size of ``a`` and ``d`` must be identical
			* ensure a file extension is present so that filetype is correctly set\
			(choose one of png, pdf, ps, eps or svg).
		"""
	
	def prob_plot_p3(psiX,psiY,psiZ,savefile,t,initstate,d,a,N):
	
		# convert vectors to arrays
		probX = _np.real(_np.asarray(psiX))
		probY = _np.real(_np.asarray(psiY))
		probZ = _np.real(_np.asarray(psiZ))

		# create plot
		lbl = [0,1,2];
	
		fig = _plt.figure()
		lbl[0], = _plt.plot(x, probX, 'b-')
		lbl[1], = _plt.plot(x, probY, 'r--')
		lbl[2], = _plt.plot(x, probZ, 'g-.')
	
		_plt.ylabel("$|\langle j|\psi\\rangle|^2$",rotation=90)
		_plt.grid(b=None, which='major', axis='both', linestyle='-', alpha=0.3)
		_plt.xlabel("$j$")
	
		# legend
		leg = _plt.legend(lbl, ['P1','P2','P3'], loc='center right', bbox_to_anchor=(1.015, 0.93))
		_plt.setp(leg.get_texts(), fontsize='small') 
	
		# Plot titles
		if initstate[0]=='f':
			disp = "\nInitial state: {0}"
			IS_disp = [initstate]
		else:
			IS_disp = ()
			disp = "\nInitial state: $|\psi(0)\\rangle="

			for i in range(len(initstate)):
				IS_disp = IS_disp + ("({3: .3g})|{0.real: .0f},{1.real: .0f},{2.real: .0f}\\rangle".format(*(initstate[i])),)
	
				if i == 0:
					disp = disp + "{" + str(i) + "} "
				else:
					disp = disp + "+ {" + str(i) + "} "	
	
		_plt.suptitle("CTQW probability distribution at time $t={}$".format(t))
	
		if (len(list(set(a))) == 1) and (list(set(a))[0] == 0.0):
			_plt.title(disp.format(*IS_disp) + "$",
						fontsize=11,multialignment='left')
		else:
			def_disp = "\nDefects: $"
			for i in range(len(d)):
				def_disp += "{1: .3}|{0}\\rangle +".format(d[i],a[i])	
		
			_plt.title(disp.format(*IS_disp) + "$" +  def_disp[:-2] + "$",
							fontsize=11,multialignment='left')

		# save plot
		#_plt.ylim((0,0.3))
		_plt.subplots_adjust(top=0.85)
		_pl.savefig(savefile)
		_plt.close()
	
	# scatter psiX to process 0
	commX = psiX.getComm()
	scatterX, psiX0 = _PETSc.Scatter.toZero(psiX)
	scatterX.scatter(psiX, psiX0, False, _PETSc.Scatter.Mode.FORWARD)
	
	# scatter psiY to process 0
	commY = psiY.getComm()
	scatterY, psiY0 = _PETSc.Scatter.toZero(psiY)
	scatterY.scatter(psiY, psiY0, False, _PETSc.Scatter.Mode.FORWARD)

	# scatter psiZ to process 0
	commZ = psiY.getComm()
	scatterZ, psiZ0 = _PETSc.Scatter.toZero(psiZ)
	scatterY.scatter(psiZ, psiZ0, False, _PETSc.Scatter.Mode.FORWARD)
	
	# use process 0 to create the plot
	if rank==0:
		prob_plot_p3(psiX0,psiY0,psiZ0,savefile,t,init_state,d,amp,N)
	
	# deallocate	
	commX.barrier()
	scatterX.destroy()
	psiX0.destroy()
	
	commY.barrier()
	scatterY.destroy()
	psiY0.destroy()

	commZ.barrier()
	scatterZ.destroy()
	psiZ0.destroy()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------- Graph plot functions -----------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getGraphNodes(adj,layout='spring'):
	""" Returns arrays containing cartesian coordinates of the nodes/connections \
		contained in the input adjacency matrix.

		Args:
			adj (petsc4py.PETSc.Vec) : an :math:`N\\times N` PETSc-type adjacency matrix.
			layout (str):  the format to store the position of the nodes (only used
							when running :func:`plotGraph`).

							* ``spring`` *(default)* - spring layout.
							* ``circle`` - nodes are arranged in a circle.
							* ``spectral`` - nodes are laid out according to the \
											spectrum of the graph.
							* ``random`` - nodes are arranged in a random pattern.

		:rtype:	tuple of arrays

		.. important::
			Requires `NetworkX <http://networkx.github.io/>`_

		.. admonition:: Example

				>>> nodePos, lineX, lineY = pyCTQW.MPI.plots.getGraphNodes(adj,layout='spring')

			where
				* :attr:`nodePos` contains the :math:`(x,y)` coordinates of the vertices
				* :attr:`lineX` contains the :math:`x` coordinates of edges connecting vertices
				* :attr:`lineY` contains the :math:`y` coordinates of edges connecting vertices
		"""

	try:
		import networkx as _nx
	except:
		print '\nNetworkX Python module required for graph plotting!'
		return

	graph = _nx.from_scipy_sparse_matrix(_io.matToSparse(adj).real)
	
	if layout == 'circle':
		pos = _nx.circular_layout(graph)
	elif layout == 'spectral':
		pos = _nx.spectral_layout(graph)
	elif layout == 'random':
		pos = _nx.random_layout(graph)
	elif layout == 'shell':
		pos = _nx.shell_layout(graph)
	else:
		pos = _nx.spring_layout(graph,dim=2)
	
	testpos = []
	for i in pos.itervalues():
		testpos.append(i.tolist())
	testpos = _np.array(testpos)

	lineX = []
	lineY = []
	for i in _nx.edges(graph):
		lineX.append([testpos[i[0]][0], testpos[i[1]][0]])
		lineY.append([testpos[i[0]][1], testpos[i[1]][1]])
		
	return testpos, lineX, lineY

def plotGraph(ax,pos,lineX,lineY,prob=None,prob2=None,nodesize=None,barscale=1,
	barcolor='green',baralpha=0.25,barcolor2='blue',baralpha2=0.25,
	bartext=False,bartextcolor='black',bartextbg=None,btoffset=[-0.025,-0.025,0.05],
	bartext2=False,bartextcolor2='black',bartextbg2=None,btoffset2=[-0.025,-0.025,-0.05],
	nodecolor='red',nodealpha=0.5,
	nodetext=True,nodetextcolor='black',nodetextbg=None,nodetextbg2=None,ntoffset=[0,0,-0.15]):

	""" Creates a plot of probability vs node superimposed on a 3D visualisation of the graph vertices.

			Args:
				ax (mpl_toolkits.mplot3d.plt3d.Axes3D): a matplotlib 3D axes to plot the graph on.
				pos (array) : :math:`(x,y)` coordinates of the graph vertices.
				lineX (array): :math:`x` coordinates of graph edges connecting vertices.
				lineY (array): :math:`y` coordinates of graph edges connecting vertices.
				probP (petsc4py.PETSc.Vec) : the probability to be represented
								as bars placed above/below each graph vertex for ``P=(1|2)``
								repectively. If ``None``, the graph	is plotted without any probability.
			
			Keyword Args:
				nodesize (float) :  size of the vertices in the plot. If left blank, this is
									determined automatically.
				nodecolor (str)  :  vertex color (*default* ``'red'``).

									For more details on how to specify a color,
									see the `matplotlib documentation
									<http://matplotlib.org/api/colors_api.html#module-matplotlib.colors>`_.

				nodealpha (float)  :  value between 0 and 1 specifying the vertex opacity (*default* 0.25)

				nodetext (bool) :  if set ``True``, the vertices are labelled by number.
				nodetextcolor (str) : vertex label text color (*default* ``'black'``).
				nodetextbg (str) : vertex label background color (*default* ``'None'``).                
				ntofffset (array of floats): the :math:`(x,y,z)` vertex label offset relative \
											 to the vertex (*default* ``[0.,0.,-0.15]``).

				barscaleP (float) :  scaled height of the probability bars (*default* ``1``).
				barcolorP (str)   :  probability bar color (*default* ``'green'``).
				baralphaP (float) : value between 0 and 1 specifying the opacity (*default* 0.25)
				bartext (bool) :  if set ``True``, the probability bars are labelled with their value.
				bartextcolorP (str) : probability label text color (*default* ``'black'``).
				bartextbgP (str) : probability label background color (*default* ``'None'``).
				btoffsetP (array of floats): the :math:`(x,y,z)` probability label offset relative to the top of
											the probability bars (*default* ``[-0.025,-0.025,0.05]``)

			.. important::
				* Where an argument ends with ``P`` above, substitute ``P=(1|2)`` to \
					access the property for particle 1 and 2 respectively.

				* Only MPI node **0** will run this function; all others will just ``Pass``. \
					**Thus if the** :attr:`probP` **vectors are distributed over multiple nodes, \
					they must first be gathered to node 0**

				* This function does not initialize a 3D matplotlib figure/axes; that **must** be done \
					independently, and passed through the :attr:`ax` argument.
			"""


	# use process 0 to create the plot
	rank = _PETSc.COMM_WORLD.Get_rank()
	if rank==0:

		from matplotlib.colors import ColorConverter as cl
		from matplotlib.patches import Circle
		import mpl_toolkits.mplot3d.art3d as art3d

		for i in range(len(lineX)):
			ax.plot(lineX[i], lineY[i],zs=-0.01,color='black',alpha=0.8, linewidth=2)
			
		#ax.scatter3D(pos.T[0], pos.T[1], color = 'orange', marker = "o", s=200)
		#ax.bar3d([x-0.04 for x in pos.T[0]],[x-0.04 for x in pos.T[1]],
		#         _np.zeros_like(pos.T[0]),0.08,0.08,0, color='orange',alpha=0.3,edgecolor='gray')
		
		if nodesize is None:
			nodesize=[]
			for i in range(len(lineX)):
			   nodesize.append(_np.sqrt(_np.sum(_np.square([-_np.subtract(*lineX[i]), -_np.subtract(*lineY[i])]))))
			
			nodesize=min(_np.min(nodesize)*0.6/2,0.05)
		
		for i in range(len(pos)):
			p = Circle((pos.T[0][i], pos.T[1][i]), nodesize, color=nodecolor, alpha=nodealpha)
			ax.add_patch(p)
			art3d.pathpatch_2d_to_3d(p, z=0.0, zdir="z")
			
			if nodetext:
				ax.text(pos.T[0][i]+ntoffset[0], pos.T[1][i]+ntoffset[1],ntoffset[2],str(i),color=nodetextcolor,
						bbox=(None if nodetextbg is None else dict(facecolor=nodetextbg, alpha=0.2)))
		
		if prob is not None:
			# probability bars
			for i,val in enumerate(prob):
				if val != 0:
					ax.bar3d(pos.T[0][i]-0.7*nodesize/2,pos.T[1][i]-0.7*nodesize/2,0,0.7*nodesize,0.7*nodesize,val*barscale,
							 color=barcolor,edgecolor=cl.to_rgba(cl(),barcolor,alpha=baralpha),alpha=baralpha)
					
					if bartext:
						ax.text(pos.T[0][i]+btoffset[0], pos.T[1][i]+btoffset[1],val*barscale+btoffset[2],'{0: .4f}'.format(val),
								color=bartextcolor,bbox=(None if bartextbg is None else dict(facecolor=nodetextbg, alpha=0.1)))
		
		if prob2 is not None:
			# probability bars
			for i,val in enumerate(prob2):
				if val != 0:
					ax.bar3d(pos.T[0][i]-0.7*nodesize/2,pos.T[1][i]-0.7*nodesize/2,0,0.7*nodesize,0.7*nodesize,-val*barscale,
							 color=barcolor2,edgecolor=cl.to_rgba(cl(),barcolor2,alpha=baralpha2),alpha=baralpha2)
					
					if bartext:
						ax.text(pos.T[0][i]+btoffset[0], pos.T[1][i]+btoffset[1],-val*barscale-btoffset[2],'{0: .4f}'.format(val),
								color=bartextcolor2,bbox=(None if bartextbg2 is None else dict(facecolor=nodetextbg2, alpha=0.1)))

		
		if prob2 is None:
			ax.set_zlim3d([0,1])
		else:
			ax.set_zlim3d([-1,1])

		ax.set_xlim3d([pos.T[0].min(),pos.T[0].max()])
		ax.set_ylim3d([pos.T[1].min(),pos.T[1].max()])
		ax.set_axis_off()

def plotEntanglement(t,entArray,savefile,initstate,d,a):
	""" Creates a plot of entanglement vs time.

		Args:
			t (array of floats) : array containing the time data.
			entArray (array or numpy.array of floats) : array of length matching :attr:`t` \
											containing the entanglement data to be plotted.
			savefile (str) : the absolute/relative path to the desired output file.
			initstate (array) : an :math:`n\\times 3` array, containing the initial \
						state of the quantum walker in the format ``[[x1,y1,amp1],[x2,y2,amp2],...]``.
			d (array of ints): an array containing *integers* indicating the nodes
						where diagonal defects are to be placed (e.g. ``d=[0,1,4]``).
			a (array of floats):   an array containing *floats* indicating the diagonal defect
						amplitudes corresponding to each element in ``d`` (e.g. ``amp=[0.5,-1,4.2]``).

		Warning:
			* The size of ``a`` and ``d`` must be identical
			* ensure a file extension is present so that filetype is correctly set\
			(choose one of png, pdf, ps, eps or svg).
		"""

	# create plot
	fig = _plt.figure()
	_plt.plot(t, entArray)
	_plt.ylabel("Entanglement",rotation=90)
	_plt.grid(b=None, which='major', axis='both', linestyle='-', alpha=0.3)
	_plt.xlabel("$t$")

	# Plot titles
	if initstate[0]=='f':
		disp = "\nInitial state: {0}"
		IS_disp = [initstate]
	else:
		IS_disp = ()
		disp = "\nInitial state: $|\psi(0)\\rangle="

		for i in range(len(initstate)):
			IS_disp = IS_disp + ("({2: .3g})|{0.real: .0f},{1.real: .0f}\\rangle".format(*(initstate[i])),)

			if i == 0:
				disp = disp + "{" + str(i) + "} "
			else:
				disp = disp + "+ {" + str(i) + "} "	

	_plt.suptitle("CTQW two particle entanglement")

	if (len(list(set(a))) == 1) and (list(set(a))[0] == 0.0):
		_plt.title(disp.format(*IS_disp) + "$",
			multialignment='left', fontsize=11)
	else:
		def_disp = "\nDefects: $"
		for i in range(len(d)):
			def_disp += "{1: .3}|{0}\\rangle +".format(d[i],a[i])	
	
		_plt.title(disp.format(*IS_disp) + "$" +  def_disp[:-2] + "$",
			multialignment='left', fontsize=11)

	# save plot
	_plt.subplots_adjust(top=0.85)
	_pl.savefig(savefile)
	_plt.close()

def plotNodes(time,nodes,probArray,savefile,p=0):
	""" Creates a plot of probability vs time for specified nodes.

		Args:
			time (array or numpy.array of floats) : array containing the nodes to be plotted
			nodes (array) : the node numbers to be plotted.
			probArray (numpy.array of floats) : 2D array containing probabilities over time for each node.
			savefile (str) : the absolute/relative path to the desired output file.
			p (int) : (``0``-``3``) - the particle that is being plotted. If ``p=0``, then
						it is assumed that this is a 1 particle system.

		Warning:
			* ensure a file extension is present so that filetype is correctly set\
			(choose one of png, pdf, ps, eps or svg).
		"""
	from itertools import cycle
	
	nodeNum = probArray.shape[0]

	# create plot
	lbl = range(nodeNum)
	plotStyles = cycle(["-","--","-.",":"])

	fig = _plt.figure()
	for i in range(nodeNum):
		lbl[i], = _plt.plot(time, probArray[i], next(plotStyles))

	_plt.xlabel("$t$")
	_plt.ylabel("$|\langle j|\psi\\rangle|^2$",rotation=90)
	_plt.grid(b=None, which='major', axis='both', linestyle='-', alpha=0.3)

	_plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)

	# legend
	leg = _plt.legend(lbl, [str(x) for x in nodes], loc='lower center', bbox_to_anchor=(0.5, -0.3),ncol=8, fancybox=True)
	_plt.setp(leg.get_texts(), fontsize='small') 

	if p == 0:
		_plt.title("CTQW probability distribution over time")
	else:
		_plt.title("Particle {} CTQW probability distribution over time".format(p))

	# save plot
	#_plt.ylim((0,0.3))
	_pl.savefig(savefile)
	_plt.close()

def plotNodes2P(time,node,probXArray,probYArray,savefile):
	""" Creates a 2 particle plot of probability vs time for specified nodes.

		Args:
			time (array or numpy.array of floats) : array containing the nodes to be plotted
			nodes (array) : the node numbers to be plotted.
			probXArray (numpy.array of floats) : 2D array containing particle 1 probabilities \
												over time for each node.
			probYArray (numpy.array of floats) : 2D array containing particle 2 probabilities \
												over time for each node.
			savefile (str) : the absolute/relative path to the desired output file.

		Warning:
			* ensure a file extension is present so that filetype is correctly set\
			(choose one of png, pdf, ps, eps or svg).
		"""
	from itertools import cycle

	# create plot
	lbl = range(2)
	plotStyles = cycle(["-","--","-.",":"])

	fig = _plt.figure()
	lbl[0], = _plt.plot(time, probXArray, next(plotStyles))
	lbl[1], = _plt.plot(time, probYArray, next(plotStyles))

	_plt.xlabel("$t$")
	_plt.ylabel("$|\langle j|\psi\\rangle|^2$",rotation=90)
	_plt.grid(b=None, which='major', axis='both', linestyle='-', alpha=0.3)

	_plt.subplots_adjust(left=0.1, right=0.85, top=0.9, bottom=0.1)

	# legend
	leg = _plt.legend(lbl, ['P1','P2'], loc='lower center', bbox_to_anchor=(1.1, 0.5), fancybox=True)
	_plt.setp(leg.get_texts(), fontsize='small') 

	_plt.title("CTQW probability of node {} over time".format(node))

	# save plot
	#_plt.ylim((0,0.3))
	_pl.savefig(savefile)
	_plt.close()

def plotNodes3P(time,node,probXArray,probYArray,probZArray,savefile):
	""" Creates a 3 particle plot of probability vs time for specified nodes.

		Args:
			time (array or numpy.array of floats) : array containing the nodes to be plotted
			nodes (array) : the node numbers to be plotted.
			probXArray (numpy.array of floats) : 2D array containing particle 1 probabilities \
												over time for each node.
			probYArray (numpy.array of floats) : 2D array containing particle 2 probabilities \
												over time for each node.
			probZArray (numpy.array of floats) : 2D array containing particle 3 probabilities \
												over time for each node.
			savefile (str) : the absolute/relative path to the desired output file.

		Warning:
			* ensure a file extension is present so that filetype is correctly set\
			(choose one of png, pdf, ps, eps or svg).
		"""
	from itertools import cycle

	# create plot
	lbl = range(3)
	plotStyles = cycle(["-","--","-.",":"])

	fig = _plt.figure()
	lbl[0], = _plt.plot(time, probXArray, next(plotStyles))
	lbl[1], = _plt.plot(time, probYArray, next(plotStyles))
	lbl[2], = _plt.plot(time, probZArray, next(plotStyles))

	_plt.xlabel("$t$")
	_plt.ylabel("$|\langle j|\psi\\rangle|^2$",rotation=90)
	_plt.grid(b=None, which='major', axis='both', linestyle='-', alpha=0.3)

	_plt.subplots_adjust(left=0.1, right=0.85, top=0.9, bottom=0.1)

	# legend
	leg = _plt.legend(lbl, ['P1','P2','P3'], loc='lower center', bbox_to_anchor=(1.1, 0.5), fancybox=True)
	_plt.setp(leg.get_texts(), fontsize='small') 

	_plt.title("CTQW probability of node {} over time".format(node))

	# save plot
	#_plt.ylim((0,0.3))
	_pl.savefig(savefile)
	_plt.close()