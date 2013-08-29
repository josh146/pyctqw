#!/usr/bin/python
"""
Plotting functions.
"""
from petsc4py import PETSc as _PETSc
from matplotlib import pyplot as _plt
import numpy as _np
import pylab as _pl

import io as _io

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------- Plotting functions -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def plot(x,prob,savefile,t,init_state,d,amp,N,rank):

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

def plotGraph(ax,pos,lineX,lineY,prob=None,prob2=None,nodesize=None,barscale=1,output=None,
	barcolor='green',baralpha=0.25,barcolor2='blue',baralpha2=0.25,
	bartext=False,bartextcolor='black',bartextbg=None,btoffset=[-0.025,-0.025,0.05],
	bartext2=False,bartextcolor2='black',bartextbg2=None,btoffset2=[-0.025,-0.025,-0.05],
	nodecolor='red',nodealpha=0.5,
	nodetext=True,nodetextcolor='black',nodetextbg=None,nodetextbg2=None,ntoffset=[0,0,-0.15]):

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