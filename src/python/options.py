#!/usr/bin/env python2.7
import os
import argparse
import ConfigParser
from math import *
import re

def parse_args():
	# program info
	progname = "PyCTQW"
	libname = "libpyctqw"
	configFile = "input.pyqw"
	progversion = "PyCTQW v0.1\n Josh Izaac\n josh.izaac@uwa.edu.au"
	info = '{0} is a Python interface to the Fortran library {1}.\
		\nThis library contains functions that allow continuous-time\
		\nquantum walks (CTQWs) to be simulated fast and efficiently.\
		\n\nIn particular, pyctqw allows you to:\
			\n\t * Propagate a single particle CTQW on an\
			\n\t   infinite line in the presence of defects\
			\n\t * Propagate a two particle non-interacting\
			\n\t   CTQW on an infinite line in the presence of defects'
			 	 
	# look for config file		 	 
	conf_parser = argparse.ArgumentParser(add_help=False)

	conf_parser.add_argument("-c", "--conf",
		help="specify input file", metavar="FILE")
		
	args, remaining_argv = conf_parser.parse_known_args()
	
	# config settings
	config = ConfigParser.SafeConfigParser()
	if args.conf:
		config.read([args.conf])
		defaults = dict(config.items("DEFAULTS"))
		
	elif os.path.exists(configFile):
		config.read([configFile])
		
		# add configuration to the dictionary, and convert to proper type
		defaults = {}
		for sec in config.sections():
			for option,value in config.items(sec):
				if value.lower()=='true' or value.lower()=='false':
					defaults[option] = config.getboolean(sec,option)
				else:
					try:
						defaults[option] = eval(value)
					except:
						defaults[option] = value
								
	else:
		# No config defaults
		defaults =	{ "output" 	: "./out",
				  "statespace"	: False,
				  "particles"	: 2,
				  "grid_length"	: 50,
				  "time"	: 1.0,
				  "expm"	: 'chebyshev',
				  "evalgen"	: False,
				  "evallib"	: 'lapack',
				  "defect_nodes": [0],
				  "defect_amp"	: [0],
				  "p1_initial_state": [(0,1/sqrt(2)), (1,1j/sqrt(2))],
				  "p2_initial_state": [(0,0,1/sqrt(2)), (1,1,1/sqrt(2))]
				}
	
	# command line arguments
	parser = argparse.ArgumentParser(parents=[conf_parser],description=info.format(progname,libname),
			version=progversion, add_help=True, formatter_class=argparse.RawDescriptionHelpFormatter)
	
	parser.set_defaults(**defaults)
	
	parser.add_argument('-o', dest='output', metavar='DIR',
		help='specify an output directory', type=str)
	
	parser.add_argument('-s', "--statespace", action="store_true",
		help="specify whether the statespace is exported")
	
	parser.add_argument('-p', '--particles',
		help='set the number of quantum walkers',			
		type=int, metavar='1|2', nargs=1)
		
	parser.add_argument('-t', '--time',
		help='QW propagation time',			
		type=float, metavar='T', nargs=1)
	
	parser.add_argument('-N','--grid-length',
		help='set the number of vertices (must be an even & > 2)',			
		type=int, metavar='NUMBER', nargs=1)

	parser.add_argument('-d', '--defect-nodes',
		help='location of defects',			
		type=float, metavar='D', nargs='+')
	
	parser.add_argument('-a', '--defect-amp',
		help='defects amplitudes',			
		type=float, metavar='A', nargs='+')
	
	parser.add_argument('-2p0', '--p2-initial-state',
		help='initial state of the quantum walker',			
		type=float, metavar='(X,Y,C)', nargs='+')
	
	parser.add_argument('-1p0', '--p1-initial-state',
		help='initial state of the quantum walker',			
		type=float, metavar='(J,C)', nargs='+')
		
	parser.add_argument('--expm', dest='expm',
		help='choose the algorithm used to calculate the matrix exponential,\
		\'Chebyshev\' (default), or \'Burkadt\'',			
		type=str, metavar='METHOD', nargs=1)
		
	parser.add_argument('-evg', '--eval-gen', dest='evalgen', action="store_true",
		help='when calculating the eigenvalues, force the linear algebra library\
		to treat the Hamiltonian as a general complex matrix, rather than as a\
		real, symmetric	band matrix (the default)')
		
	parser.add_argument('--evallib', dest='evallib',
		help='choose the library used to calculate the eigenvalues,\
		\'lapack\' (default) or \'numpy\'. Note that lapack is significantly\
		faster, however requires either intel-mkl or LAPACK be installed',			
		type=str, metavar='LIB', nargs=1)

	return parser.parse_args(remaining_argv)
	
	
	
	
	
	
