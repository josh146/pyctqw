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
	
	# config defaults
	defaults =	{ "input_state"	: "",
			  "output" 	: "./out",
			  "fortran"	: "intel",
			  "sparse"	: False,
			  "statespace"	: False,
			  "particles"	: 1,
			  "grid_length"	: 160,
			  "time"	: 30.0,
			  "expm"	: 'chebyshev',
			  "eig_method"	: False,
			  "eig_lib"	: 'lapack',
			  "defect_nodes": [0],
			  "defect_amp"	: [0],
			  "p1_initial_state": [(0,1)],
			  "p2_initial_state": [(0,5,1/sqrt(2))]
			}
			 	 
	# look for config file		 	 
	conf_parser = argparse.ArgumentParser(add_help=False)

	conf_parser.add_argument("-c", "--conf",
		help="specify input file", metavar="FILE")
		
	args, remaining_argv = conf_parser.parse_known_args()
	
	# config settings
	config = ConfigParser.SafeConfigParser()
	if args.conf or os.path.exists(configFile):
		if args.conf:
			print "Input script " + args.conf + " found!"
			config.read([args.conf])
		elif os.path.exists(configFile):
			print "Input script " + configFile + " found!"
			config.read([configFile])
			
		# add configuration to the dictionary, and convert to proper type
		for sec in config.sections():
			for option,value in config.items(sec):
				if value.lower()=='true' or value.lower()=='false':
					defaults[option] = config.getboolean(sec,option)
				else:
					try:
						defaults[option] = eval(value)
					except:
						defaults[option] = value

	# command line arguments
	parser = argparse.ArgumentParser(parents=[conf_parser],description=info.format(progname,libname),
			version=progversion, add_help=True, formatter_class=argparse.RawDescriptionHelpFormatter)
	
	parser.set_defaults(**defaults)
	
	parser.add_argument('-i', '--input-state', metavar='FILE',
		help='specify an input state', type=str)
	
	parser.add_argument('-o', dest='output', metavar='DIR',
		help='specify an output directory', type=str)
		
	parser.add_argument('-fc', '--fortran',
		help='specify the fortran library to use (either \'intel\' or \'gcc\')', type=str)
	
	parser.add_argument('-s', "--statespace", action="store_true",
		help="specify whether the statespace is exported")
	
	parser.add_argument('-p', '--particles',
		help='set the number of quantum walkers',			
		type=int, metavar='1|2')
		
	parser.add_argument('-t', '--time',
		help='QW propagation time',			
		type=float, metavar='T')
	
	parser.add_argument('-N','--grid-length',
		help='set the number of vertices (must be an even & > 2)',			
		type=int, metavar='NUMBER')

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
		\'chebyshev\' (default, fortran), \'burkadt\' (fortran), or \'python-chebyshev\'',			
		type=str, metavar='METHOD')
		
	parser.add_argument('-eg', '--eig-general', action="store_true",
		help='when calculating the eigenvalues, force the linear algebra library\
		to treat the Hamiltonian as a general complex matrix, rather than as a\
		real, symmetric	band matrix (the default)')
		
	parser.add_argument('--sparse', action="store_true",
		help='use SciPy and ARPACK to create sparse matrices. Note that this overrides\
		the \'expm\'  and \'eig-lib\' setting, as it requires expm=\'python-chebyshev\'\
		and eig-lib=\'scipy-arpack\'')
		
	parser.add_argument('--eig-lib',
		help='choose the library used to calculate the eigenvalues,\
		\'lapack\' (default), \'python-scipy\' (with ARPACK bindings)\
		or \'python-numpy\'. Note that lapack is significantly\
		faster, however requires either intel-mkl or LAPACK be installed',			
		type=str, metavar='LIB')

	return parser.parse_args(remaining_argv)
	
	
	
	
	
	
