#!/usr/bin/env python
#PBS -W group_list=director565
#PBS -q workq
#PBS -l walltime=00:20:00
#PBS -l select=20:ncpus=12:mem=64gb
#PBS -j oe
#PBS -o /dev/null
#PBS -e /dev/null

import os
import sys
import re
import shutil
import subprocess

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~Set properties here~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
program = 'mpirun -np 2 example -N '
environment = 'intel'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def module_add(modulename):
	p = subprocess.Popen("/usr/bin/modulecmd python add " + modulename, 
		stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	stdout,stderr = p.communicate()
	exec stdout

# Checking to see if running on PBS or not
if (os.getenv("PBS_ENVIRONMENT") == None):
	print "in a shell"
	workdir = os.getcwd()
	myEnv = os.environ.copy()
else:
	print "in PBS"
	workdir = os.getenv('PBS_O_WORKDIR')
	arrayIndex = os.getenv('PBS_ARRAY_INDEX')
	jobid = os.getenv('PBS_JOBID')
	
	# set additional environment variables
	myEnv = os.environ.copy()
	if environment == 'gcc':
		FORNAX_MODULES = 'gcc atlas openmpi cuda valgrind'
		petsc_arch = 'linux-gnu'
	else:
		FORNAX_MODULES = 'intel intel-mkl openmpi cuda valgrind'
		petsc_arch = 'linux-gnu-intel'
		
	#load modules
	module_add(FORNAX_MODULES)
	print os.environ["LOADEDMODULES"]
	
	myEnv = os.environ.copy()
	myEnv['PETSC_ARCH'] = petsc_arch

# change to working directory and find input script
os.chdir(workdir)

# run the simulation
for i in range(100000,1000000,10000):
	program2 = program + str(i)
	simulation = subprocess.Popen(program2.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=myEnv)
	stdout,stderr = simulation.communicate()
	print(stdout)
	print >> sys.stderr, stderr
