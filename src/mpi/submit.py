#!/usr/bin/env python
#PBS -W group_list=intern2012
#PBS -q workq
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=12:mem=64gb
#PBS -j oe

import os
import sys
import re
import shutil
import subprocess

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~Set properties here~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
program = 'mpirun example -n 1000 -log_summary'
FORNAX_MODULES = 'gcc atlas openmpi cuda valgrind'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def module_add(modulename):
	p = subprocess.Popen("/usr/bin/modulecmd python add " + modulename, 
		stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	stdout,stderr = p.communicate()
	exec stdout

#load modules
module_add(FORNAX_MODULES)
print os.environ["LOADEDMODULES"]

# Checking to see if running on PBS or not
if (os.getenv("PBS_ENVIRONMENT") == None):
	print "in a shell"
else:
	print "in PBS"
	workdir = os.getenv('PBS_O_WORKDIR')
	arrayIndex = os.getenv('PBS_ARRAY_INDEX')
	jobid = os.getenv('PBS_JOBID')

# change to working directory and find input script
os.chdir(workdir)

# compile mumax2
simulation = subprocess.Popen(program.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout,stderr = simulation.communicate()
print(stdout)
print >> sys.stderr, stderr
