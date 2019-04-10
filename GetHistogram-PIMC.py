#!/usr/bin/env python
import math
from math import *
import sys
import numpy as np
from numpy import *
import matplotlib
#matplotlib.use('pdf')
import matplotlib.pyplot as plt
from itertools import cycle
import itertools
from scipy.optimize import curve_fit
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import ScalarFormatter
import support_without_parallel as support
import os
 
import matplotlib.mlab as mlab
fig = plt.figure()
plt.grid(True)

#Input specifications
natom = int(sys.argv[1])
naxis = int(sys.argv[2])
ncol = naxis + (natom-1)*3 # We have considered costheta, phi, chi
print(ncol)

#Import the data files here
#file1 = "/Users/tsahoo/linear-rotors-PIMC/PIMC-RotDOFs-Rpt10.05Angstrom-DipoleMoment1.0Debye-beta0.05Kinv-Blocks20000-Passes200-System2HF-e0vsbeads4/results/output_instant.dof"
file2 = "/Users/tsahoo/linear-rotors-PIMC/PIMC-RotDOFs-Rpt10.05Angstrom-DipoleMoment4.0Debye-beta0.05Kinv-Blocks20000-Passes200-System2HF-e0vsbeads32/results/output_instant.dof"

#Data processiong by Numpy 
#x1 = loadtxt(file1, unpack=True, usecols=[ncol])
x2 = loadtxt(file2, unpack=True, usecols=[ncol])

if (naxis != 0): 
	x2 = np.fabs(np.fmod(x2,2.0*pi))

# Normalize desity
num_bins = 20
kwargs = dict(alpha=0.5, bins='auto', density=True, stacked=True)

# Lebels for plotting
label1 = r'$\mu = 2.0$'+'Debye'
label2 = r'$\mu = 4.0$'+'Debye'
label3 = r'$\mu = 6.0$'+'Debye'
label4 = r'$\mu = 8.0$'+'Debye'

#Histogram plot
#plt.hist(x1, **kwargs, color='r', label=label1)
plt.hist(x2, **kwargs, color='b', label=label2)
#plt.hist(x3, **kwargs, color='y', label=label3)
#plt.hist(x4, **kwargs, color='r', label=label4)
#plt.hist(x5, **kwargs, color='g', label='P=64')

#X and Y labels
if (naxis != 0): 
	plt.xlabel('bins of '+r'$\phi$')
else:
	plt.xlabel('bins of '+r'$\cos(\theta)$')

plt.ylabel('Density')

# Tweak spacing to prevent clipping of ylabel
#plt.gca().set(title='Probability Histogram of Diamond Depths', ylabel='Probability')
#plt.xlim(50,75)
plt.legend();
plt.subplots_adjust(left=0.15)
plt.legend(bbox_to_anchor=(0.35, 0.20), loc=2, borderaxespad=1., shadow=True )

#Saving the plot in a file
outfile = "histgramOfCosTheta.eps"
plt.savefig(outfile, dpi = 100, format = 'eps')
plt.show()
