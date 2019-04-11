#!/usr/bin/env python
import math
from math import *
import sys
import numpy as np
from numpy import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages
import os
plt.rcParams['font.size'] = 20
 
#Input specifications
natom = int(sys.argv[1])
#nbeads = int(sys.argv[2])
nbeads = 1
naxis = int(sys.argv[2])
gfact = float(sys.argv[3])
NumBeads = int(sys.argv[4])

ncol1 = nbeads + (natom-1)*3 # We have considered 0, M, N-1 beads
ncol = naxis + ncol1*3 # We have considered costheta, phi, chi
ncol = ncol+1
print(str(ncol+1)+'th column')
beadposition={0:'first-bead', 1:'middle-bead', 2:'last-bead'}

#Import the data files here
file1 = "/work/tapas/test-ratio-trick/ENT-RotDOFs-Rpt10.05Angstrom-gFactor"+str(gfact)+"-beta0.2Kinv-Blocks21000-Passes100-System2HF-ParticleA1-e0vsbeads-SWAPTOUNSWAP"+str(NumBeads)+"/results/output_instant.dof"

#Data processiong by Numpy 
x1 = loadtxt(file1, unpack=True, usecols=[ncol])

if (naxis != 0): 
	x1 = np.fabs(np.fmod(x1,2.0*pi))

#pyplot calling first
fig = plt.figure()
plt.grid(True)

font=28
purpose = "article"
if (purpose == "article"):
	font      = 28
elif (purpose == "ppt"):
	font      = 40
	rc('axes', linewidth=2)
fontlegend    = font/2.0

#X and Y labels
if (naxis != 0): 
	plt.xlabel('bins of '+r'$\phi$')
	fname='Phi'
else:
	plt.xlabel('bins of '+r'$\cos(\theta)$')
	fname='CosTheta'

plt.ylabel('Density')

# Normalize density
num_bins = 20
kwargs = dict(alpha=0.5, bins='auto', normed=True, stacked=True)

#Histogram plot
plt.hist(x1, **kwargs, color='b', label='""')

# Lebels for plotting
label1 = str(gfact)
label1 += ';  P = '+str(NumBeads)
print(label1)

# Tweak spacing to prevent clipping of ylabel
plt.gca().set(title='g = '+label1)
#plt.xlim(50,75)
#plt.legend();
#plt.subplots_adjust(left=0.15)
#plt.legend(bbox_to_anchor=(0.25, 2.20), loc=2, borderaxespad=1., shadow=True )
plt.subplots_adjust(top=0.92, bottom=0.14, left=0.18, right=0.98, hspace=0.0, wspace=0.0)

#Saving the plot in a file
outfile = "/home/tapas/ResultsOfENT/ENT-RotDOFs-Rpt10.05Angstrom-gFactor"+str(gfact)+"-beta0.2Kinv-Blocks21000-Passes100-System2HF-ParticleA1-e0vsbeads-SWAPTOUNSWAP"+str(NumBeads)+"-histgramOf"+fname+"-for-"+str(natom)+"th-rotor-and-"+beadposition[nbeads]+".eps"
print(outfile)
plt.savefig(outfile, dpi = 100, format = 'eps')
