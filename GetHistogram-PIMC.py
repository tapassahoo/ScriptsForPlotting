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
naxis = int(sys.argv[2])
DipoleMoment = float(sys.argv[3])
NumBeads = int(sys.argv[4])

ncol = naxis + (natom-1)*3 # We have considered costheta, phi, chi
print(str(ncol+1)+'th column')

#Import the data files here
file1 = "/work/tapas/linear-rotors-PIMC/PIMC-RotDOFs-Rpt10.05Angstrom-DipoleMoment"+str(DipoleMoment)+"Debye-beta0.05Kinv-Blocks20000-Passes200-System2HF-e0vsbeads"+str(NumBeads)+"/results/output_instant.dof"

#Data processiong by Numpy 
x1 = loadtxt(file1, unpack=True, usecols=[ncol])

if (naxis != 0): 
	x1 = np.fabs(np.fmod(x1,2.0*pi))

#pyplot calling first
fig = plt.figure()
plt.grid(True)

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
kwargs = dict(alpha=0.5, bins=20, normed=True, stacked=True)

#Histogram plot
plt.hist(x1, **kwargs, color='b', label='""')

# Lebels for plotting
label1 = str(DipoleMoment)+' Debye'
label1 += ';  P = '+str(NumBeads)
print(label1)

# Tweak spacing to prevent clipping of ylabel
plt.gca().set(title=r'$\mu = $'+label1)
#plt.xlim(50,75)
#plt.legend();
#plt.subplots_adjust(left=0.15)
#plt.legend(bbox_to_anchor=(0.25, 2.20), loc=2, borderaxespad=1., shadow=True )
plt.subplots_adjust(top=0.92, bottom=0.14, left=0.18, right=0.98, hspace=0.0, wspace=0.0)

#Saving the plot in a file
outfile = "/home/tapas/ResultsOfPIMC/PIMC-RotDOFs-Rpt10.05Angstrom-DipoleMoment"+str(DipoleMoment)+"Debye-beta0.05Kinv-Blocks20000-Passes200-System2HF-e0vsbeads"+str(NumBeads)+"-histgramOf"+fname+"-for-"+str(natom)+"th-rotor.eps"
print(outfile)
plt.savefig(outfile, dpi = 100, format = 'eps')
