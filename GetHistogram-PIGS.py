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

natom = int(sys.argv[1])
nbeads = int(sys.argv[2])
naxis = int(sys.argv[3])

ncol1 = nbeads + (natom-1)*3 # We have considered 0, M, N-1 beads
ncol = naxis + ncol1*3 # We have considered costheta, phi, chi
ncol = ncol+1

num_bins = 20
file1 = "/Users/tsahoo/linear-rotors-PIGS/PIGS-RotDOFs-Rpt10.05Angstrom-DipoleMoment8.0Debye-beta0.2Kinv-Blocks20000-Passes200-System2HF-e0vsbeads17/results/output_instant.dof"
#file2 = "/Users/tsahoo/linear-rotors-PIGS/PIGS-RotDOFs-Rpt10.05Angstrom-DipoleMoment3.0Debye-beta0.2Kinv-Blocks20000-Passes200-System2HF-e0vsbeads17/results/output_instant.dof"
#file3 = "/Users/tsahoo/linear-rotors-PIGS/PIGS-RotDOFs-Rpt10.05Angstrom-DipoleMoment3.0Debye-beta0.2Kinv-Blocks20000-Passes200-System2HF-e0vsbeads33/results/output_instant.dof"
x1 = loadtxt(file1, unpack=True, usecols=[ncol])
#x2 = loadtxt(file2, unpack=True, usecols=[ncol])
#x3 = loadtxt(file3, unpack=True, usecols=[ncol])
#x = np.fmod(x,2.0*pi)
#y = np.fabs(np.fmod(y,2.0*pi))

# Normalize
kwargs = dict(alpha=0.5, bins='auto', density=True, stacked=True)

# Plot
plt.hist(x1, **kwargs, color='g', label='P=5')
#plt.hist(x2, **kwargs, color='b', label='P=17')
#plt.hist(x3, **kwargs, color='y', label='P=33')
#plt.hist(x4, **kwargs, color='r', label='P=32')
#plt.hist(x5, **kwargs, color='g', label='P=64')
plt.gca().set(title='Probability Histogram of Diamond Depths', ylabel='Probability')
#plt.xlim(50,75)
plt.legend();

#plt.hist(x, bins='auto', density=True, facecolor='blue', alpha=0.7, label = "PIGS-ENT")
#plt.hist(y, bins='auto', normed=1, facecolor='green', alpha=0.7, label = "PIMC")
plt.xlabel('bins')
plt.ylabel('Probability Distribution')

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.legend(bbox_to_anchor=(0.35, 0.20), loc=2, borderaxespad=1., shadow=True )

outfile = "hist-test.eps"
plt.savefig(outfile, dpi = 200, format = 'eps')
plt.show()
#call(["open", outfile])
