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

natom = int(sys.argv[1])
naxis = int(sys.argv[2])

ncol = naxis + (natom-1)*3 # We have considered 0, M, N-1 beads
print(ncol)

num_bins = 50
file1 = "/Users/tsahoo/data-PIMC/PIMC-RotDOFs-Rpt10.05Angstrom-DipoleMoment2.0Debye-beta0.05Kinv-Blocks20000-Passes100-System2HF-e0vsbeads4_instant.dof"
file2 = "/Users/tsahoo/data-PIMC/PIMC-RotDOFs-Rpt10.05Angstrom-DipoleMoment2.0Debye-beta0.05Kinv-Blocks20000-Passes100-System2HF-e0vsbeads8_instant.dof"
file3 = "/Users/tsahoo/data-PIMC/PIMC-RotDOFs-Rpt10.05Angstrom-DipoleMoment2.0Debye-beta0.05Kinv-Blocks20000-Passes100-System2HF-e0vsbeads16_instant.dof"
file4 = "/Users/tsahoo/data-PIMC/PIMC-RotDOFs-Rpt10.05Angstrom-DipoleMoment2.0Debye-beta0.05Kinv-Blocks20000-Passes100-System2HF-e0vsbeads32_instant.dof"
file5 = "/Users/tsahoo/data-PIMC/PIMC-RotDOFs-Rpt10.05Angstrom-DipoleMoment2.0Debye-beta0.05Kinv-Blocks20000-Passes100-System2HF-e0vsbeads64_instant.dof"
#file2 = "/Users/tsahoo/mount_graham/PIMC-RotDOFs-Rpt10.05Angstrom-DipoleMoment2.0Debye-beta0.05Kinv-Blocks20000-Passes100-System2HF-e0vsbeads8/results/output_instant.dof"
#file3 = "/Users/tsahoo/mount_graham/PIMC-RotDOFs-Rpt10.05Angstrom-DipoleMoment2.0Debye-beta0.05Kinv-Blocks20000-Passes100-System2HF-e0vsbeads16/results/output_instant.dof"
#file4 = "/Users/tsahoo/mount_graham/PIMC-RotDOFs-Rpt10.05Angstrom-DipoleMoment2.0Debye-beta0.05Kinv-Blocks20000-Passes100-System2HF-e0vsbeads32/results/output_instant.dof"
#file5 = "/Users/tsahoo/mount_graham/PIMC-RotDOFs-Rpt10.05Angstrom-DipoleMoment2.0Debye-beta0.05Kinv-Blocks20000-Passes100-System2HF-e0vsbeads64/results/output_instant.dof"
x1 = loadtxt(file1, unpack=True, usecols=[ncol])
x2 = loadtxt(file2, unpack=True, usecols=[ncol])
x3 = loadtxt(file3, unpack=True, usecols=[ncol])
x4 = loadtxt(file4, unpack=True, usecols=[ncol])
x5 = loadtxt(file5, unpack=True, usecols=[ncol])
#x = np.fmod(x,2.0*pi)
#y = np.fabs(np.fmod(y,2.0*pi))

# Normalize
kwargs = dict(alpha=0.9, bins='auto', density=True, stacked=True)

# Plot
#plt.hist(x1, **kwargs, color='g', label='P=4')
plt.hist(x2, **kwargs, color='b', label='P=8')
#plt.hist(x3, **kwargs, color='y', label='P=16')
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
