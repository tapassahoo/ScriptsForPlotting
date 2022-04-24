import os
from collections import OrderedDict
from subprocess import call

import matplotlib as mpl
import matplotlib.axes
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.ticker as mticker
import numpy as np
from matplotlib.ticker import ScalarFormatter
from numpy import *
from pylab import *

# Global settings for plot
size=24
params = {'legend.fontsize': size*0.6,
    'figure.figsize': (8,6),
    'axes.labelsize': size,
    'axes.titlesize': size,
    'xtick.labelsize': size*0.75,
    'ytick.labelsize': size*0.75,
	'font.family': 'sans-serif',
	'mathtext.fontset':'dejavusans',
	'font.size': size*0.75,
    'axes.titlepad': size}
plt.rcParams.update(params)

# The directory where all the input files are stored
src_dir       = '/Users/tsahoo/academic-project/outputs/rotational_density_matrix/'

# Name of the inputs
file_spinless = src_dir+'rotational_propagator_spinless_linear_rotor.out'
file_para     = src_dir+'rotational_propagator_para_linear_rotor.out'
file_ortho    = src_dir+'rotational_propagator_ortho_linear_rotor.out'

# Reading data from the files
cost_spinless, rho_spinless = loadtxt(file_spinless, unpack=True, usecols=[0, 1])
cost_para, rho_para         = loadtxt(file_para, unpack=True, usecols=[0, 1])
cost_ortho, rho_ortho       = loadtxt(file_ortho, unpack=True, usecols=[0, 1])

# Plotting environment is created
fig = plt.figure()

# plot function is called to plot date files
plt.plot(cost_spinless, rho_spinless, color="black", ls="-", linewidth=2, marker="None", markersize=9, label="spinless")
plt.plot(cost_para, rho_para, color="blue", ls="--", linewidth=2, marker="None", markersize=9, label="para")
plt.plot(cost_ortho, rho_ortho, color="red", ls="-.", linewidth=2, marker="None", markersize=9, label="ortho")

# plt.ylim(0.0, 0.7)
# Tweaking spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.legend(bbox_to_anchor=(0.35, 0.20), loc=2, borderaxespad=2.0, shadow=True)

plt.ylabel(r'$\rho(\omega_{t, i}, \omega_{t+1, i}; \tau)$',labelpad=5)
plt.xlabel(r'$\cos (\gamma_{t, i})$', labelpad=5)
plt.subplots_adjust(top=0.98,bottom=0.14,left=0.16,right=0.99,hspace=0.0,wspace=0.0)
plt.legend(numpoints=1,loc='upper center')

#For ticks manipulating
plt.minorticks_on()
plt.tick_params(axis="both", direction="in", which="minor", right=False, top=False, length=2)
plt.tick_params(axis="both", direction="in", which="major", right=True, top=True, length=6)

outfile = src_dir+'Fig_rotational_propagator.pdf'
plt.savefig(outfile, dpi = 200, format = 'pdf')
plt.show()
