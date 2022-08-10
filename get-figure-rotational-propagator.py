#!/usr/bin/env python
""" 
This code is used for generating contour plot. 
"""
__author__ = "Dr. Tapas Sahoo"
__contact__ = "tapascuchem@gmail.com"
__copyright__ = "Copyright 2022, GITHUB"
__credits__ = ["Dr. Tapas Sahoo"]
__date__ = "2022/08/06"
__deprecated__ = False
__email__ =  "tapascuchem@gmail.com"
__license__ = "GPLv3"
__maintainer__ = "Dr. Tapas Sahoo"
__status__ = "Production"
__version__ = "0.0.1"

import os
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--verbose", help="increase output verbosity", action="store_true")
args = parser.parse_args()

# Global setup for each figure
size=24
params = {'legend.fontsize': size*0.6,
    'figure.figsize': (8,6),
    'axes.labelsize': size,
    'axes.titlesize': size,
    'xtick.labelsize': size*0.75,
    'ytick.labelsize': size*0.75,
    'font.size': size*0.75,
    'axes.titlepad': size}
plt.rcParams.update(params)
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

# Reading the file, we would like to plot.
file_str="den030_120K_spinless"
home=os.path.expanduser("~")
inputFile=home+"/academic-project/outputs/rotational-propagator-for-water/rho."+file_str
rho = np.genfromtxt(inputFile, unpack=True, usecols=[3])

# Data file generation for the contour plot.
size_grid = 361
theta=np.arange(361)
phi=theta
theta2d, phi2d = np.meshgrid(theta, phi)
rho2d=np.reshape(rho, (size_grid, size_grid))

# Plot env is started here.
fig,ax=plt.subplots()
cp = ax.contourf(theta2d, phi2d, rho2d)
fig.colorbar(cp, pad=0.01)

# Set axes labeling
ax.set_xlabel(r'$\phi$', labelpad=10)
ax.set_ylabel(r'$\chi$', labelpad=10)

# Ticks manipulating
tick_position = np.arange(0, 361, 90)
tick_label = [r"$0$", r"$\dfrac{\pi}{2}$", r"$\pi$", r"$\dfrac{3\pi}{2}$", r"$2\pi$"]

ax.set_xticks(tick_position)
ax.set_xticklabels(tick_label)
ax.set_yticks(tick_position)
ax.set_yticklabels(tick_label)

# Adjustment the plot
plt.subplots_adjust(top=0.98,bottom=0.2,left=0.13,right=0.98,hspace=0.0,wspace=0.0)

# Saving the figure
outputFile=home+"/academic-project/outputs/rotational-propagator-for-water/Figure-"+file_str+".pdf"
plt.savefig(outputFile, format='pdf')
plt.show()
