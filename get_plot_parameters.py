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

import matplotlib.pyplot as plt

# Global setup for each figure
def get_plt_params():
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
