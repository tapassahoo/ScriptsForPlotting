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

# Import Library
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
import get_plot_parameters as pm

parser = argparse.ArgumentParser()
parser.add_argument("--verbose", help="increase output verbosity", action="store_true")
args = parser.parse_args()

# Import the function for the global parameters used in all the figure
pm.plot_parameters()

# Reading the file, we would like to plot.
home=os.path.expanduser("~")

# Define data
theta_array = np.linspace(0, 2.0*np.pi, 51)
phi_array=np.array([0, np.pi/2])
angle_dict = {0:"0",1:"90"}
#
# The generic function for the dipole-dipole interaction potential function
for count, phi in enumerate(phi_array):
	raw_output_file=home+"/academic-project/output/final-pigs-outputs-for-plotting/raw-data-dipole-dipole-interaction-potential-for-phi-" + angle_dict[count] + "-degree.txt"
	file = open(raw_output_file,"w")
	for theta1 in theta_array:
		for theta2 in theta_array:
			func = np.sin(theta1)*np.sin(theta2)*np.cos(phi)-2.0*np.cos(theta1)*np.cos(theta2)
			file.write(str(theta1) + "  " + str(theta2) + "   " + str(func) + "\n")
	file.close()

# Define data
theta1 = np.linspace(0, 2.0*np.pi, 201)
theta2 = np.linspace(0, 2.0*np.pi, 201)
theta1_new = theta1[:,np.newaxis]
theta2_new = theta2[np.newaxis,:]
#
phi_array=np.array([0, np.pi/2])
#
# The generic function for the dipole-dipole interaction potential function
for count, phi in enumerate(phi_array):
	fig,ax=plt.subplots()
	func = np.sin(theta1_new)*np.sin(theta2_new)*np.cos(phi)-2.0*np.cos(theta1_new)*np.cos(theta2_new)

	# Data file generation for the contour plot.
	theta1_2d, theta2_2d = np.meshgrid(theta1, theta2)

	# Plot env is started here.
	cp = ax.contourf(theta1_2d, theta2_2d, func)
	colorbar_format = '% 2.1f'
	fig.colorbar(cp, pad=0.02, format=colorbar_format)

	# Set axes labels
	ax.set_xlabel(r'$\theta_1$', labelpad=10)
	ax.set_ylabel(r'$\theta_2$', labelpad=10)

	# Set ticklabels
	tick_position = [0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]
	tick_label = [r"$0$", r"$\dfrac{1}{2}\pi$", r"$\pi$", r"$\dfrac{3}{2}\pi$", r"$2\pi$"]

	ax.set_xticks(tick_position)
	ax.set_xticklabels(tick_label)
	ax.set_yticks(tick_position)
	ax.set_yticklabels(tick_label)
	 
	if (phi == 0):
		title_name = r'$\phi=0$'
		str_name = "0"
	elif (phi == np.pi/2):
		title_name = r'$\phi=90$'
		str_name = "90"
	elif (phi == np.pi):
		title_name = r'$\phi=180$'
		str_name = "180"

	plt.rcParams['axes.titley'] = 1.0	
	plt.rcParams['axes.titlepad'] = -24
	# displaying the title 
	plt.title(title_name, color='blue')
	#raw_output_file=home+"/academic-project/outputs/final-pigs-outputs-for-plotting/Figure-dipole-dipole-interaction-potential-for-phi-" + str(angle) + "degree.pdf"

# Adjustment the plot
plt.subplots_adjust(top=0.98,bottom=0.2,left=0.13,right=1,hspace=0.0,wspace=0.0)

# Saving the figure
output_file=home+"/academic-project/output/final-pigs-outputs-for-plotting/Figure-dipole-dipole-interaction-potential.pdf"
print(output_file)
plt.savefig(output_file, format='pdf')
plt.show()
