#!/usr/bin/env python
""" 
This code is exploited for a 3-dimensional surface plot. 
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

# Import modules
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os

# Reading the file, we would like to plot.
home=os.path.expanduser("~")
inputFile=home+"/academic-project/outputs/rotational-propagator-for-water/rho.den018_ortho"
print(inputFile)
exit()
rho = np.genfromtxt(inputFile, unpack=True, usecols=[3])

size_grid = 361
theta=np.arange(361)
phi=theta
theta2d, phi2d = np.meshgrid(theta, phi)
rho2d=np.reshape(rho, (size_grid, size_grid))

# Plot env is started here
fig = plt.figure(figsize = (12,10))
ax = plt.axes(projection='3d')
surf = ax.plot_surface(theta2d, phi2d, rho2d, cmap = plt.cm.cividis)

# Set axes label
ax.set_xlabel('x', labelpad=20)
ax.set_ylabel('y', labelpad=20)
ax.set_zlabel('z', labelpad=20)

fig.colorbar(surf, shrink=0.5, aspect=8)
plt.show()
