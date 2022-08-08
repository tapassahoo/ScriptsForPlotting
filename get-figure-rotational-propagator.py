import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Reading the file, we would like to plot.
inputFile="/home/tapassahoo/academic-project/outputs/rotational-propagator-for-water/rho.den018_1500K_spinless"
inputFile1="/home/tapassahoo/academic-project/outputs/rotational-propagator-for-water/rho.den018"
rho = np.genfromtxt(inputFile, unpack=True, usecols=[3])
rho1 = np.genfromtxt(inputFile1, unpack=True, usecols=[3])
print(max(rho))
print(max(rho1))
exit()

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
