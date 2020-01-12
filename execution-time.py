from numpy import *
import numpy as np
import matplotlib
#matplotlib.use('eps')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from subprocess import call
#from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import ScalarFormatter
import mypkg.pkgMoribs.support_without_parallel as support
import os
from pylab import *
import matplotlib.axes
##rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
from collections import OrderedDict

nthreads, exetime = genfromtxt("execution-time-N10-P16.txt", unpack=True, usecols=[0, 1], skip_header=0, skip_footer=0)
# create some data to use for the plot
print(nthreads)
print(exetime)

fig = plt.figure(figsize=(8, 6))
#plt.grid(True)
myfont = 20
fontlegend = myfont/2.0

# the main axes is subplot(111) by default
plt.title('PIGS: N=10, P=17, Blocks=1000', fontsize=myfont)
plt.plot(nthreads, exetime, marker='o',markersize=10)
#plt.axis([0, 1, 1.1*np.amin(s), 2*np.amax(s)])
plt.xlabel('Number of threads', fontsize=myfont, labelpad=10)
plt.ylabel('Execution time in seconds', fontsize=myfont, labelpad=10)
plt.xticks(fontsize=myfont, rotation=0)
plt.yticks(fontsize=myfont, rotation=0)

plt.subplots_adjust(top=0.92, bottom=0.12, left=0.12, right=0.99, hspace=0.0, wspace=0.0)
plt.savefig("execution-time-vs-nthreads.pdf", dpi = 200, format = 'pdf')

plt.show()
