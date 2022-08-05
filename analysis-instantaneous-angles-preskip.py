import os
import sys
from subprocess import call
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes
#matplotlib.use('eps')
from pylab import *
import mypkg.pkgMoribs.support_without_parallel as support

rc('text', usetex=True)
size=24
params = {'legend.fontsize': size*0.5,
	'figure.figsize': (8,18),
	'axes.labelsize': size,
	'axes.titlesize': size,
	'xtick.labelsize': size*0.75,
	'ytick.labelsize': size*0.75,
	'axes.titlepad': size}
plt.rcParams.update(params)
matplotlib.rcParams.update({'font.size': size*0.75})
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"


###  Program starts from here  ########

#Read some parameters from the command-line
dir_read = sys.argv[1]
mc_read = sys.argv[2]
numb_particle = int(sys.argv[3])
dofs_read1 = int(sys.argv[4])
axis_read= sys.argv[5]
numb_beads= int(sys.argv[6])
particle_index = int(sys.argv[7])

#
#particle_index = int(sys.argv[2])
axis_index = {"cost":0, "phi":1, "chi":2}
dofs_read = {0:"TransAndRot", 1:"Rot", 2:"Trans"} 

#
file_read = dir_read + 'results/output.xyz'
string1 =dir_read[dir_read.find(mc_read):-1]
numb_blocks = int(string1[string1.find("Blocks")+6:string1.find("-Passes")])
system_replaced = string1[string1.find("System"):string1.find("p-H2O")+5]
rpt_read=string1[string1.find("Rpt"):string1.find("Angstrom")]

first_fragment = dir_read[:dir_read.find("Rpt")]
last_fragment = dir_read[dir_read.find("Angstrom")+8:]

ndofs=3
if (mc_read != "PIMC"):
	beads_pos = int((numb_beads-1)/2)

fig, ax = plt.subplots()
isubplot = 1
system_replaced_by = {2:"System2-p-H2O", 11:"System11-p-H2O"}
label_panel_list = {0:"(a)", 1:"(b)", 2:"(c)"}
rptdict = {0:[3.0, 3.2], 1:[3.4, 3.6], 2:[4.0, 5.0]}

for isubplot in range(3):
	label_panel = label_panel_list[isubplot]
	plt.subplot(3, 1, isubplot+1)
	preskip = 0
	postskip = 0

	rptList=rptdict[isubplot]
	print(rptList)
	colorList=["yellow", "green", "red", "magenta", "red", "cyan"]
	j=axis_index[axis_read]
	ig=0
	for rpt in rptList:
		file1=file_read.replace(system_replaced, system_replaced_by[numb_particle])
		file2=file1.replace(rpt_read,"Rpt"+str(rpt))
		dir_rpt=dir_read.replace(rpt_read,"Rpt"+str(rpt))
		final_dir_in_work = dir_rpt

		file_old = final_dir_in_work+"results/output.xyz_old"
		if (os.path.isfile(file_old) == True):
			file_new = final_dir_in_work+"results/output.xyz"

			if (os.path.isfile(file_new) == True):
				print(" -- Restarted data")

				if "H2O1" in open(file_new).read():
					rmstr = int(numb_particle*numb_beads+3)
					cmd1="tail -n +"+str(rmstr)+" "+final_dir_in_work+"results/output.xyz>bb"
					os.system(cmd1)
					col_data_new = np.genfromtxt("bb")
					call(["rm", "bb"])
				else:
					col_data_new = np.genfromtxt(final_dir_in_work+"results/output.xyz")
				index = int(col_data_new[0,0])
				col_data_old = np.genfromtxt(final_dir_in_work+"results/output.xyz_old")
				marged_data  = np.concatenate((col_data_old[:index-1], col_data_new), axis=0)
				aa = col_data_new[:,0]
				final_data_set = marged_data[0:int(aa[-1]),:]
			else:
				final_data_set = np.genfromtxt(final_dir_in_work+"results/output.xyz_old", skip_header=0, skip_footer=0)
		else:
			final_data_set = np.genfromtxt(final_dir_in_work+"results/output.xyz", skip_header=0, skip_footer=0)

		data_len = len(final_data_set[:,0])
		nlen = int(len(final_data_set[:,0]))
	
		workingNdim = int(math.log(data_len)/math.log(2))
		trunc = int(data_len-2**workingNdim)
		data_len=int(data_len-(preskip+trunc))
		print(data_len)
	
		save_data = np.zeros((numb_particle,data_len))

		print(file2)
		for i in range(numb_particle):
			ncol1 = beads_pos+i*numb_beads
			ncol = j+ncol1*ndofs
			ncol = ncol+1
			print(str(ncol)+'th column')
			save_data[i,:] = final_data_set[(preskip+trunc):(nlen-postskip),ncol] 
			

		label_str=r'$r$='+str(rpt)	

		if (numb_particle == 11):
			vec_plot = np.reshape(save_data[2:(numb_particle-2),:], (numb_particle-4)*data_len)
		if (numb_particle == 2):
			vec_plot = np.reshape(save_data, numb_particle*data_len)
		data_plot = vec_plot
		plt.hist(data_plot, bins=50, density=True, alpha=0.5, stacked=True, color=colorList[ig], edgecolor='black', label=label_str)
		ig=ig+1

	#
	numb_label=int(numb_particle)
	plt.ylabel(r'$\mathrm{Density}$',labelpad=5)
	#plt.xlim(-1.01,1.01)
	xmin,xmax=plt.xlim()
	ymin,ymax=plt.ylim()
	plt.text(xmin+(xmax-xmin)*0.04,ymax-(ymax-ymin)*0.1,label_panel)
	plt.text(xmin+0.48*(xmax-xmin),ymin+0.9*(ymax-ymin),r'$N$='+str(numb_label))
	isubplot = isubplot+1
	if (isubplot == 3):
		plt.xlabel(r'$\mathrm{bins \ of \ \cos(\theta)}$', labelpad=5)
	if (isubplot != 3):
		frame1 = plt.gca()
		frame1.axes.xaxis.set_ticklabels([])
	if (numb_particle == 2):
		plt.subplots_adjust(top=0.99,bottom=0.05,left=0.11,right=0.98,hspace=0.06,wspace=0.0)
	if (numb_particle == 11):
		plt.subplots_adjust(top=0.99,bottom=0.05,left=0.10,right=0.98,hspace=0.06,wspace=0.0)
	ax.tick_params(right=True,top=True,left=True,bottom=True)
	plt.legend(numpoints=1,loc='center')
	plt.minorticks_on()
	plt.tick_params(axis="both", direction="out", which="minor", right=True, top=True, length=2)
	plt.tick_params(axis="both", direction="out", which="major", right=True, top=True, length=5)

#
index_cut=dir_read.find(mc_read)
home = os.path.expanduser("~")
final_results_path = home + "/ResultsOf" + mc_read + "/"
FilePlotDensity=final_results_path+first_fragment[index_cut:-1]+last_fragment[:-1]+"-histogram-of-"+axis_read+"-vs-gFactor-preskip"+str(preskip)+"-postskip"+str(postskip)+".pdf"
print(FilePlotDensity)

plt.savefig(FilePlotDensity, dpi=100, format='pdf')
plt.show()

#python analysis_instantaneous_angles_preskip.py /Users/tsahoo/nonlinear-rotors/PIGS-qTIP4P-RotDOFs-Rpt3.0Angstrom-beta0.1Kinv-Blocks20000-Passes200-System2-p-H2O-e0vsbeads51 PIGS 2 1 cost 51 0 

