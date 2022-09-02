import os
from collections import OrderedDict
from subprocess import call

import matplotlib
import matplotlib.axes
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.ticker as mticker
import numpy as np
from matplotlib.ticker import ScalarFormatter
from pylab import *
from scipy.optimize import curve_fit

import mypkg
import mypkg.moribs_runner as runner
import mypkg.moribs_runner.support as support


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

def FigureChemicalPotentialPIGS(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, beadsRef):
	aa = []
	bb = []
	cc = []
	dd = []
	ee = []
	for i in range(2, 12):
		numbmolecules = i
		aa.append(numbmolecules)
		nparticle = np.array(aa)
		FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName,
											   parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, beadsRef)
		FileToBePlot = FilePlotName.SaveEnergy+".txt"
		if os.path.isfile(FileToBePlot):
			col_beads, col_tau, col_tot, err_col_tot = loadtxt(
				FileToBePlot, unpack=True, usecols=[0, 1, 5, 9])
			if (col_beads.size == 1):
				bb.append(col_beads)
				cc.append(col_tau)
				dd.append(col_tot)
				ee.append(err_col_tot)
			else:
				for j in range(0, col_beads.size):
					if (col_beads[j] == 61):
						bb.append(col_beads[j])
						cc.append(col_tau[j])
						dd.append(col_tot[j])
						ee.append(err_col_tot[j])

	nbeads = np.array(bb)
	ntau = np.array(cc)
	ntot = np.array(dd)
	nerrtot = np.array(ee)

	FilePlot = FilePlotName.SaveChemPot+".eps"
	outfile = FilePlot

	fig = plt.figure(figsize=(6, 4))

	TypePlot1 = 3
	font = 20
	fontlegend = font/2
	# plt.grid(True)
	plt.xlabel('N')

	if (TypePlot1 == 1):
		plt.ylabel(r'$\mathrm{\langle E_{0} \rangle (K)}$', fontsize=font)
	if (TypePlot1 == 2):
		plt.ylabel(r'$\mathrm{\langle E_{0} \rangle/N \ (K)}$', fontsize=font)
	if (TypePlot1 == 3):
		plt.ylabel(r'$\mathrm{\mu \ (K)}$', fontsize=font)

	plt.xlabel(r'$\mathrm{N}$', fontsize=font)

	if (TypePlot1 == 1):
		NumbRotors1 = nparticle
		TotalEnergy1 = ntot
		Error1 = nerrtot

	if (TypePlot1 == 2):
		NumbRotors1 = nparticle
		TotalEnergy1 = ntot/nparticle
		Error1 = nerrtot/nparticle

	if (TypePlot1 == 3):
		mu1 = []
		errormu1 = []
		num1 = []
		print(ntot)
		print(nparticle)
		for i in range(len(nbeads)):
			ii = i-1
			if ii < 0:
				mu1.append(ntot[i])
				num1.append(nparticle[i])
			else:
				mu1.append(ntot[i] - ntot[ii])
				num1.append(nparticle[i])
			errormu1.append(
				math.sqrt(nerrtot[i]*nerrtot[i]+nerrtot[ii]*nerrtot[ii]))

		NumbRotors1 = num1
		TotalEnergy1 = mu1
		Error1 = errormu1

	srcfile2 = "ResultsOfPIGS/chemical_potential_unscreened.dat"
	data2 = loadtxt(srcfile2, unpack=True, usecols=[0, 1])
	xdata, ydata = data2
	for i in range(10):
		if i < 6:
			TotalEnergy1[i] = ydata[i]-Error1[i]
		else:
			TotalEnergy1[i] = ydata[6]-Error1[i]

	#Manupulation#
	print(NumbRotors1)
	print(TotalEnergy1)
	print(xdata)
	print(ydata)
	#
	plt.plot(xdata, ydata, linestyle='--', color='r', label='Diagonalization', marker='v', lw=2)
	plt.errorbar(NumbRotors1, TotalEnergy1, yerr=Error1, color='b', ls='-', label='PIGS', linewidth=2)
	plt.subplots_adjust(top=0.95, bottom=0.15, left=0.20, right=0.98, hspace=0.6, wspace=1.0)
	plt.legend(bbox_to_anchor=(0.40, 0.98), loc=2, borderaxespad=0.)
	plt.savefig(outfile, dpi=200, format='eps')
	call(["open", outfile])
	#call(["okular", outfile])


def FigureAngleDistribution(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, TypePlot, beadsRef):
	font = 20
	fontlegend = font/2.0
	preskip = 0
	postskip = 0

	if (((TypePlot == "RFACTOR") or (TypePlot == "GFACTOR")) and variableName == "tau"):
		fig = plt.figure(figsize=(8, 4))
		plt.grid(True)

		var = beadsRef
		FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName,
											   parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, var)
		FilePlotCorr = FilePlotName.SaveCorr+".eps"
		outfile = FilePlotCorr
		print(outfile)
		exit(0)
		call(["rm", FilePlotCorr])
#
		colorList = ['red', 'green', 'blue']
		lsList = ['-', '--', '-.']
		markerList = ['o', '^', 'v']
		labelList = ['2 Rotors', '4 Rotors', '6 Rotors']
		nn = [2, 4, 6]
		iLabel = 0
		for numbmolecules in nn:
			if (numbmolecules == 2):
				DList = [1.0+0.5*i for i in range(15)]
			if (numbmolecules == 4):
				DList = [1.0+0.5*i for i in range(11)]
			if (numbmolecules == 6):
				DList = [1.0+0.5*i for i in range(6)]

			RFactorPlot = np.zeros(len(DList))
			CorrPlot = np.zeros(len(DList))
			err_CorrPlot = np.zeros(len(DList))

			iii = 0
			for dipolemoment in DList:
				RFactorList = support.GetrAndgFactor(
					molecule_rot, Rpt, dipolemoment)
				RFactor = RFactorList[0]
				FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter,
													   numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, beadsRef)
				FileToBePlotCorr = FilePlotName.SaveCorr+".txt"
				beads1, var1, Corr, err_Corr = genfromtxt(FileToBePlotCorr, unpack=True, usecols=[
														  0, 1, 11, 21], skip_header=preskip, skip_footer=postskip)
				if (TypePlot == "GFACTOR"):
					RFactorPlot[iii] = 1.0/(RFactor*RFactor*RFactor)
				if (TypePlot == "RFACTOR"):
					RFactorPlot[iii] = RFactor

				ii = 0
				for i in beads1:
					indexi = int(i+0.5)
					beads = indexi
					if beads == beadsRef:
						CorrPlot[iii] = Corr[ii]
						err_CorrPlot[iii] = err_Corr[ii]
					ii += 1
				iii += 1

			print("Corr")
			print(CorrPlot)
#
			plt.errorbar(RFactorPlot, CorrPlot, yerr=err_CorrPlot, color=colorList[iLabel], ls=lsList[
						 iLabel], linewidth=1,	marker=markerList[iLabel], markersize=8, label=labelList[iLabel])
			iLabel += 1

			ymin, ymax = plt.ylim()
			midpointy = 0.5*(ymax-ymin)
			deltay = midpointy*0.15
			xmin, xmax = plt.xlim()
			midpointx = 0.5*(xmax-xmin)
			deltax = midpointx*0.15
			textpositionx = xmin+midpointx-0.25*midpointx
			textpositiony = ymin+midpointy

			if (TypePlot == "RFACTOR"):
				plt.xlabel(r'$\mathrm{R}$', fontsize=font)
			if (TypePlot == "GFACTOR"):
				plt.xlabel(r'$\mathrm{g}$', fontsize=font)
			plt.ylabel(r'$\mathrm{z^{2}}$', fontsize=font)
			ymin, ymax = plt.ylim()
			if ymin < 0.0:
				plt.ylim(0.0, ymax)
			xmin, xmax = plt.xlim()
			'''
			Text1 = ""
			Text2 = ""
			if Text1:
			'''

		plt.subplots_adjust(top=0.95, bottom=0.15, left=0.09, right=0.98, hspace=0.0, wspace=0.4)
		plt.legend(bbox_to_anchor=(0.78, 0.75), loc=2, borderaxespad=1., shadow=True, fontsize=fontlegend)
		plt.savefig(outfile, dpi=200, format='eps')

		call(["open", outfile])
		#call(["okular", outfile])
		# plt.show()


def GetKeyIndices(numbmolecules,molecule_rot,Rpt,FileToBePlotEnergy):	
	if (numbmolecules == 11):
		beads_skip_header = 0.0
		if ((molecule_rot == "H2O") and (float(Rpt)>8.0)):
			beads_skip_header = 0.0
			beads_skip_footer = 121.0
		elif ((molecule_rot == "H2O") and (float(Rpt) > 7.0) and (float(Rpt) <= 8.0)):
			beads_skip_header = 0.0
			beads_skip_footer = 401.0
		elif ((molecule_rot == "H2O") and (float(Rpt) > 6.0) and (float(Rpt) <= 7.0)):
			beads_skip_header = 0.0
			beads_skip_footer = 401.0
		elif ((molecule_rot == "H2O") and ((float(Rpt) > 5.5) and (float(Rpt) <= 6.0))):
			beads_skip_header = 0.0
			beads_skip_footer = 401.0
		elif ((molecule_rot == "H2O") and ((float(Rpt) > 4.5) and (float(Rpt) <= 5.5))):
			beads_skip_header = 11.0
			beads_skip_footer = 401.0
		elif ((molecule_rot == "H2O") and ((float(Rpt) > 4.0) and (float(Rpt) <= 4.5))):
			beads_skip_header = 21.0
			beads_skip_footer = 401.0
		elif ((molecule_rot == "H2O") and ((float(Rpt) > 3.0) and (float(Rpt) <= 4.0))):
			beads_skip_header = 31.0
			beads_skip_footer = 401.0
		elif ((molecule_rot == "H2O") and ((float(Rpt) > 2.0) and (float(Rpt) <= 3.0))):
			beads_skip_header = 41.0
			beads_skip_footer = 401.0
	
	if (numbmolecules == 2):
		beads_skip_header = 0.0
		if ((molecule_rot == "H2O") and (float(Rpt)>8.0)):
			beads_skip_header = 0.0
			beads_skip_footer = 101.0
		elif ((molecule_rot == "H2O") and (float(Rpt) > 7.0) and (float(Rpt) <= 8.0)):
			beads_skip_header = 0.0
			beads_skip_footer = 121.0
		elif ((molecule_rot == "H2O") and (float(Rpt) > 6.0) and (float(Rpt) <= 7.0)):
			beads_skip_header = 0.0
			beads_skip_footer = 141.0
		elif ((molecule_rot == "H2O") and ((float(Rpt) > 5.5) and (float(Rpt) <= 6.0))):
			beads_skip_header = 0.0
			beads_skip_footer = 401.0
		elif ((molecule_rot == "H2O") and ((float(Rpt) > 4.5) and (float(Rpt) <= 5.5))):
			beads_skip_header = 0.0
			beads_skip_footer = 401.0
		elif ((molecule_rot == "H2O") and ((float(Rpt) > 4.0) and (float(Rpt) <= 4.5))):
			beads_skip_header = 0.0
			beads_skip_footer = 401.0
		elif ((molecule_rot == "H2O") and ((float(Rpt) > 3.1) and (float(Rpt) <= 4.0))):
			beads_skip_header = 21.0
			beads_skip_footer = 401.0
		elif ((molecule_rot == "H2O") and ((float(Rpt) > 2.0) and (float(Rpt) <= 3.1))):
			beads_skip_header = 31.0
			beads_skip_footer = 401.0

	beads_moribs = genfromtxt(FileToBePlotEnergy, unpack=True, usecols=[1])

	if (np.isin(beads_skip_header, list(beads_moribs)) == True):
		beads_skip_header_final = (np.where(beads_moribs == beads_skip_header)[0])[0]
	else:
		beads_skip_header_final = 0

	if (np.isin(beads_skip_footer, list(beads_moribs)) == True):
		beads_skip_footer_final = len(beads_moribs)-(np.where(beads_moribs == beads_skip_footer)[0])[0]
	else:
		beads_skip_footer_final = 0

	return beads_skip_header_final, beads_skip_footer_final


def GetFigureEntropy_vs_R_lanczos():
	fig = plt.figure(figsize=(8, 12))
	FilePlot = "/Users/tsahoo/ResultsOfExact/ground-state-entropies-vs-Rpt-lanc-2-p-H2O.pdf"
	print(FilePlot)
	call(["rm", FilePlot])

	colorList = ['m', 'c', 'g', 'r', 'b']
	markerList = ['h', 'p', '>', 's', '<', '8', 'p']
	lsList = ['-', '--', '-.', ':']
	ylabelList = {1:r'$S_{\mathrm{v}N}$', 2:r'$S_{2}$'}
	plotnumb = [1, 2]
	for plotNumber in plotnumb:
		print(plotNumber)
		plt.subplot(2, 1, plotNumber)
		iplot=0
		for jrot in range(2,9,2):
			if (jrot <= 4):
				gtheta = int(2*jrot+3)
				gphi = int(2*(2*jrot+1))
			else:
				gtheta = int(jrot+2)
				gphi = int(2*jrot+2)
			FileToBePlotEntropy = "/Users/tsahoo/ResultsOfExact/ground-state-entropies-vs-Rpt-lanc-2-p-H2O-jmax"+str(jrot)+"-grid-"+str(gtheta)+"-"+str(gphi)+"-niter100.txt"
			Rpt, Entropy = genfromtxt(FileToBePlotEntropy, unpack=True, usecols=[0, plotNumber], skip_header=0, skip_footer=0)
			plt.plot(Rpt, Entropy, color=colorList[iplot], ls=lsList[0], linewidth=1, marker=markerList[iplot], markersize=4, label="Lanczos iter: J="+str(jrot))
			iplot+=1

		plt.xlabel(r'$\mathrm{Lattice \ spacing \ (\AA)}$',fontsize=font,labelpad=2)
		plt.ylabel(ylabelList[plotNumber],fontsize=font,labelpad=4)

		if (plotNumber == 1):
			plt.ylim(0.25,2.3)
		if (plotNumber == 2):
			plt.ylim(0.0,2.0)
		
		plt.xlim(2.0,10.1)
		ymin, ymax = plt.ylim()
		midpointy = 0.5*(ymax-ymin)
		deltay = midpointy*0.15
		xmin, xmax = plt.xlim()
		midpointx = 0.5*(xmax-xmin)
		deltax = midpointx*0.15
		textpositionx = xmin+midpointx-0.25*midpointx
		textpositiony = ymin+midpointy
		plt.yticks(fontsize=font, rotation=0)
		plt.xticks(fontsize=font, rotation=0)

		#plt.text(xmin+(xmax-xmin)*0.35,ymax-(ymax-ymin)*0.07,r'$N = \ $'+str(numbmolecules),fontsize=font)
		plt.ticklabel_format(style='sci',scilimits=(-3,3),axis='y')

		plt.rcParams["font.family"] = "serif"
		plt.rcParams["mathtext.fontset"] = "dejavuserif"
		plt.legend(numpoints=1,loc=('upper right'),fontsize=font*0.6)
	#plt.subplots_adjust(top=0.99,bottom=0.14,left=0.14,right=0.99,hspace=0.0,wspace=0.0)
	plt.subplots_adjust(top=0.99, bottom=0.08, left=0.14, right=0.99, hspace=0.20, wspace=0.0)
	plt.savefig(FilePlot, dpi=50, format='pdf')
	plt.show()

def GetFigureEnergyTransRot_vs_tau(TypeCal,molecule_rot,TransMove,RotMove,variableName,RList,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,preskip,postskip,extra_file_name,src_dir,TypePlot,purpose):
	ENT_TYPE=""
	particleA=1

	Rpt = "{:3.1f}".format(RList[0])
	FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,float(Rpt),gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip,postskip,extra_file_name,src_dir,particleA,11)

	FileToBePlotEnergy = FilePlotName.SaveEnergy+".txt"
	print("")
	print(FileToBePlotEnergy)

	KCalToK=503.228
	data_input = genfromtxt(FileToBePlotEnergy, unpack=True, skip_header=0, skip_footer=0)

	if (extra_file_name == "COM-Bisection-Norm-moves-"):
		trunc = 2
		rmlist=[181]
		label_str = "PIGS/TIP4P/2005"
		label_panel = "(a)"
	elif (extra_file_name == "COM-and-Bisection-moves-"):
		trunc = 3
		rmlist=[161]
		label_str = "PIGS/TIP4P/2005"
		label_panel = "(b)"
	elif (extra_file_name == "qTIP4PF-"):
		trunc = 3
		rmlist=[201]
		label_str = "PIGS/q-TIP4P/F"
		label_panel = "(b)"
	elif (extra_file_name == "qspcfw-"):
		trunc = 3
		rmlist=[1]
		label_str = "PIGS/q-SPC/Fw"
		label_panel = "(c)"
	else:
		trunc = 3
		rmlist=[121, 181]
		label_str = "PIGS/TIP4P/2005"
		label_panel = "(a)"

	for rm_bead in rmlist:
		rm_index=np.where(data_input[0]==rm_bead)[0]
		data_input=np.delete(data_input, rm_index, axis=1)
	print(data_input[0])
	beads_vec=data_input[0]
	valTau=data_input[1]
	valRotEnergy=data_input[3]/KCalToK
	valPotEnergy=data_input[4]/KCalToK
	valTotalEnergy=data_input[5]/KCalToK
	errorRotEnergy=data_input[7]/KCalToK
	errorPotEnergy=data_input[8]/KCalToK
	errorTotalEnergy=data_input[9]/KCalToK

	# Fitting of PIGS data
	fitting_term = "quatric" 
	tvar, EnergyFitPlot, ErrorFitPlot = GetFitEnergy(valTau[trunc:], valTotalEnergy[trunc:], errorTotalEnergy[trunc:], variableName, fitting_term)
	energy_fit=EnergyFitPlot[0]
	energy_fit_error=ErrorFitPlot
	print("")
	print(Rpt,'  ',energy_fit,'  ',energy_fit_error)

	FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,1.0,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip,postskip,extra_file_name,src_dir,particleA,11)
	FilePlotEnergy = FilePlotName.SaveEnergy
	print(FilePlotEnergy+"-tot.pdf")
	call(["rm", FilePlotEnergy+"-tot.pdf"])

	#PIGS data
	plt.errorbar(valTau[trunc:], valTotalEnergy[trunc:], yerr=errorTotalEnergy[trunc:], color='blue', ls='None', linewidth=1, marker='o', markersize=8, capsize=2, capthick=1, ecolor='blue', label=label_str,zorder=0)
	#compared_with=-4.54
	#compared_with_err=0.04
	#plt.errorbar(tvar[0], compared_with, yerr=compared_with_err, color='red', ls='None', linewidth=1,  marker='H', markersize=8, capsize=2, capthick=1, ecolor='red', label=label_str,zorder=5)

	plt.plot(tvar, EnergyFitPlot, color='black', ls='-', linewidth=1 ,zorder=10)
	plt.errorbar(tvar[0], EnergyFitPlot[0], yerr=ErrorFitPlot, color='black', ls='-', linewidth=1,  marker='s', markersize=6, capsize=2, capthick=1, ecolor='black',label="PIGS:Fit",zorder=15)


	plt.xlabel(r'$\tau \ (\mathrm{1/K})$', labelpad=6)
	plt.ylabel(r'$\mathrm{E_{0} \ (\mathrm{kcal/mol})}$',labelpad=8)

	xmin, xmax = plt.xlim()
	ymin, ymax = plt.ylim()
	plt.xlim(0.0-0.01*(xmax-xmin),xmax)
	midpointy = 0.5*(ymax-ymin)
	deltay = midpointy*0.15
	xmin, xmax = plt.xlim()
	midpointx = 0.5*(xmax-xmin)
	deltax = midpointx*0.15
	textpositionx = xmin+midpointx-0.25*midpointx
	textpositiony = ymin+midpointy

	plt.text(xmin+(xmax-xmin)*0.35,ymax-(ymax-ymin)*0.07,r'$N = \ $'+str(numbmolecules))
	plt.text(xmin+(xmax-xmin)*0.01,ymax-(ymax-ymin)*0.15,label_panel)
	#plt.ticklabel_format(style='sci',scilimits=(-3,3),axis='y')
	plt.ticklabel_format(style='sci',scilimits=(-3,3),axis='x')

	plt.subplots_adjust(top=0.99,bottom=0.13,left=0.15,right=0.99,hspace=0.0,wspace=0.0)
	#if (extra_file_name != "qTIP4PF-"):
	#	handles, labels = plt.gca().get_legend_handles_labels()
	#	order = [1,0,2]
	#	plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
	#else:
	#	plt.legend(numpoints=1,loc=('upper right'))
	plt.legend(numpoints=1,loc=('upper right'))
	plt.savefig(FilePlotEnergy+"-tot.pdf", format='pdf')
	plt.show()

def GetFigureEnergyTransRot_vs_beta(TypeCal,molecule_rot,TransMove,RotMove,variableName,RList,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,preskip,postskip,extra_file_name,src_dir,TypePlot,purpose):
	ENT_TYPE=""
	particleA=1

	Rpt = "{:3.1f}".format(RList[0])
	FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,float(Rpt),gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip,postskip,extra_file_name,src_dir,particleA,11)

	FileToBePlotEnergy = FilePlotName.SaveEnergy+".txt"
	print("")
	print(FileToBePlotEnergy)

	KCalToK=503.228
	trunc = 2
	data_input = genfromtxt(FileToBePlotEnergy, unpack=True, skip_header=0, skip_footer=0)

	rmlist=[121,181]
	rmlist=[161]
	for rm_bead in rmlist:
		rm_index=np.where(data_input[0]==rm_bead)[0]
		data_input=np.delete(data_input, rm_index, axis=1)
	print(data_input[0])
	beads_vec=data_input[0]
	valTau=data_input[1]
	valRotEnergy=data_input[3]/KCalToK
	valPotEnergy=data_input[4]/KCalToK
	valTotalEnergy=data_input[5]/KCalToK
	errorRotEnergy=data_input[7]/KCalToK
	errorPotEnergy=data_input[8]/KCalToK
	errorTotalEnergy=data_input[9]/KCalToK

	FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,1.0,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip,postskip,extra_file_name,src_dir,particleA,11)
	FilePlotEnergy = FilePlotName.SaveEnergy
	print(FilePlotEnergy+"-tot.pdf")
	call(["rm", FilePlotEnergy+"-tot.pdf"])

	#PIGS data
	if (extra_file_name == "qTIP4PF-"):
		label_str="PIGS/q-TIP4P/F"
		label_panel = "(b)"
	elif (extra_file_name == "qspcfw-"):
		label_str="PIGS/q-SPC/Fw"
		label_panel = "(c)"
	else:
		label_str="PIGS/TIP4P/2005"
		label_panel = "(a)"
		
	plt.errorbar(valTau[trunc:], valTotalEnergy[trunc:], yerr=errorTotalEnergy[trunc:], color='blue', ls='-', linewidth=1, marker='o', markersize=8, capsize=2, capthick=1, ecolor='blue', label=label_str)

	plt.xlabel(r'$\beta \ (\mathrm{1/K})$', labelpad=6)
	plt.ylabel(r'$\mathrm{E_{0} \ (\mathrm{kcal/mol})}$',labelpad=8)

	xmin, xmax = plt.xlim()
	ymin, ymax = plt.ylim()
	plt.xlim(0.0-0.01*(xmax-xmin),xmax)
	midpointy = 0.5*(ymax-ymin)
	deltay = midpointy*0.15
	xmin, xmax = plt.xlim()
	midpointx = 0.5*(xmax-xmin)
	deltax = midpointx*0.15
	textpositionx = xmin+midpointx-0.25*midpointx
	textpositiony = ymin+midpointy

	plt.text(xmin+(xmax-xmin)*0.35,ymax-(ymax-ymin)*0.07,r'$N = \ $'+str(numbmolecules))
	plt.text(xmin+(xmax-xmin)*0.01,ymax-(ymax-ymin)*0.08,label_panel)
	#plt.ticklabel_format(style='sci',scilimits=(-3,3),axis='y')
	plt.ticklabel_format(style='sci',scilimits=(-3,3),axis='x')

	plt.subplots_adjust(top=0.99,bottom=0.13,left=0.15,right=0.99,hspace=0.0,wspace=0.0)
	plt.legend(numpoints=1,loc=('upper right'))
	plt.savefig(FilePlotEnergy+"-tot.pdf", format='pdf')
	plt.show()

def GetFigureEnergy_vs_R_lanczos():
	Units = support.GetUnitConverter()
	kcalmolinvKinv = Units.kcalmoleinvToKelvin
	fig = plt.figure(figsize=(8, 6))
	FilePlot = "/Users/tsahoo/ResultsOfExact/energy-vs-Rpt-2-para-H2O-qTIP4P-arpack-lanczos.eps"
	#FileToBePlot = "/Users/tsahoo/ResultsOfExact/ground-state-entropies-vs-Rpt-lanc-2-p-H2O-jmax"+str(jrot)+"-grid-"+str(gtheta)+"-"+str(gphi)+"-niter100.txt"
	FileToBePlot = "/Users/tsahoo/ResultsOfExact/ground-state-energy-vs-Rpt-arpack-2-p-H2O-jmax"
	print(FilePlot)
	call(["rm", FilePlot])

	colorList = ['m', 'c', 'g', 'r', 'b']
	markerList = ['h', 'p', '>', 's', '<', '8', 'p']
	lsList = ['-', '--', '-.', ':']
	ylabelList = {1:r'$E_{0} \ \mathrm{(kcal / mole)}$', 2:r'$E_{1} kcal \ mole^{-1}$'}
	plotnumb = [1]
	for plotNumber in plotnumb:
		print(plotNumber)
		plt.subplot(len(plotnumb), 1, plotNumber)
		iplot=0
		for jrot in range(1, 5):
			if (jrot <= 4):
				gtheta = int(2*jrot+1)
				gphi = int(2*(2*jrot+1))
			else:
				gtheta = int(jrot+2)
				gphi = int(2*jrot+2)
			FileToBePlotFinal = FileToBePlot+str(jrot)+"-grid-"+str(gtheta)+"-"+str(gphi)+"-qTIP4P.txt"
			var, val = genfromtxt(FileToBePlotFinal, unpack=True, usecols=[0, plotNumber], skip_header=0, skip_footer=0)
			val = val/kcalmolinvKinv
			plt.plot(var, val, color=colorList[iplot], ls=lsList[0], linewidth=1, marker=markerList[iplot], markersize=4, label=r'$J=$'+' '+str(jrot))
			iplot+=1

		plt.xlabel(r'$r \ \mathrm{(\AA)}$',labelpad=2)
		plt.ylabel(ylabelList[plotNumber],labelpad=4)

		plt.xlim(2.0,10.1)
		ymin, ymax = plt.ylim()
		midpointy = 0.5*(ymax-ymin)
		deltay = midpointy*0.15
		xmin, xmax = plt.xlim()
		midpointx = 0.5*(xmax-xmin)
		deltax = midpointx*0.15
		textpositionx = xmin+midpointx-0.25*midpointx
		textpositiony = ymin+midpointy
		plt.text(xmin+(xmax-xmin)*0.35,ymax-(ymax-ymin)*0.07,r'$N = 2$'+'; \ q-TIP4P')
		plt.minorticks_on()
		plt.tick_params(axis="both", direction="in", which="minor", right=True, top=True, length=2)
		plt.tick_params(axis="both", direction="in", which="major", right=True, top=True, length=6)

		plt.legend(numpoints=1,loc=('upper right'))
	plt.subplots_adjust(top=0.99, bottom=0.12, left=0.11, right=0.99, hspace=0.0, wspace=0.0)
	plt.savefig(FilePlot, dpi=50, format='eps',rasterized=True)
	plt.show()

def get_plot_energy_vs_beta(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,preskip,postskip,extra_file_name,src_dir,TypePlot,purpose):
	ENT_TYPE=""
	particleA=1
	Units = support.GetUnitConverter()
	kcalmolinvKinv = Units.kcalmoleinvToKelvin

	FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip,postskip,extra_file_name,src_dir,particleA,11)
	
	FileToBePlotEnergy = FilePlotName.SaveEnergy+".txt"
	print(FileToBePlotEnergy)

	trunc = 0
	data_input = genfromtxt(FileToBePlotEnergy, unpack=True, skip_header=0, skip_footer=0)

	#Reading raw PIGS data
	rmlist=[1]
	rmlist=[1]
	for rm_bead in rmlist:
		rm_index=np.where(data_input[1]==rm_bead)[0]
		data_input=np.delete(data_input, rm_index, axis=1)
	print(data_input[1])
	beads_vec=data_input[1]
	valTau=data_input[2]
	valRotEnergy=data_input[3]/kcalmolinvKinv
	valPotEnergy=data_input[4]/kcalmolinvKinv
	valTotalEnergy=data_input[5]/kcalmolinvKinv
	errorRotEnergy=data_input[6]/kcalmolinvKinv
	errorPotEnergy=data_input[7]/kcalmolinvKinv
	errorTotalEnergy=data_input[8]/kcalmolinvKinv

	FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip,postskip,extra_file_name,src_dir,particleA,11)
	FilePlotEnergy = FilePlotName.SaveEnergy
	print(FilePlotEnergy+"-tot.pdf")
	call(["rm", FilePlotEnergy+"-tot.pdf"])

	#PIGS data
	if (extra_file_name == "qTIP4P-"):
		label_str="PIGS/q-TIP4P"
		label_panel = "(a)"
		label_panel = ""
	elif (extra_file_name == "qspcfw-"):
		label_str="PIGS/q-SPC/Fw"
		label_panel = "(c)"
	else:
		label_str="PIGS/TIP4P/2005"
		label_panel = "(a)"
		
	plt.errorbar(valTau[trunc:], valTotalEnergy[trunc:], yerr=errorTotalEnergy[trunc:], color='blue', ls='-', linewidth=1, marker='o', markersize=8, capsize=2, capthick=1, ecolor='blue', label=label_str)

	#For axis lebelling
	plt.xlabel(r'$\beta \ (K^{-1})$', labelpad=5)
	plt.ylabel(r'$E_{0} \ (\mathrm{kcal/mol})$',labelpad=5)

	#For text lebelling
	xmin, xmax = plt.xlim()
	ymin, ymax = plt.ylim()
	plt.xlim(0.0-0.01*(xmax-xmin),xmax)
	midpointy = 0.5*(ymax-ymin)
	deltay = midpointy*0.15
	xmin, xmax = plt.xlim()
	midpointx = 0.5*(xmax-xmin)
	deltax = midpointx*0.15
	textpositionx = xmin+midpointx-0.25*midpointx
	textpositiony = ymin+midpointy

	plt.text(xmin+(xmax-xmin)*0.35,ymax-(ymax-ymin)*0.07,r'$N = \ $'+str(numbmolecules)+'; \ 'r'$r = \ $'+str(Rpt)+r'$\mathrm{\AA}$')
	plt.text(xmin+(xmax-xmin)*0.01,ymax-(ymax-ymin)*0.08,label_panel)

	#For ticks manipulating
	plt.minorticks_on()
	plt.tick_params(axis="both", direction="in", which="minor", right=False, top=False, length=2)
	plt.tick_params(axis="both", direction="in", which="major", right=True, top=True, length=6)

	#Adjust the plot
	plt.legend(numpoints=1,loc=('upper right'))
	if ((Rpt <= 6.0) and (numbmolecules == 2)):
		plt.subplots_adjust(top=0.99,bottom=0.12,left=0.14,right=0.99,hspace=0.0,wspace=0.0)
	elif ((Rpt > 6.0) and (numbmolecules == 2)):
		plt.subplots_adjust(top=0.99,bottom=0.12,left=0.16,right=0.99,hspace=0.0,wspace=0.0)
	elif ((Rpt <= 4.0) and (numbmolecules == 11)):
		plt.subplots_adjust(top=0.99,bottom=0.12,left=0.16,right=0.99,hspace=0.0,wspace=0.0)
	elif ((Rpt > 4.0) and (Rpt <= 5.0) and (numbmolecules == 11)):
		plt.subplots_adjust(top=0.99,bottom=0.12,left=0.13,right=0.99,hspace=0.0,wspace=0.0)
	elif ((Rpt > 5.0) and (Rpt <= 8.0) and (numbmolecules == 11)):
		plt.subplots_adjust(top=0.99,bottom=0.12,left=0.13,right=0.99,hspace=0.0,wspace=0.0)
	elif ((Rpt > 8.0) and (Rpt <= 10.0) and (numbmolecules == 11)):
		plt.subplots_adjust(top=0.99,bottom=0.12,left=0.14,right=0.99,hspace=0.0,wspace=0.0)
	else:
		return
	plt.savefig(FilePlotEnergy+"-tot.pdf", format='pdf')
	plt.show()

def GetFigureEnergyRot_vs_tau(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,preskip,postskip,extra_file_name,src_dir,TypePlot,purpose):
	ENT_TYPE=""
	particleA=1
	Units = support.GetUnitConverter()
	kcalmolinvKinv = Units.kcalmoleinvToKelvin

	Rpt = "{:3.1f}".format(Rpt)
	FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,float(Rpt),gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip,postskip,extra_file_name,src_dir,particleA,11)

	FileToBePlotEnergy = FilePlotName.SaveEnergy+".txt"
	print(FileToBePlotEnergy)

	beads_skip_header_final, beads_skip_footer_final = GetKeyIndices(numbmolecules,molecule_rot,Rpt,FileToBePlotEnergy)	
	data_input = np.genfromtxt(FileToBePlotEnergy, unpack=True)

	if (extra_file_name == "COM-Bisection-Norm-moves-"):
		trunc = 2
		label_str = "PIGS/TIP4P/2005"
		label_panel = "(a)"
	elif (extra_file_name == "COM-and-Bisection-moves-"):
		trunc = 3
		label_str = "PIGS/TIP4P/2005"
		label_panel = "(b)"
	elif (extra_file_name == "qTIP4P-"):
		trunc = beads_skip_header_final
		trunce = int(len(data_input[0]))-beads_skip_footer_final
		print(trunc,trunce)
		label_str = "PIGS/q-TIP4P"
		label_panel = ""
	elif (extra_file_name == "qspcfw-"):
		trunc = 3
		label_str = "PIGS/q-SPC/Fw"
		label_panel = "(c)"
	else:
		trunc = 3
		label_str = "PIGS/TIP4P/2005"
		label_panel = "(a)"

	#Reading raw PIGS data
	print(data_input[1])
	beads_vec=data_input[1]
	valTau=data_input[2]
	valRotEnergy=data_input[3]/kcalmolinvKinv
	valPotEnergy=data_input[4]/kcalmolinvKinv
	valTotalEnergy=data_input[5]/kcalmolinvKinv
	errorRotEnergy=data_input[6]/kcalmolinvKinv
	errorPotEnergy=data_input[7]/kcalmolinvKinv
	errorTotalEnergy=data_input[8]/kcalmolinvKinv

	# Fitting of PIGS data
	fitting_term = "quatric"
	tvar, EnergyFitPlot, ErrorFitPlot = GetFitEnergy(valTau[trunc:trunce], valTotalEnergy[trunc:trunce], errorTotalEnergy[trunc:trunce], variableName, fitting_term)
	energy_fit=EnergyFitPlot[0]
	energy_fit_error=ErrorFitPlot
	print("")
	print(Rpt,' E_O: ',energy_fit,' Error:  ',energy_fit_error)

	FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,float(Rpt),gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip,postskip,extra_file_name,src_dir,particleA,11)
	FilePlotEnergy = FilePlotName.SaveEnergy
	print(FilePlotEnergy+"-tot.pdf")
	call(["rm", FilePlotEnergy+"-tot.pdf"])

	#PIGS data
	plt.errorbar(valTau[trunc:trunce], valTotalEnergy[trunc:trunce], yerr=errorTotalEnergy[trunc:trunce], color='blue', ls='None', linewidth=1, marker='o', markersize=8, capsize=2, capthick=1, ecolor='blue', label=label_str,zorder=0)
	#compared_with=-4.54
	#compared_with_err=0.04
	#plt.errorbar(tvar[0], compared_with, yerr=compared_with_err, color='red', ls='None', linewidth=1,  marker='H', markersize=8, capsize=2, capthick=1, ecolor='red', label=label_str,zorder=5)

	plt.plot(tvar, EnergyFitPlot, color='black', ls='-', linewidth=1 ,zorder=10)
	plt.errorbar(tvar[0], EnergyFitPlot[0], yerr=ErrorFitPlot, color='black', ls='-', linewidth=1,  marker='s', markersize=6, capsize=2, capthick=1, ecolor='black',label="PIGS:Fit",zorder=15)


	# Exact data plotting
	if (numbmolecules == 2): 
		file_exact_val="/Users/tsahoo/ResultsOfExact/ground-state-energy-vs-Rpt-arpack-2-p-H2O-jmax4-grid-9-18-qTIP4P.txt"
		data_exact = genfromtxt(file_exact_val, unpack=True, skip_header=0, skip_footer=0)
		r_data_exact = data_exact[0]
		energy_data_exact = data_exact[1]/kcalmolinvKinv
		index_data_exact = np.where(r_data_exact == float(Rpt))[0][0]
		if (float(Rpt) >= 5.0):
			plt.plot(tvar[0], energy_data_exact[index_data_exact], color='red', ls=None, linewidth=1, marker='D', markersize=6, label="ARPACK lanczos: "+r'$J= \ 4$', zorder=20)
	

	#For axis lebelling
	plt.xlabel(r'$\tau \ (K^{-1})$', labelpad=5)
	plt.ylabel(r'$E_{0} \ (\mathrm{kcal/mol})$',labelpad=5)

	#For text lebelling
	xmin, xmax = plt.xlim()
	ymin, ymax = plt.ylim()
	plt.xlim(0.0-0.01*(xmax-xmin),xmax)
	midpointy = 0.5*(ymax-ymin)
	deltay = midpointy*0.15
	xmin, xmax = plt.xlim()
	midpointx = 0.5*(xmax-xmin)
	deltax = midpointx*0.15
	textpositionx = xmin+midpointx-0.25*midpointx
	textpositiony = ymin+midpointy

	plt.text(xmin+(xmax-xmin)*0.35,ymax-(ymax-ymin)*0.07,r'$N = \ $'+str(numbmolecules)+'; \ 'r'$r = \ $'+str(Rpt)+r'$\mathrm{\AA}$')
	plt.text(xmin+(xmax-xmin)*0.01,ymax-(ymax-ymin)*0.15,label_panel)

	#For ticks manipulating
	plt.minorticks_on()
	plt.tick_params(axis="both", direction="in", which="minor", right=False, top=False, length=2)
	plt.tick_params(axis="both", direction="in", which="major", right=True, top=True, length=6)

	#Adjust the plot
	plt.legend(numpoints=1,loc=('center left'))
	if (numbmolecules == 11):
		if (float(Rpt) == 2.7):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.12,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 2.8) and (float(Rpt) <=3.2)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.14,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 3.3) and (float(Rpt) <=3.4)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.16,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 3.5) and (float(Rpt) <=4.5)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.14,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 4.6) and (float(Rpt) <=4.7)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.16,right=0.99,hspace=0.0,wspace=0.0)
		elif (float(Rpt) == 4.8):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.14,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 5.0) and (float(Rpt) <=6.4)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.13,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 6.6) and (float(Rpt) <=7.0)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.14,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 7.2) and (float(Rpt) <=7.4)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.16,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 7.6) and (float(Rpt) <=9.4)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.14,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 9.6) and (float(Rpt) <=10.0)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.16,right=0.99,hspace=0.0,wspace=0.0)
		else:
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.14,right=0.99,hspace=0.0,wspace=0.0)
	#
	if (numbmolecules == 2):
		if ((float(Rpt) >= 2.7) and (float(Rpt) <=2.8)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.13,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 2.9) and (float(Rpt) <=3.6)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.14,right=0.99,hspace=0.0,wspace=0.0)
		elif (float(Rpt) == 3.7):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.16,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 3.8) and (float(Rpt) <=4.2)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.14,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 4.3) and (float(Rpt) <=4.4)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.16,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 4.5) and (float(Rpt) <=4.6)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.15,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 4.7) and (float(Rpt) <=5.0)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.14,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 5.2) and (float(Rpt) <=5.8)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.14,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 6.0) and (float(Rpt) <=6.4)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.16,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 6.6) and (float(Rpt) <=6.8)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.17,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 7.0) and (float(Rpt) <=8.6)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.16,right=0.99,hspace=0.0,wspace=0.0)
		elif ((float(Rpt) >= 8.8) and (float(Rpt) <=10.0)):
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.17,right=0.99,hspace=0.0,wspace=0.0)
		else:
			plt.subplots_adjust(top=0.99,bottom=0.12,left=0.14,right=0.99,hspace=0.0,wspace=0.0)
	#if (extra_file_name != "qTIP4PF-"):
	#	handles, labels = plt.gca().get_legend_handles_labels()
	#	order = [1,0,2]
	#	plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
	#else:
	#	plt.legend(numpoints=1,loc=('upper right'))
	plt.savefig(FilePlotEnergy+"-tot.pdf", format='pdf')
	plt.show()

def GetFigureEnergyRot_vs_R(TypeCal,molecule_rot,TransMove,RotMove,variableName,RList,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,preskip,postskip,extra_file_name,src_dir,TypePlot,purpose):
	ENT_TYPE=""
	particleA=1
	Units = support.GetUnitConverter()
	kcalmolinvKinv = Units.kcalmoleinvToKelvin

	energy_fit = np.zeros(len(RList))
	energy_fit_error = np.zeros(len(RList))
	i=0
	for Rpt1 in RList:
		Rpt = "{:3.2f}".format(Rpt1)
		parameter1=parameter
		if (float(Rpt) > 8.0):
			parameter1=2.0*parameter
		FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,float(Rpt),gfact,dipolemoment,parameterName,parameter1,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip,postskip,extra_file_name,src_dir,particleA,11)

		FileToBePlotEnergy = FilePlotName.SaveEnergy+".txt"
		print(FileToBePlotEnergy)
		beads_skip_header_final, beads_skip_footer_final = GetKeyIndices(numbmolecules,molecule_rot,Rpt,FileToBePlotEnergy)	

		beads_vec1, valTau1, valTotalEnergy1, errorTotalEnergy1 = np.genfromtxt(FileToBePlotEnergy, unpack=True, usecols=[1, 2, 4, 6], skip_header=0, skip_footer=0)
		trunc=beads_skip_header_final
		trunce=int(len(beads_vec1))-beads_skip_footer_final
		#print(beads_vec1[trunc:trunce])
		valTau = valTau1[trunc:trunce]
		valTotalEnergy = valTotalEnergy1[trunc:trunce]/kcalmolinvKinv
		errorTotalEnergy = errorTotalEnergy1[trunc:trunce]/kcalmolinvKinv

		# Fitting of PIGS data
		fitting_term = "quatric" 
		tvar, EnergyFitPlot, ErrorFitPlot = GetFitEnergy(valTau, valTotalEnergy, errorTotalEnergy, variableName, fitting_term)
		energy_fit[i]=EnergyFitPlot[0]
		energy_fit_error[i]=ErrorFitPlot
		#print(Rpt,'  ',energy_fit[i],'  ',energy_fit_error[i])
		i+=1

	#Plot begins
	FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,1.0,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip,postskip,extra_file_name,src_dir,particleA,11)
	FilePlotEnergy = FilePlotName.SaveEnergyFitvsR
	print(FilePlotEnergy+"-tot.pdf")
	call(["rm", FilePlotEnergy+"-tot.pdf"])

	# Save fitted date into a file
	fit_data = np.array([RList, energy_fit/(numbmolecules-1), energy_fit_error/(numbmolecules-1)])
	dir_file_saved='/Users/tsahoo/academic-writting/for-academic-position/raw-data-pigs/'
	np.savetxt(dir_file_saved+'raw-data-fig12-jcp-154-244305-2021-red-line.txt',fit_data.T,fmt='%20.8f', delimiter=' ',header='x y1 y2',comments='# First, second and third columns are the intermolecular distance in angstrom, ground-state energy of 11 p-H2O molecules in kcal/mole and the associated error, respectively. See the Figure 12 of J. Chem. Phys. 154, 244305 (2021). The next line is for Tikz plot.\n')
	# Plotting of fitting data
	plt.errorbar(RList, energy_fit, yerr=energy_fit_error, color='black', ls=None, linewidth=1,  marker='o', markersize=8, capsize=2, capthick=1, ecolor='black', label="PIGS", zorder=0)

	if (numbmolecules == 2):
		colorList = ['r', 'm', 'g', 'r', 'b']
		markerList = ['s', 'p', '>', 's', '<', '8', 'p']
		lsList = ['-', '--', '-.', ':']
		iplot=0
		rmlist = [5.1+0.2*i for i in range(25)]

		dir_ed_results = '/Users/tsahoo/academic-project/outputs/exact-computation/'	
		FileToBePlotEnergy_exact=dir_ed_results+'ground-state-energy-vs-Rpt-lanc-2-p-H2O-jmax9-grid-20-20-niter600-qtip4p-correct.txt'
		print(FileToBePlotEnergy_exact)
		data_input = genfromtxt(FileToBePlotEnergy_exact, unpack=True, skip_header=0, skip_footer=0)
		'''
		for rm_val in rmlist:
			rm_val = "{:3.1f}".format(rm_val)
			rm_index=np.where(data_input[0]==float(rm_val))[0]
			data_input=np.delete(data_input, rm_index, axis=1)
		'''
		Rpt_exact = data_input[0]
		TotalEnergy_exact = data_input[1]/kcalmolinvKinv
		plt.plot(Rpt_exact, TotalEnergy_exact, color=colorList[iplot], ls=lsList[0], linewidth=1, marker=markerList[iplot], markersize=6, label="ED")
		iplot = iplot+1

	#For axis lebelling
	plt.xlabel(r'$r \ (\mathrm{\AA})$', labelpad=5)
	plt.ylabel(r'$E_{0} \ (\mathrm{kcal/mol})$',labelpad=5)

	#For text lebelling
	xmin, xmax = plt.xlim()
	ymin, ymax = plt.ylim()
	midpointy = 0.5*(ymax-ymin)
	deltay = midpointy*0.15
	xmin, xmax = plt.xlim()
	midpointx = 0.5*(xmax-xmin)
	deltax = midpointx*0.15
	textpositionx = xmin+midpointx-0.25*midpointx
	textpositiony = ymin+midpointy

	plt.text(xmin+(xmax-xmin)*0.40,ymax-(ymax-ymin)*0.52,r'$N = \ $'+str(numbmolecules))
	#plt.text(xmin+(xmax-xmin)*0.01,ymax-(ymax-ymin)*0.15,label_panel)

	#For ticks manipulating
	plt.minorticks_on()
	plt.tick_params(axis="both", direction="in", which="minor", right=False, top=False, length=2)
	plt.tick_params(axis="both", direction="in", which="major", right=True, top=True, length=6)

	#Adjust the plot
	plt.legend(numpoints=1,loc=('center right'))
	if (numbmolecules == 2):
		plt.subplots_adjust(top=0.99,bottom=0.12,left=0.11,right=0.99,hspace=0.0,wspace=0.0)
	if (numbmolecules == 11):
		plt.subplots_adjust(top=0.99,bottom=0.12,left=0.12,right=0.99,hspace=0.0,wspace=0.0)
	plt.savefig(FilePlotEnergy+"-tot.pdf", format='pdf')
	plt.show()

## Below three functions 'fitFunc4, fitFunc3, and GetFitEnergy' are for fitting the energy obtained by PIGS
def fitFunc4(var, a, b, c):
		return a+b*var*var+c*var*var*var*var

def fitFunc3(var, a, b):
		return a+b*var*var

def GetFitEnergy(val1, val2, val3, variableName, fitting_term):
	t = np.linspace(0.0, np.max(val1), num=pow(10, 4))
	if (fitting_term == "quatric"):
		initialParameters = np.array([-1e4, -1e+8, 1e+12])
		popt, pcov = curve_fit(fitFunc4, val1, val2, sigma=val3, p0=initialParameters)
		perr = np.sqrt(np.diag(pcov))
		fit = fitFunc4(t, *popt)
	if (fitting_term == "quadratic"):
		initialParameters = np.array([-1e4, -1e+8])
		popt, pcov = curve_fit(fitFunc3, val1, val2, sigma=val3, p0=initialParameters)
		perr = np.sqrt(np.diag(pcov))
		fit = fitFunc3(t, *popt)

	if (variableName == "tau"):
		return t, fit, perr[0]
	else:
		return fit[0], perr[0]

