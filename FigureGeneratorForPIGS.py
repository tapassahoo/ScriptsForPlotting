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
from numpy import *
from pylab import *
from scipy.optimize import curve_fit

import mypkg.pkgMoribs.support_without_parallel as support

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

def FigureENT(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, TypePlot, beadsRef):
	font = 28
	fontlegend = font/2.0
	preskip = 0
	postskip = 0

	plotnum = 2

	BConstant = support.GetBconst(molecule_rot)  # in wavenumber
	Units = support.GetUnitConverter()
	BConstantK = BConstant*Units.CMRECIP2KL

	if (((TypePlot == "RFACTOR") or (TypePlot == "GFACTOR")) and variableName == "tau"):
		font = 28
		fontlegend = font/2
		fig = plt.figure(figsize=(8, 6))

		var = beadsRef
		FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, -1.0, -1.0, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, var)

		if (TypePlot == "RFACTOR"):
			FilePlotEntropy = FilePlotName.SaveEntropyRFAC+str(plotnum)+".eps"
		if (TypePlot == "GFACTOR"):
			FilePlotEntropy = FilePlotName.SaveEntropyGFAC+str(plotnum)+".eps"
		outfileEntropy = FilePlotEntropy
		call(["rm", FilePlotEntropy])
		print(outfileEntropy)
#
		if (plotnum != 2):
			plt.axhline(y=math.log(2.0), color='blue', lw=2.0, linestyle='-', label='ln(2)')
#
		if (plotnum != 2):
			nn = [2, 4, 6]
		else:
			nn = [16, 32]
		for numbmolecules in nn:
			particleA = int(numbmolecules/2)
			if (numbmolecules == 2):
				beadsRef = 101
				DList = [1.0+0.5*i for i in range(7)]
				DList += [4.5, 5.0, 5.5, 6.0]
				numbblocks = 50000
				postskip1 = 30000
			if (numbmolecules == 4):
				DList = [1.0+0.5*i for i in range(7)]
				numbblocks = 20000
			if (numbmolecules == 6):
				DList = [1.0+0.25*i for i in range(10)]
				numbblocks = 20000
			if (numbmolecules == 8):
				DList = [0.9106, 1.1152, 1.2878, 1.4398, 1.5772, 1.7036, 1.8212, 1.9317, 2.0362]
				#DList	= [1.5772, 1.7036, 1.8212, 1.9317, 2.0362]
				numbblocks = 20000
				preskip1 = 8000
				postskip1 = 10000
			if (numbmolecules == 16):
				#DList	= [1.5772, 1.7036, 1.8212, 1.9317, 2.0362]
				DList = [0.9106, 1.1152, 1.2878, 1.4398, 1.5772, 1.7036, 1.8212, 1.9317, 2.0362]
				numbblocks = 10000
				preskip1 = 8000
				postskip1 = 0
			if (numbmolecules == 32):
				#DList	= [1.5772, 1.7036, 1.8212, 1.9317, 2.0362]
				DList = [0.9106, 1.1152, 1.2878, 1.4398, 1.5772, 1.7036, 1.8212, 1.9317, 2.0362]
				numbblocks = 10000
				preskip1 = 6000
				postskip1 = 0
			RFactorPlot = np.zeros(len(DList))
			entropy1Plot = np.zeros(len(DList))
			purity1Plot = np.zeros(len(DList))
			err_entropy1Plot = np.zeros(len(DList))
			err_purity1Plot = np.zeros(len(DList))
			entropy2Plot = np.zeros(len(DList))
			entropy3Plot = np.zeros(len(DList))
			iii = 0
			for dipolemoment in DList:
				RFactorList = support.GetrAndgFactor(
					molecule_rot, Rpt, dipolemoment)
				RFactor = RFactorList[0]
				if (dipolemoment > 4.0):
					if (numbmolecules == 2):
						numbblocks = 20000
						postskip1 = 0
					else:
						numbblocks = 20000
				FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, -1.0, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, beadsRef)
				FileToBePlotEntropy = FilePlotName.SaveEntropy+".txt"
				#FileToBePlotED    = FilePlotName.SaveEntropyED+".txt"

				beads1, var1, purity1, entropy1, err_purity1, err_entropy1 = genfromtxt(FileToBePlotEntropy, unpack=True, usecols=[0, 1, 4, 5, 8, 9], skip_header=preskip, skip_footer=postskip)
				if (numbmolecules <= 6):
					RFactor, energy3, entropy3 = loadtxt(FileToBePlotED, unpack=True, usecols=[0, 1, 2])

				if (TypePlot == "GFACTOR"):
					RFactorPlot[iii] = 1.0/(RFactor*RFactor*RFactor)
				if (TypePlot == "RFACTOR"):
					RFactorPlot[iii] = RFactor
				if (plotnum != 2):
					entropy3Plot[iii] = entropy3
				if ((numbmolecules == 4) and (dipolemoment == 4.0)):
					beadsRef = 61
					beadsRef1 = 81
				elif ((numbmolecules == 6) and (dipolemoment == 3.25)):
					beadsRef = 21
					beadsRef1 = 81
				elif ((numbmolecules == 6) and (dipolemoment == 3.5)):
					beadsRef = 101
					beadsRef1 = 81
				else:
					beadsRef = 81
					beadsRef1 = 81
				if (numbmolecules == 2):
					beadsRef = 101
					beadsRef1 = beadsRef
				if ((numbmolecules == 2) and (dipolemoment == 6.0)):
					beadsRef = 41
					beadsRef1 = 41
				if (numbmolecules == 8):
					beadsRef = 21
					beadsRef1 = beadsRef
				if (numbmolecules == 16):
					beadsRef = 21
					beadsRef1 = beadsRef
				if (numbmolecules == 32):
					beadsRef = 21
					beadsRef1 = beadsRef

				if (np.isscalar(entropy1) == True):
					entropy1Plot[iii] = entropy1
					err_entropy1Plot[iii] = err_entropy1
					purity1Plot[iii] = purity1
					err_purity1Plot[iii] = err_purity1
				else:
					ii = 0
					for i in beads1:
						indexi = int(i+0.5)
						beads = indexi
						if beads == beadsRef:
							entropy1Plot[iii] = entropy1[ii]
							purity1Plot[iii] = purity1[ii]
							err_purity1Plot[iii] = err_purity1[ii]

						ii += 1
					ii1 = 0
					for i in beads1:
						indexi = int(i+0.5)
						beads = indexi
						if beads == beadsRef1:
							err_entropy1Plot[iii] = err_entropy1[ii1]

						ii1 += 1
				iii += 1

			print("S2:	PIGS "+str(numbmolecules))
			print(entropy1Plot)
			print("S2:	ED ")
			print(entropy3Plot)
#
			plotEntropyENT1(numbmolecules, RFactorPlot, entropy1Plot, err_entropy1Plot, variableName, RFactorPlot, entropy2Plot, entropy3Plot, font, TypePlot)
		ymin, ymax = plt.ylim()
		xmin, xmax = plt.xlim()
		if plotnum == 1:
			plt.ylabel(r'$S_{2}$', fontsize=font)
			if ymin < 0.0:
				plt.ylim(0.0, 0.82)
			if xmin < 0.0:
				plt.xlim(0.0, xmax)
			plt.xticks(np.arange(0, 9, step=2))
			plt.yticks(np.arange(0, 0.9, step=0.2))
			Text1 = "(a)"
		if plotnum == 2:
			plt.xlim(0.192, 1.01)
			plt.ylim(-0.46, 0.64)
			plt.xticks(np.arange(0.2, 1.1, step=0.2))
			plt.yticks(np.arange(-0.4, 0.6, step=0.2))
			Text1 = r'$\textrm{(b)}$'
			Text1 = ''
			plt.ylabel(r'$S_{2}$', fontsize=font, labelpad=-14)
		Text2 = ""
		if Text1:
			PlotLabel(Text1, Text2, font, xmin, xmax, ymin, ymax, variableName, parameter, numbmolecules, molecule, Rpt, dipolemoment)

		plt.subplots_adjust(top=0.97, bottom=0.14, left=0.14, right=0.98, hspace=0.0, wspace=0.0)
		if (plotnum != 2):
			plt.legend(bbox_to_anchor=(0.68, 0.75), loc=2, borderaxespad=1., shadow=True, fontsize=fontlegend)
		else:
			plt.legend(bbox_to_anchor=(0.98, 0.35), borderaxespad=1., shadow=True, fontsize=fontlegend)
		plt.savefig(outfileEntropy, dpi=200, format='eps')
		plt.show()

	if (TypePlot == "S2"):
		font = 28
		fontlegend = font/2
		iFrame = 1
		iFigLabel = 0
		FigureLabel = [r'$\mathrm{(a)}$', r'$\mathrm{(b)}$']
		variableList = ["beta", "tau"]
		for variableName in variableList:
			if variableName == "beta":
				parameterName = "tau"
				parameter = 0.005
				postskip = 8
			if variableName == "tau":
				parameterName = "beta"
				parameter = 0.2
				postskip = 0
##
			fig = plt.figure(figsize=(8, 6))
			plt.grid(True)
##
			FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, beadsRef)
			FileToBePlot = FilePlotName.SaveEntropy+".txt"
			FileToBePlotMM = FilePlotName.SaveEntropyMM+".txt"
			FileToBePlotED = FilePlotName.SaveEntropyED+".txt"
			print(FileToBePlot)
			print(FileToBePlotMM)
			print(FileToBePlotED)
##
			FilePlot = FilePlotName.SaveEntropy+".eps"
			outfile = FilePlot
			call(["rm", FilePlot])
			print(FilePlot)
##
			var2 = 0.0
			entropy2 = [0.0, 0.0]
			entropy3 = 0.0
###================data reading begin==================###
			var1, purity1, entropy1, err_purity1, err_entropy1 = genfromtxt(FileToBePlot, unpack=True, usecols=[1, 4, 5, 8, 9], skip_header=preskip, skip_footer=postskip)
			if (numbmolecules <= 4):
				beads, entropy2 = genfromtxt(FileToBePlotMM, unpack=True, usecols=[0, 2], skip_header=preskip, skip_footer=postskip)
				if (variableName == "tau"):
					var2 = parameter/(beads-1.0)

				if (variableName == "beta"):
					var2 = parameter*(beads-1.0)

			if (numbmolecules <= 6):
				RFactor, entropy3 = loadtxt(
					FileToBePlotED, unpack=True, usecols=[0, 2])
###================data reading end====================###
			print(entropy1)
			print(entropy2)
			if variableName == "tau":
				var2 = var2[1:]
				entropy2 = entropy2[1:]

			#plt.subplot(1, 2, iFrame)
			iFrame = iFrame + 1
			plotEntropyENT(var1, entropy1, err_entropy1, variableName, var2, entropy2, entropy3, font, font)
			plt.ylabel(r'$S_{2}$', fontsize=font)
			ymin, ymax = plt.ylim()
			xmin, xmax = plt.xlim()
			Text1 = FigureLabel[iFigLabel]
			iFigLabel = iFigLabel + 1
			RGFactor = support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
			gFactor = RGFactor[1]
			arg = "%3.2f" % gFactor
			Text2 = str(arg)
			PlotLabel(Text1, Text2, font, xmin, xmax, ymin, ymax, variableName, parameter, numbmolecules, molecule, Rpt, dipolemoment)
			plt.legend(bbox_to_anchor=(0.80, 0.50), loc=2, borderaxespad=0., shadow=True, fontsize=fontlegend)
			plt.subplots_adjust(top=0.98, bottom=0.16, left=0.18, right=0.95, hspace=0.0, wspace=0.)
			plt.savefig(outfile, dpi=200, format='eps')

			call(["open", outfile])
			# plt.show()


def FigureENTCOMBINE(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, TypePlot, beadsRef):
	font = 28
	fontlegend = font/2.0
	preskip = 0
	postskip = 0

	if (variableName == "tau"):
		fig = plt.figure(figsize=(8, 6))
		plt.grid(True)

		var = beadsRef
		FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, var)
		FilePlotEntropy = FilePlotName.SaveEntropyCOMBO+".eps"
		outfileEntropy = FilePlotEntropy
		call(["rm", FilePlotEntropy])
		print(outfileEntropy)
#
		DList = [1.0+0.5*i for i in range(7)]
		DList += [4.5, 5.0, 5.5, 6.0]
		labelIndex = 0
		TypeList = ["SWAPTOUNSWAP", "BROKENPATH"]
		for ENT_TYPE in TypeList:
			beadsRef1 = beadsRef
			if (ENT_TYPE == "SWAPTOUNSWAP"):
				numbblocks = 50000
				postskip1 = 30000
			if (ENT_TYPE == "BROKENPATH"):
				numbblocks = 20000
				postskip1 = 0

			RFactorPlot = np.zeros(len(DList))
			entropy1Plot = np.zeros(len(DList))
			err_entropy1Plot = np.zeros(len(DList))
			entropy2Plot = np.zeros(len(DList))

			iii = 0
			for dipolemoment in DList:
				RFactorList = support.GetrAndgFactor(
					molecule_rot, Rpt, dipolemoment)
				if (dipolemoment > 4.0):
					numbblocks = 20000
					postskip1 = 0
					if ((ENT_TYPE == "SWAPTOUNSWAP") and (dipolemoment == 6.0)):
						beadsRef1 = 41

				FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, beadsRef1)
				FileToBePlotEntropy = FilePlotName.SaveEntropy+".txt"

				if (ENT_TYPE == "SWAPTOUNSWAP"):
					FileToBePlotED = FilePlotName.SaveEntropyED+".txt"
					beads1, var1, purity1, entropy1, err_purity1, err_entropy1 = genfromtxt(FileToBePlotEntropy, unpack=True, usecols=[0, 1, 4, 5, 8, 9], skip_header=preskip, skip_footer=postskip)
					entropyED = loadtxt(FileToBePlotED, unpack=True, usecols=[2])
				if (ENT_TYPE == "BROKENPATH"):
					beads1, var1, entropy1, err_entropy1 = genfromtxt(FileToBePlotEntropy, unpack=True, usecols=[0, 1, 4, 7], skip_header=preskip, skip_footer=postskip)
				RFactorPlot[iii] = RFactorList[1]

				ii = 0
				for i in beads1:
					indexi = int(i+0.5)
					beads = indexi
					if beads == beadsRef1:
						entropy1Plot[iii] = entropy1[ii]
						entropy2Plot[iii] = entropyED
						err_entropy1Plot[iii] = err_entropy1[ii]
					ii += 1
				iii += 1

			print("S2:	PIGS ")
			print(entropy1Plot)
			print(RFactorPlot)
#
			colorList = ['red', 'blue']
			lsList = ['--', '-.']
			markerList = ['^', 'v']
			labelList = ['Swap+Unswap grand ensemble', 'Broken path ensemble']

			if (ENT_TYPE == "SWAPTOUNSWAP"):
				plt.plot(RFactorPlot, entropy2Plot, color='black', ls='-', linewidth=1,  marker='o', markersize=9, label='ED')
			plt.errorbar(RFactorPlot, entropy1Plot, yerr=err_entropy1Plot, color=colorList[labelIndex], ls=lsList[labelIndex], linewidth=1,	marker=markerList[labelIndex], markersize=8, label=labelList[labelIndex])

			labelIndex += 1

			ymin, ymax = plt.ylim()
			plt.ylim(-0.001, 0.901)
			xmin, xmax = plt.xlim()
			plt.xlim(0, 9)
			Text1 = ""
			Text2 = ""
			if Text1:
				PlotLabel(Text1, Text2, font, xmin, xmax, ymin, ymax, variableName, parameter, numbmolecules, molecule, Rpt, dipolemoment)
			plt.xticks(np.arange(0, 10, step=1), fontsize=font, rotation=0)
			plt.yticks(np.arange(0.0, 0.91, step=0.1), fontsize=font, rotation=0)

		plt.ylabel(r'$S_{2}$', fontsize=font)
		plt.xlabel(r'$g$', fontsize=font, labelpad=-3)

		plt.subplots_adjust(top=0.97, bottom=0.14, left=0.14, right=0.98, hspace=0.0, wspace=0.4)
		plt.legend(bbox_to_anchor=(0.40, 0.45), loc=2, borderaxespad=1., shadow=True, fontsize=fontlegend)
		plt.savefig(outfileEntropy, dpi=200, format='eps')

		call(["open", outfileEntropy])


def plotEntropyENT1(numbmolecules, var, val, err_val, variableName, var1, val1, val2, font, TypePlot):
	if (numbmolecules == 2):
		plt.errorbar(var, val, yerr=err_val, color="red", ls='-', linewidth=1,  marker="o", markersize=8, label='PIGS: '+r'$N=2$')
		plt.plot(var, val2, color='black', ls='None', linewidth=1, marker="o", markersize=10, label='ED: N=2')
		#plt.plot(var, val2, color = 'black', ls = '-', linewidth=3, marker = "o", markersize = 12, label = 'ED: 2 HF')
		#plt.plot(DMRGdatagFac, DMRGdataPlot, linestyle = 'None', color = 'darkcyan', marker = "o", markersize = 10, label = 'DMRG: 2 HF')
		#plt.plot(var, val1, linestyle = 'None', color = 'mediumpurple', marker = "o", markersize = 10, label = 'MM: 2 HF')

	if (numbmolecules == 4):
		plt.errorbar(var, val, yerr=err_val, color="red", ls='-', linewidth=1,  marker="s", markersize=8, label='PIGS: N=4')
		plt.plot(var, val2, color='black', ls='None', linewidth=1, marker="s", markersize=10, label='ED: N=4')
		#plt.plot(DMRGdatagFac, DMRGdataPlot, linestyle = 'None', color = 'darkcyan', marker = "s", markersize = 8, label = 'DMRG: 4 HF')

	if (numbmolecules == 6):
		plt.errorbar(var, val, yerr=err_val, color="red", ls='-', linewidth=1,  marker="v", markersize=8, label='PIGS: N=6')
		plt.plot(var, val2, color='black', ls='None', linewidth=1, marker="v", markersize=10, label='ED: N=6')
		#plt.plot(DMRGdatagFac, DMRGdataPlot, linestyle = 'None', color = 'darkcyan', marker = "8", markersize = 7, label = 'DMRG: 6 HF')

	if (numbmolecules == 8):
		plt.errorbar(var, val, yerr=err_val, color="red", ls='-', linewidth=1,  marker="p", markersize=8, label='PIGS: N=8')
		#plt.plot(var, val2, color = 'black', ls = 'None', linewidth=1, marker = "v", markersize = 10, label = 'ED: N=6')
		#plt.plot(DMRGdatagFac, DMRGdataPlot, linestyle = 'None', color = 'darkcyan', marker = "8", markersize = 7, label = 'DMRG: 6 HF')

	if (numbmolecules == 16):
		plt.errorbar(var, val, yerr=err_val, color="blue", ls='-', linewidth=1,  marker="8", markersize=8, label='PIGS: '+r'$N=16$')
		#plt.plot(var, val2, color = 'black', ls = 'None', linewidth=1, marker = "v", markersize = 10, label = 'ED: N=6')
		#plt.plot(DMRGdatagFac, DMRGdataPlot, linestyle = 'None', color = 'darkcyan', marker = "8", markersize = 7, label = 'DMRG: 6 HF')

	if (numbmolecules == 32):
		plt.errorbar(var, val, yerr=err_val, color="green", ls='-', linewidth=1,  marker="^", markersize=8, label='PIGS: '+r'$N=32$')
		#plt.plot(var, val2, color = 'black', ls = 'None', linewidth=1, marker = "v", markersize = 10, label = 'ED: N=6')
		#plt.plot(DMRGdatagFac, DMRGdataPlot, linestyle = 'None', color = 'darkcyan', marker = "8", markersize = 7, label = 'DMRG: 6 HF')

	ymin, ymax = plt.ylim()
	midpointy = 0.5*(ymax-ymin)
	deltay = midpointy*0.15
	xmin, xmax = plt.xlim()
	midpointx = 0.5*(xmax-xmin)
	deltax = midpointx*0.15
	textpositionx = xmin+midpointx-0.25*midpointx
	textpositiony = ymin+midpointy
	plt.xticks(fontsize=font, rotation=0)
	plt.yticks(fontsize=font, rotation=0)

	if (TypePlot == "RFACTOR"):
		plt.xlabel(r'$R$', fontsize=font, labelpad=0)
	if (TypePlot == "GFACTOR"):
		plt.xlabel(r'$g$', fontsize=font, labelpad=-3)


def plotEntropyENT(var, val, err_val, variableName, var1, val1, val2, font, fontlegend):
	plt.errorbar(var, val, yerr=err_val, color='red', ls='-', linewidth=1,  marker="o", markersize=7, label='PIGS')
	if (val1[0] != 0.0):
		plt.plot(var1, val1, linestyle='--', linewidth=2, color='blue', marker="v", markersize=7, label='MM')

	if (variableName == "tau"):
		label_xtics = [0.00]
		plt.xlim(0, 0.0201)
		if (val2 != 0.0):
			plt.axhline(y=val2, color='black', lw=2.0, linestyle='-', label='ED')
		plt.xticks(np.arange(0.0, 0.021, step=0.005), fontsize=fontlegend, rotation=0)
		plt.ylim(0.0399, 0.11001)
		plt.yticks(np.arange(0.04, 0.11, step=0.02), fontsize=fontlegend, rotation=0)
	else:
		plt.xlim(0.0195, 0.12001)
		plt.xticks(np.arange(0.02, 0.13, step=0.02), fontsize=fontlegend, rotation=0)
		plt.ylim(0.02499, 0.06001)
		plt.yticks(np.arange(0.025, 0.060, step=0.01), fontsize=fontlegend, rotation=0)

	if (variableName == "beta"):
		plt.xlabel(r'$\beta \ \  (\mathrm{K^{-1}})$', fontsize=font)
	if (variableName == "tau"):
		plt.xlabel(r'$\tau \ \	(\mathrm{K^{-1}})$', fontsize=font)


def FigureCorrelation(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, beadsRef, RefPointList):

	font = 40
	fontlegend = font/2.0

	colorList = ['r', 'g', 'b', 'm']
	markerList = ['o', 's', 'p', '<', '>', 's', '8', 'p']
	lsList = ['-', '--', '-.', ':']

	TextLabel = [r'$\mathrm{(a)}$']
	TextLabel += [r'$\mathrm{(b)}$']
	TextLabel += [r'$\mathrm{(c)}$']
	TextLabel += [r'$\mathrm{(d)}$']
	TextLabel += [r'$\mathrm{(e)}$']
	TextLabel += [r'$\mathrm{(f)}$']
	TextLabel += [r'$\mathrm{(g)}$']
	TextLabel += [r'$\mathrm{(h)}$']
	TextLabel += [r'$\mathrm{(i)}$']
	TextLabel += [r'$\mathrm{(j)}$']
	DList = [1.0+i for i in range(4)]
	DList = [1.25, 2.25,  3.25]
	iTextLabel = 0
	for dipolemoment in DList:
		RGFactor = support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
		gFactor = RGFactor[1]
		arg = "%3.2f" % gFactor

		TypeCorrList = ["TotalCorr", "ZCorr", "XYCorr"]
		for TypeCorr in TypeCorrList:
			FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName,
												   parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, beadsRef)
			if TypeCorr == "TotalCorr":
				FileToBePlot = FilePlotName.SaveTotalCorr+".txt"
				FilePlot = FilePlotName.SaveTotalCorr + \
					"-ithRotor"+str(RefPointList[0])+".eps"
			elif TypeCorr == "ZCorr":
				FileToBePlot = FilePlotName.SaveZCorr+".txt"
				FilePlot = FilePlotName.SaveZCorr + \
					"-ithRotor"+str(RefPointList[0])+".eps"
			elif TypeCorr == "XYCorr":
				FileToBePlot = FilePlotName.SaveXYCorr+".txt"
				FilePlot = FilePlotName.SaveXYCorr + \
					"-ithRotor"+str(RefPointList[0])+".eps"
			print(FilePlot)
			call(["rm", FilePlot])
			datacorr = genfromtxt(FileToBePlot)
			if (datacorr.ndim == 2):
				index = 0
				for i in range(len(datacorr)):
					if datacorr[i, 0] == beadsRef:
						findex = index
						break
					index = index+1

			FuncCorr = np.zeros((numbmolecules, numbmolecules))
			ErrorFuncCorr = np.zeros((numbmolecules, numbmolecules))
			ii = 0
			for i in range(numbmolecules):
				for j in range(i, numbmolecules):

					nc = 2+(2*ii)
					nec = nc+1
					if datacorr.ndim == 2:
						FuncCorr[i, j] = datacorr[findex, nc]
						ErrorFuncCorr[i, j] = datacorr[findex, nec]
					else:
						FuncCorr[i, j] = datacorr[nc]
						ErrorFuncCorr[i, j] = datacorr[nec]
					if (j != i):
						FuncCorr[j, i] = FuncCorr[i, j]
						ErrorFuncCorr[j, i] = ErrorFuncCorr[i, j]
					ii = ii+1
			fig = plt.figure(figsize=(8, 6))
			iRef = 0
			for RefPoint in RefPointList:
				val1 = np.arange(1, numbmolecules+1)
				val2 = np.zeros(numbmolecules)
				val3 = np.zeros(numbmolecules)
				for j in range(numbmolecules):
					val2[j] = FuncCorr[RefPoint, j]
					val3[j] = ErrorFuncCorr[RefPoint, j]

				value = "%3.2f" % val2[0]
				# plt.errorbar(val1[1:], val2[1:], yerr=val3[1:], color = colorList[iRef], ls = lsList[iRef], linewidth=3, marker = markerList[iRef], markersize = 20)	#excluding selfcorrelation
				# excluding selfcorrelation
				plt.errorbar(
					val1, val2, yerr=val3, color=colorList[iRef], ls=lsList[iRef], linewidth=3, marker=markerList[iRef], markersize=20)

				ymin, ymax = plt.ylim()
				midpointy = 0.5*(ymax-ymin)
				deltay = midpointy*0.15
				xmin, xmax = plt.xlim()
				plt.xlim(xmin-0.1, xmax+0.1)
				midpointx = 0.5*(xmax-xmin)
				deltax = midpointx*0.15
				textpositionx = xmin+midpointx-0.25*midpointx
				textpositiony = ymin+midpointy

				Text1 = TextLabel[iTextLabel]
				Text2 = str(arg)
				if (TypeCorr == "TotalCorr"):
					plt.ylabel(
						r'$\mathrm{C}_{'+str(RefPoint+1)+',j}$', fontsize=font)
					#plt.text(textpositionx-2*deltax, textpositiony-6*deltay, r'$\mathrm{C}_{'+str(RefPoint+1)+',1} = $'+str(value), fontsize=font)
					if Text1:
						plt.text(textpositionx-(4.75*deltax), textpositiony+(5.5*deltay), Text1, fontsize=font)

					if Text2:
						plt.text(textpositionx+(3.5*deltax), textpositiony + 5.5*deltay, r'$g = $'+Text2, fontsize=font)

				if (TypeCorr == "ZCorr"):
					plt.ylabel(
						r'$\mathrm{C}^{\parallel}_{'+str(RefPoint+1)+',j}$', fontsize=font)
					# if dipolemoment == 1.25:
					#	plt.text(textpositionx-0*deltax, textpositiony-2*deltay,  r'$\mathrm{C}^{\parallel}_{'+str(RefPoint+1)+',1} = $'+str(value), fontsize=font)
					# else:
					#	plt.text(textpositionx-2*deltax, textpositiony-6*deltay,  r'$\mathrm{C}^{\parallel}_{'+str(RefPoint+1)+',1} = $'+str(value), fontsize=font)
					if Text1:
						'''
						if dipolemoment == 1.25:
								plt.text(textpositionx-(4.05*deltax), textpositiony+(5.5*deltay), Text1, fontsize=font)
						else:
								plt.text(textpositionx-(4.75*deltax), textpositiony+(3.5*deltay), Text1, fontsize=font)
						'''
						plt.text(textpositionx-(4.75*deltax), textpositiony+(5.5*deltay), Text1, fontsize=font)

					if Text2:
						plt.text(textpositionx+(3.5*deltax), textpositiony + 5.5*deltay, r'$g = $'+Text2, fontsize=font)

				if (TypeCorr == "XYCorr"):
					plt.ylabel(
						r'$\mathrm{C}^{\bot}_{'+str(RefPoint+1)+',j}$', fontsize=font)
					plt.axhline(y=0.0, color='black', lw=3.0, linestyle='--')
					#plt.text(textpositionx-2*deltax, textpositiony-6*deltay, r'$\mathrm{C}^{\bot}_{'+str(RefPoint+1)+',1} = $'+str(value), fontsize=font)
					if Text1:
						plt.text(textpositionx-(4.75*deltax), textpositiony+(5.5*deltay), Text1, fontsize=font)
						#plt.text(textpositionx-(4.75*deltax), textpositiony+(3.5*deltay), Text1, fontsize=font)

					if Text2:
						#plt.text(textpositionx+(3.5*deltax), textpositiony+3.0*deltay, r'$g = $'+Text2, fontsize=font)
						plt.text(textpositionx+(3.5*deltax), textpositiony + 5.5*deltay, r'$g = $'+Text2, fontsize=font)

				iRef += 1

			# plt.grid(True)

			plt.xlabel(r'$j$', fontsize=font)

			iTextLabel += 1

			plt.xticks(fontsize=font, rotation=0)
			plt.yticks(fontsize=font, rotation=0)
			plt.subplots_adjust(top=0.96, bottom=0.23, left=0.30, right=0.95, hspace=0.6, wspace=1.0)
			#plt.subplots_adjust(top=0.96, bottom=0.18, left=0.22, right=0.95, hspace=0.6, wspace=1.0)
			#plt.legend(bbox_to_anchor=(0.58, 0.99), loc=2, borderaxespad=0., shadow=True, fontsize = fontlegend)
			plt.savefig(FilePlot, dpi=200, format='eps')

			#call(["okular", FilePlot])
			call(["open", FilePlot])
			# plt.show()

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


def PlotLabel(Text1, Text2, font, xmin, xmax, ymin, ymax, variableName, parameter, numbmolecules, molecule, Rpt, dipolemoment):
	midpointy = 0.5*(ymax-ymin)
	deltay = midpointy*0.15
	midpointx = 0.5*(xmax-xmin)
	deltax = midpointx*0.15
	textpositionx = xmin+midpointx-0.25*midpointx
	textpositiony = ymin+midpointy

	if Text1:
		plt.text(textpositionx-(4.45*deltax), textpositiony + (5.0*deltay), Text1, fontsize=font)
		#plt.text(textpositionx-(4.70*deltax), textpositiony+(5.5*deltay), Text1, fontsize=font)

	if Text2:
		plt.text(textpositionx+(0.0*deltax), textpositiony + 0*deltay, r'$g = '+Text2+'$', fontsize=font)
		#plt.text(textpositionx+(4.3*deltax), textpositiony+5*deltay, r'$g$ = '+Text2, fontsize=font)

	'''
	if (variableName == "beta" and Text2):
		plt.text(textpositionx, textpositiony+0*deltay, 'Parameters:', fontsize=font)
		plt.text(textpositionx, textpositiony-1*deltay, r'$\mathrm{N} =$'+str(numbmolecules)+" "+molecule, fontsize=font)
		plt.text(textpositionx, textpositiony-2*deltay, r'$\mathrm{R} =$'+str(Rpt)+r'$\mathrm{\AA}$', fontsize=font)
		plt.text(textpositionx, textpositiony-3*deltay, r'$\mathrm{\mu} =$'+str(dipolemoment)+"Debye", fontsize=font)
		plt.text(textpositionx, textpositiony-4*deltay, r'$\mathrm{\tau} =$' +str(parameter) +' '+"K"+r'$^{-1}$', fontsize=font)

	if (variableName == "tau" and Text2):
		#plt.text(textpositionx, textpositiony+5*deltay, 'Parameters:', fontsize=font)
		#plt.text(textpositionx, textpositiony+4*deltay, r'$\mathrm{N} =$'+str(numbmolecules)+" "+molecule, fontsize=font)
		#plt.text(textpositionx, textpositiony+3*deltay, r'$\mathrm{R} =$'+str(Rpt)+r'$\mathrm{\AA}$', fontsize=font)
		#plt.text(textpositionx, textpositiony+1*deltay, r'$\mathrm{\beta} =$' +str(parameter) +' '+"K"+r'$^{-1}$', fontsize=font)
	'''


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
				PlotLabel(Text1, Text2,fontlegend,xmin,xmax,ymin,ymax,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
			'''

		plt.subplots_adjust(top=0.95, bottom=0.15, left=0.09, right=0.98, hspace=0.0, wspace=0.4)
		plt.legend(bbox_to_anchor=(0.78, 0.75), loc=2, borderaxespad=1., shadow=True, fontsize=fontlegend)
		plt.savefig(outfile, dpi=200, format='eps')

		call(["open", outfile])
		#call(["okular", outfile])
		# plt.show()


def FigureAngleDistributionGfactor(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, TypePlot, beadsRef):
	font = 28
	fontlegend = font/2.0
	preskip = 0
	postskip = 0

	if (((TypePlot == "RFACTOR") or (TypePlot == "GFACTOR")) and variableName == "tau"):
		fig = plt.figure(figsize=(8, 6))
		plt.grid(True)

		var = beadsRef
		gFact = -1.0
		FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, gFact, dipolemoment, parameterName,
											   parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, var)
		FilePlotCorr = FilePlotName.SaveCorrGFAC+".eps"
		outfile = FilePlotCorr
		print("--------------------------------------------------")
		print("Name of the path and the Figure is given below - ")
		print(outfile)
		print("--------------------------------------------------")
		call(["rm", FilePlotCorr])
#
		colorList = ['red', 'green', 'blue', 'magenta']
		lsList = ['-', '--', '-.', '-']
		markerList = ['o', '^', 'v', 's']
		labelList = ['N = 24', 'N = 48', 'N = 64']
		nn = [24, 48]
		iLabel = 0
		gFactorList = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
		DList = np.zeros(len(gFactorList))
		ig = 0
		for gFactor in gFactorList:
			DipoleMoment = support.GetDipoleMomentFromGFactor(
				molecule, Rpt, gFactor)
			output = '{:1.4f}'.format(DipoleMoment)
			DList[ig] = output
			ig = ig+1

		# print(DList)
		for numbmolecules in nn:
			gFactorPlot = np.zeros(len(DList))
			CorrPlot = np.zeros(len(DList))
			err_CorrPlot = np.zeros(len(DList))

			iii = 0
			for dipolemoment in DList:
				gFact = -1.0
				FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, gFact, dipolemoment, parameterName,
													   parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, beadsRef)
				FileToBePlotCorr = FilePlotName.SaveCorr+".txt"
				if os.path.isfile(FileToBePlotCorr):
					beads1, var1, Corr, err_Corr = genfromtxt(FileToBePlotCorr, unpack=True, usecols=[
															  0, 1, 8, 15], skip_header=preskip, skip_footer=postskip)
					FactorList = support.GetrAndgFactor(
						molecule_rot, Rpt, dipolemoment)
					if (TypePlot == "GFACTOR"):
						gFactorPlot[iii] = FactorList[1]
					if (TypePlot == "RFACTOR"):
						rFactorPlot[iii] = FactorList[0]

					beadsRef1 = beadsRef
					# if ((beadsRef not in Corr) and (np.isscalar(Corr) == False)):
					#	beadsRef1 = beads1[-1]
					#	print(beadsRef1)

					if (np.isscalar(Corr) == True):
						CorrPlot[iii] = Corr
						err_CorrPlot[iii] = err_Corr
					else:
						ii = 0
						for i in beads1:
							indexi = int(i+0.5)
							beads = indexi
							if beads == beadsRef1:
								CorrPlot[iii] = Corr[ii]
								err_CorrPlot[iii] = err_Corr[ii]
							ii += 1
					iii += 1
					beadsRef1 = beadsRef

				# print("Corr")
				# print(CorrPlot)
#
			plt.errorbar(gFactorPlot,CorrPlot,yerr=err_CorrPlot,color=colorList[iLabel],ls=lsList[iLabel],linewidth=1,marker=markerList[iLabel],markersize=8,label=labelList[iLabel])
			iLabel += 1

			plt.xlim(0.4, 1.6)
			plt.ylim(0.0, 0.8)
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
				plt.xlabel(r'${g}$', fontsize=font)
			plt.ylabel(r'${\phi_{\mathrm{abs}}}$', fontsize=font)
			ymin, ymax = plt.ylim()
			if ymin < 0.0:
				plt.ylim(0.0, ymax)
			xmin, xmax = plt.xlim()

		plt.xticks(fontsize=font, rotation=0)
		plt.yticks(fontsize=font, rotation=0)

		plt.subplots_adjust(top=0.96, bottom=0.15, left=0.14, right=0.95, hspace=0.0, wspace=0.4)
		plt.legend(bbox_to_anchor=(0.75, 0.55), loc=2, borderaxespad=1., shadow=True, fontsize=fontlegend)
		plt.savefig(outfile, dpi=200, format='eps')

		#call(["open", outfile])
		plt.show()


def GetFigureEntropyRT_vs_gFactor(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,parameterName,parameter,numbblocks,numbpass,molecule,ENT_TYPE,preskip1,postskip1,extra_file_name,src_dir,TypePlot,beadsRef,plotNumber,purpose):
	'''
	The below 4 parameters are defined arbitrarily 
	as "support.GetFileNamePlot" function needs! 
	'''
	gfact=0.0
	dipolemoment=-1.0
	numbmolecules=2
	particleA=1
	#-------------------------------------------------#
	# Here some parameters are defined for the Figure!
	#-------------------------------------------------#
	font=24
	fontlegend=font/2.0
	fig=plt.figure(figsize=(8, 12))
	# plt.grid(True)
	colorList=['red', 'green', 'magenta', 'blue']
	lsList=['-', '--', '-.', '-']
	markerList=['o', '^', 'h', 's']
	#
	var = beadsRef
	FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip1,postskip1,extra_file_name,src_dir,particleA,var)
	FilePlotEntropy=FilePlotName.SaveEntropyGFACFit+"-"+purpose+".pdf"
	print(FilePlotEntropy)
	outfileEntropy=FilePlotEntropy
	call(["rm", FilePlotEntropy])
	#
	plotnumb = [1, 2]
	isubplot = 1
	for plotNumber in plotnumb:
		plt.subplot(2, 1, isubplot)
		isubplot = isubplot+1
		plt.axhline(y=math.log(2.0),color='green',lw=2.0,linestyle='-.',label='log(2)')
		#------------------------------------------------#
		nnDict = {1: [2, 4], 2: [8, 16]}
		#nnDict = {1:[2,4],2:[8,16]}
		FigLabel = {1: 'a', 2: 'b'}
		var = beadsRef
		#FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip1,postskip1,extra_file_name,src_dir,particleA,var)
		#FilePlotEntropy=FilePlotName.SaveEntropyGFACFit+"-"+FigLabel[plotNumber]+"-"+purpose+".pdf"
		#print(FilePlotEntropy)
		#outfileEntropy=FilePlotEntropy
		#call(["rm", FilePlotEntropy])
	#
		FileToBePlotDMRG=src_dir+"rotor_S2"
		iLabel=0
		nn = nnDict[plotNumber]
		for numbmolecules in nn:

			particleA = int(numbmolecules/2.0)

			if (numbmolecules == 2):
				#gFactorList  = [0.5+0.1*i for i in range(76)]
				gFactorList = [0.6+0.2*i for i in range(38)]
				gFactorList += [9.0+1.0*i for i in range(18)]
			if (numbmolecules == 4):
				#gFactorList  = [0.5+0.1*i for i in range(26)]
				gFactorList = [0.5+0.1*i for i in range(31)]
			if (numbmolecules == 8):
				gFactorList = [0.5+0.1*i for i in range(16)]
			if (numbmolecules == 16):
				gFactorList = [0.5+0.1*i for i in range(10)]
			if (numbmolecules == 32):
				gFactorList = [0.5+0.05*i for i in range(15)]

			gFactorPlot = gFactorList
			entropy1Plot = np.zeros(len(gFactorList))
			purity1Plot = np.zeros(len(gFactorList))
			err_entropy1Plot = np.zeros(len(gFactorList))
			err_purity1Plot = np.zeros(len(gFactorList))

			varEDPlot = np.zeros(len(gFactorList))
			EntropyEDPlot = np.zeros(len(gFactorList))
			EntropyFitPlot = np.zeros(len(gFactorList))
			ErrorFitPlot = np.zeros(len(gFactorList))
			#cutDataDMRG = {2:0, 4:27, 8:140, 16:195, 32:218}
			cutDataDMRG = {2: 0, 4: 77, 8: 140, 16: 195, 32: 218}

			iii = 0
			for gFact in gFactorList:
				gFact = '{:03.2f}'.format(gFact)
				gFact = float(gFact)
				FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,gFact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip1,postskip1,extra_file_name,src_dir,particleA,beadsRef)
				FileToBePlotEntropy=FilePlotName.SaveEntropyRT+".txt"

				beads1,var1,purity1,entropy1,err_purity1,err_entropy1=genfromtxt(FileToBePlotEntropy,unpack=True,usecols=[0, 1, 2, 3, 4, 5],skip_header=0,skip_footer=0)
				if (numbmolecules <= 32):
					EntropyFit, ErrorFit = GetFitPurity(FileToBePlotEntropy, numbmolecules, "g", gFact)
					EntropyFitPlot[iii] = -math.log(EntropyFit)
					ErrorFitPlot[iii] = abs(ErrorFit/EntropyFit)

				if (np.isscalar(entropy1) == True):
					entropy1Plot[iii] = entropy1
					err_entropy1Plot[iii] = err_entropy1
					purity1Plot[iii] = purity1
					err_purity1Plot[iii] = err_entropy1

					iii += 1
				else:
					ii = 0
					for i in beads1:
						indexi = int(i+0.5)
						beads = indexi
						if (beads == beadsRef):
							entropy1Plot[iii] = entropy1[ii]
							err_purity1Plot[iii] = err_purity1[ii]
							purity1Plot[iii] = purity1[ii]
							err_entropy1Plot[iii] = err_entropy1[ii]

						ii += 1
					iii += 1

			print("S2:	PIGS "+str(numbmolecules))
			print(entropy1Plot)
			print("S2:	PIGS Fit "+str(numbmolecules))
			print(EntropyFitPlot)

			if (numbmolecules == 64):
				labelString = 'PIGS:  '+r'$N = \ $'+str(numbmolecules)
				plt.errorbar(gFactorPlot, entropy1Plot, yerr=err_entropy1Plot, color=colorList[iLabel], ls='None', linewidth=0.5,  marker=markerList[iLabel], markersize=8, label=labelString)
			else:
				labelString = 'PIGS:  '+r'$N= \ $'+str(numbmolecules)
				plt.errorbar(gFactorPlot, EntropyFitPlot, yerr=ErrorFitPlot, color=colorList[iLabel], ls='None', linewidth=0.5,  marker=markerList[iLabel], markersize=8, label=labelString)

	# Data taken from Dmitri's DMRG
			labelStringDMRG = 'DMRG:	 '+r'$N = \ $'+str(numbmolecules)
			iRotors, rFact, EntropyFull = genfromtxt(FileToBePlotDMRG, unpack=True, usecols=[0, 1, 3])

			gFactDMRG = []
			EntropyDMRG = []
			for i in range(int(len(iRotors))):
				if (iRotors[i] == numbmolecules):
					gFactDMRG.append(1.0/(rFact[i]*rFact[i]*rFact[i]))
					EntropyDMRG.append(EntropyFull[i])

			plt.plot(gFactDMRG[cutDataDMRG[numbmolecules]:-110],EntropyDMRG[cutDataDMRG[numbmolecules]:-110],color='black',ls=lsList[iLabel],linewidth=2,marker='None',markersize=8,label=labelStringDMRG)
	#
			iLabel += 1

		plt.xlabel(r'$g$', fontsize=font, labelpad=0)
		plt.ylabel(r'$S_{2}$', fontsize=font, labelpad=8)

		if (nn[-1] >= 8):
			plt.ylim(-0.005, 0.901)
			Text1 = r'$\mathrm{(b)}$'
		else:
			plt.ylim(0.0, 0.901)
			Text1 = r'$\mathrm{(a)}$'
		xminlim = {1: 0.0, 2: 0.0}
		xmaxlim = {1: 20.01, 2: 20.01}
		yminlim = {1: 0.0, 2: 0.0}
		ymaxlim = {1: 0.901, 2: 0.901}
		stepx = {1: 1.0, 2: 0.25}
		stepy = {1: 0.1, 2: 0.1}
		plt.xticks(np.arange(xminlim[plotNumber], xmaxlim[plotNumber], step=stepx[plotNumber]), fontsize=font, rotation=0)
		plt.yticks(np.arange(yminlim[plotNumber], ymaxlim[plotNumber], step=stepy[plotNumber]), fontsize=font, rotation=0)

		xminlim = {1: 0.0, 2: 0.49}
		xmaxlim = {1: 8.01, 2: 2.03}
		plt.xlim(xminlim[plotNumber], xmaxlim[plotNumber])
		xmin, xmax = plt.xlim()
		ymin, ymax = plt.ylim()
		if Text1:
			plt.text((xmin+(xmax-xmin)*0.01), ymax-(ymax-ymin)*0.08, Text1, fontsize=font)
		plt.legend(bbox_to_anchor=(0.635, 0.50), loc=2, borderaxespad=1., shadow=True, fontsize=0.6*font)

	plt.rcParams["font.family"] = "serif"
	plt.rcParams["mathtext.fontset"] = "dejavuserif"
	plt.subplots_adjust(top=0.99, bottom=0.08, left=0.13, right=0.98, hspace=0.20, wspace=0.0)
	plt.savefig(outfileEntropy, dpi=50, format='pdf')
	plt.show()

def GetFigureEntropy_vs_beta(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, parameterName, parameter, numbblocks, numbpass, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, TypePlot, numbmolecules):
	'''
	The below 4 parameters are defined arbitrarily 
	as "support.GetFileNamePlot" function needs! 
	'''
	gfact = -1.0
	dipolemoment = -1.0
	particleA = 1
	#-------------------------------------------------#
	# Here some parameters are defined for the Figure!
	#-------------------------------------------------#
	fig = plt.figure(figsize=(8, 12))
	colorList = ['blue', 'green', 'red', 'magenta', 'black', 'cyan']
	lsList = ['-', '--', '-.', '-']
	markerList = ['o', '^', '*', 's', 'D']
	#------------------------------------------------#
	FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip1,postskip1,extra_file_name,src_dir,particleA,21)
	#------------------------------------------------#
	# 21 in the argument list of support.GetFileNamePlot
	# is just arbitrary number here
	#------------------------------------------------#

	FigLabel = {2: 'a', 8: 'b', 16: 'c', 32: 'd'}
	#FilePlotEntropy=FilePlotName.SaveEntropyRT+"-"+FigLabel[numbmolecules]+".pdf"
	FilePlotEntropy=FilePlotName.SaveEntropyRT+".pdf"
	outfileEntropy=FilePlotEntropy
	print(outfileEntropy)
	call(["rm", FilePlotEntropy])
#
	rotnumb = [2, 16]
	isubplot = 1
	for numbmolecules in rotnumb:
		if (numbmolecules == 2):
			#gFactorList  = [1.0, 2.0, 4.0, 6.0, 8.0]
			gFactorList = [1.0, 2.0, 4.0]
		if (numbmolecules == 8):
			#gFactorList  = [1.0, 1.5, 2.0]
			gFactorList = [1.0, 1.5]
		if (numbmolecules == 16):
			gFactorList = [1.0, 1.3]
		if (numbmolecules == 32):
			gFactorList = [1.0, 1.1]
	#
		iLabel = 0
		plt.subplot(2, 1, isubplot)
		isubplot = isubplot+1
		for gFact in gFactorList:
			gFact = '{:03.1f}'.format(gFact)
			gFact = float(gFact)
			FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, gFact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, 11)
			FileToBePlotEntropy = FilePlotName.SaveEntropyRT+".txt"
			print(FileToBePlotEntropy)

			beads1, var1, purity1, entropy1, err_purity1, err_entropy1 = genfromtxt(FileToBePlotEntropy, unpack=True, usecols=[0, 1, 2, 3, 4, 5], skip_header=0, skip_footer=0)

			print("S2:	PIGS "+str(gFact))
			print(entropy1)
	#
			labelString = r'$g = $'+str(gFact)
			plt.errorbar(var1,entropy1,yerr=err_entropy1,color=colorList[iLabel],ls=lsList[1],linewidth=1,marker=markerList[iLabel],markersize=10,label=labelString)
			iLabel += 1

		plt.xlabel(r'$\beta \ (K^{-1})$', fontsize=font, labelpad=5)
		plt.ylabel(r'$S_{2}$', fontsize=font, labelpad=5)

		xminlim = {0.005: 0.0, 0.02: 0.05, 0.04: 0.05}
		xmaxlim = {0.005: 0.2001, 0.02: 0.4001, 0.04: 0.4001}
		stepx = {0.005: 0.04, 0.02: 0.05, 0.04: 0.5}
		yminlim = {0.005: 0.0, 0.02: 0.3, 0.04: 0.3}
		ymaxlim = {0.005: 0.701, 0.02: 1.001, 0.04: 1.001}
		stepy = {0.005: 0.1, 0.02: 0.1, 0.04: 0.2}

		plt.xticks(np.arange(xminlim[parameter],xmaxlim[parameter],step=stepx[parameter]),fontsize=font,rotation=0)
		plt.yticks(np.arange(yminlim[parameter],ymaxlim[parameter],step=stepy[parameter]),fontsize=font,rotation=0)
		xmin, xmax = plt.xlim()
		ymin, ymax = plt.ylim()

		TextLabel = {2: r'$\mathrm{(a)}$', 8: r'$\mathrm{(b)}$', 16: r'$\mathrm{(b)}$', 32: r'$\mathrm{(d)}$'}

		if (numbmolecules == 2):
			plt.text((xmin+(xmax-xmin)*0.001), ymax-(ymax-ymin)*0.06, TextLabel[numbmolecules], fontsize=font)
			plt.text((xmin+(xmax-xmin)*0.4), ymax-(ymax-ymin)*0.15, r'$N =  $'+' '+str(numbmolecules), fontsize=font)
			plt.legend(bbox_to_anchor=(0.75, 0.75), loc=2, borderaxespad=1., shadow=True, fontsize=0.6*font)
			plt.xlim(-0.003,0.206)
			plt.ylim(0.0,0.72)
		if (numbmolecules == 8):
			plt.legend(bbox_to_anchor=(0.65, 0.50), loc=2, borderaxespad=1., shadow=True, fontsize=0.6*font)
		if (numbmolecules == 16):
			plt.text((xmin+(xmax-xmin)*0.001), ymax-(ymax-ymin)*0.1, TextLabel[numbmolecules], fontsize=font)
			plt.text((xmin+(xmax-xmin)*0.4), ymax-(ymax-ymin)*0.1, r'$N =  $'+' '+str(numbmolecules), fontsize=font)
			plt.legend(bbox_to_anchor=(0.75, 0.25), loc=2, borderaxespad=1., shadow=True, fontsize=0.6*font)
			plt.xlim(-0.003,0.206)
			plt.ylim(0.0,0.72)
		if (numbmolecules == 32):
			plt.legend(bbox_to_anchor=(0.65, 0.30), loc=2, borderaxespad=1., shadow=True, fontsize=0.6*font)
	plt.rcParams["font.family"] = "serif"
	plt.rcParams["mathtext.fontset"] = "dejavuserif"
	plt.subplots_adjust(top=0.99, bottom=0.08, left=0.12, right=0.99, hspace=0.25, wspace=0.0)
	plt.savefig(outfileEntropy, dpi=50, format='pdf')
	plt.show()

def fitFuncEntropy(var, a, b, c):
	return a+b*var*var+c*var*var*var


def fitFuncEntropy1(var, a, b):
	return a+b*var*var

# def plotfittingEntropy(val1, val2, val3):
#	markersize_fig = 10
#	xdata = np.linspace(0, np.max(val1), 1000)
#	fitParams, fitCovariances = curve_fit(fitFuncEntropy, val1, val2, sigma = val3)
#	plt.plot(xdata, fitFunc(xdata, fitParams[0], fitParams[1]), linestyle = '-', color = 'r', label = 'Fit', lw = 2)


def GetFitPurity(fname, numbmolecules, variableName, gFact):
	data = np.loadtxt(fname)
	listOfBeads = list(data[:, 0])
	cutoff = {2: 21, 4: 11, 8: 11, 16: 5, 32: 5}
	GetIndex = listOfBeads.index(cutoff[numbmolecules])
	xdata = data[GetIndex:, 1]
	ydata = data[GetIndex:, 2]
	ydata_err = data[GetIndex:, 4]
	if ((numbmolecules == 2) and (gFact < 5.0)):
		xdata = data[GetIndex+2:, 1]
		ydata = data[GetIndex+2:, 2]
		ydata_err = data[GetIndex+2:, 4]
	elif ((numbmolecules == 2) and (gFact >= 5.0) and (gFact < 7.0)):
		xdata = data[GetIndex+1:, 1]
		ydata = data[GetIndex+1:, 2]
		ydata_err = data[GetIndex+1:, 4]
	elif ((numbmolecules == 2) and (gFact >= 7.0) and (gFact < 9.0)):
		xdata = data[GetIndex:, 1]
		ydata = data[GetIndex:, 2]
		ydata_err = data[GetIndex:, 4]
	elif ((numbmolecules == 2) and (gFact >= 9.0) and (gFact < 13.0)):
		xdata = data[GetIndex-1:, 1]
		ydata = data[GetIndex-1:, 2]
		ydata_err = data[GetIndex-1:, 4]
	elif ((numbmolecules == 2) and (gFact >= 13.0)):
		xdata = data[GetIndex-2:, 1]
		ydata = data[GetIndex-2:, 2]
		ydata_err = data[GetIndex-2:, 4]

	if ((numbmolecules == 4) and (gFact < 2.2)):
		xdata = data[GetIndex+2:, 1]
		ydata = data[GetIndex+2:, 2]
		ydata_err = data[GetIndex+2:, 4]
	elif ((numbmolecules == 4) and (gFact >= 2.2) and (gFact < 3.0)):
		xdata = data[GetIndex+1:, 1]
		ydata = data[GetIndex+1:, 2]
		ydata_err = data[GetIndex+1:, 4]
	elif ((numbmolecules == 4) and (gFact >= 3.0) and (gFact < 3.5)):
		xdata = data[GetIndex:, 1]
		ydata = data[GetIndex:, 2]
		ydata_err = data[GetIndex:, 4]
	elif ((numbmolecules == 4) and (gFact >= 3.5) and (gFact < 5.0)):
		xdata = data[GetIndex-1:, 1]
		ydata = data[GetIndex-1:, 2]
		ydata_err = data[GetIndex-1:, 4]

	if ((numbmolecules == 8) and (gFact < 1.7)):
		xdata = data[GetIndex:, 1]
		ydata = data[GetIndex:, 2]
		ydata_err = data[GetIndex:, 4]
	elif ((numbmolecules == 8) and (gFact >= 1.7) and (gFact < 3.0)):
		xdata = data[GetIndex:, 1]
		ydata = data[GetIndex:, 2]
		ydata_err = data[GetIndex:, 4]

	if ((numbmolecules == 16) and ((gFact < 1.1) or (gFact == 1.4))):
		xdata = data[GetIndex+1:, 1]
		ydata = data[GetIndex+1:, 2]
		ydata_err = data[GetIndex+1:, 4]

	if ((numbmolecules == 32) and (gFact < 1.0)):
		xdata = data[GetIndex+1:, 1]
		ydata = data[GetIndex+1:, 2]
		ydata_err = data[GetIndex+1:, 4]

	popt, pcov = curve_fit(fitFuncEntropy1, xdata, ydata, sigma=ydata_err)
	perr = np.sqrt(np.diag(pcov))

	t = np.linspace(0, np.max(xdata), 10)
	fit = fitFuncEntropy1(t, *popt)
	if (variableName == "tau"):
		return t, fit, perr[0]
	else:
		return fit[0], perr[0]


def GetFigureEntropyRT_vs_tau(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, gfact, parameterName, parameter, numbblocks, numbpass, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, TypePlot, numbmolecules):
	'''
	The below 4 parameters are defined arbitrarily 
	as "support.GetFileNamePlot" function needs! 
	'''
	dipolemoment = -1.0
	particleA = 1
	#-------------------------------------------------#
	# Here some parameters are defined for the Figure!
	#-------------------------------------------------#
	font = 24
	fontlegend = font/2.0
	fig = plt.figure(figsize=(8, 12))
	colorList = ['blue', 'green', 'red', 'magenta', 'black', 'cyan']
	lsList = ['-', '--', '-.', '-']
	markerList = ['o', '^', '*', 's', 'D']
	#------------------------------------------------#
	gfact = '{:03.1f}'.format(gfact)
	gfact = float(gfact)
#
	FigLabel = {2: 'a', 4: 'b', 8: 'c', 16: 'c', 32: 'e'}
	FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip1,postskip1,extra_file_name,src_dir,particleA,21)
	FilePlotEntropy = FilePlotName.SaveEntropyRT+"-"+FigLabel[numbmolecules]+".pdf"
	outfileEntropy = FilePlotEntropy
	print(outfileEntropy)
	call(["rm", FilePlotEntropy])
#
	isubplot = 1
	#rotnumb = [2, 4]
	#for numbmolecules in rotnumb:
	gfactnumb = [0.6, 1.4]
	for gfact in gfactnumb:
		FilePlotName = support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip1,postskip1,extra_file_name,src_dir,particleA,11)
		FileToBePlotEntropy = FilePlotName.SaveEntropyRT+".txt"
		print(FileToBePlotEntropy)

		beads1, var1, purity1, entropy1, err_purity1, err_entropy1 = genfromtxt(FileToBePlotEntropy, unpack=True, usecols=[0, 1, 2, 3, 4, 5], skip_header=1, skip_footer=0)

		print("S2:	PIGS "+str(gfact))
		print(entropy1)
	#
		plt.subplot(2, 1, isubplot)
		labelString = 'PIGS'
		plt.errorbar(var1,entropy1,yerr=err_entropy1,color='blue',ls='None',linewidth=1,marker='o',markersize=8,label=labelString)

		rFactor = 1.0/(pow(gfact, 1.0/3.0))
		rfact = round(rFactor, 3)
		print('rfactor = '+str(rfact))
	# Data taken from Dmitri's DMRG
		NDMRG, RDMRG, S2DMRG = genfromtxt(src_dir+'rotor_S2', unpack=True, usecols=[0, 1, 3], skip_header=1, skip_footer=0)
		RDMRGExtract = RDMRG[(NDMRG == numbmolecules)]
		S2DMRGExtract = S2DMRG[(NDMRG == numbmolecules)]
		S2DMRGPlot = S2DMRGExtract[(RDMRGExtract == rfact)]
		plt.axhline(y=S2DMRGPlot, color='green', lw=2.0, linestyle='-.', label='DMRG')
	# DMRG finished here

		tvar, purityFit, ErrorFit1 = GetFitPurity(FileToBePlotEntropy, numbmolecules, "tau", gfact)
		ErrorFit = np.zeros(int(len(purityFit)))
		ErrorFit[0] = ErrorFit1
		print("Error  "+str(ErrorFit1))

		EntropyFitPlot = -np.log(purityFit)
		ErrorFitPlot = np.fabs(np.divide(ErrorFit, purityFit))
		labelString = 'PIGS-FIT '
		plt.errorbar(tvar, EntropyFitPlot, yerr=ErrorFitPlot, color='red', ls='-', linewidth=2,  marker='None', markersize=20)
		plt.errorbar(tvar[0], EntropyFitPlot[0], yerr=ErrorFitPlot[0], color='red', ls='-', linewidth=2,  marker='s', markersize=8, label=labelString)

		plt.xlabel(r'$\tau \ (K^{-1})$', fontsize=font, labelpad=5)
		plt.ylabel(r'$S_{2}$', fontsize=font, labelpad=5)
	#
		if (numbmolecules < 16):
			xmaxlim = {1: 0.021, 2: 0.021}
			xminlim = {1: 0.0,   2: 0.0}
			ymaxlim = {1: 0.165, 2: 0.8}
			yminlim = {1: 0.0,   2: 0.0}
			TextLabel = {1: r'$\mathrm{(a)}$', 2: r'$\mathrm{(b)}$'}
			stepx = {1: 0.002, 2: 0.002}
			stepy = {1: 0.01, 2: 0.01}
			xboxsize = {1: 0.0, 2: 0.0}
			yboxsize = {1: 0.7, 2: 0.7}
		if (numbmolecules == 16):
			xmaxlim = {1: 0.021, 2: 0.021}
			xminlim = {1: 0.0, 2: 0.0}
			ymaxlim = {1: 2.0, 2: 1.0}
			yminlim = {1: 0.0, 2: 0.0}
			TextLabel = {1: r'$\mathrm{(a)}$', 2: r'$\mathrm{(b)}$'}
			stepx = {1: 0.004, 2: 0.004}
			stepy = {1: 0.02, 2: 0.05}
			xboxsize = {1: 0.7, 2: 0.7}
			yboxsize = {1: 0.55, 2: 0.3}

		plt.xticks(np.arange(xminlim[isubplot],xmaxlim[isubplot],step=stepx[isubplot]),fontsize=font,rotation=0)
		plt.xlim(xminlim[isubplot], xmaxlim[isubplot])
		plt.yticks(np.arange(yminlim[isubplot],ymaxlim[isubplot],step=stepy[isubplot]),fontsize=font,rotation=0)
		plt.ylim(yminlim[isubplot], ymaxlim[isubplot])
	#
		if (numbmolecules < 16):
			xmaxlim = {1: 0.00701, 2: 0.0103}
			xminlim = {1: -0.0001, 2: -0.0001}
			ymaxlim = {1: 0.0921, 2: 0.16}
			yminlim = {1: 0.0698, 2: 0.115}
		if (numbmolecules == 16):
			xmaxlim = {1: 0.021, 2: 0.021}
			xminlim = {1: -0.001, 2: -0.001}
			ymaxlim = {1: 0.09, 2: 0.85}
			yminlim = {1: 0.02, 2: 0.6}
		plt.xlim(xminlim[isubplot], xmaxlim[isubplot])
		plt.ylim(yminlim[isubplot], ymaxlim[isubplot])
		ymin, ymax = plt.ylim()
		xmin, xmax = plt.xlim()
		plt.text((xmin+(xmax-xmin)*0.35),ymax-(ymax-ymin)*0.07,r'$N = \ $'+str(numbmolecules)+';'+r'$ \ g = \ $'+str(gfact),fontsize=font)
		plt.text((xmin+(xmax-xmin)*0.01),ymax-(ymax-ymin)*0.07,TextLabel[isubplot],fontsize=font)
		plt.legend(bbox_to_anchor=(xboxsize[isubplot], yboxsize[isubplot]), loc=2, borderaxespad=1., shadow=True, fontsize=0.6*font)
		isubplot = isubplot+1

	plt.subplots_adjust(top=0.99, bottom=0.08, left=0.14, right=0.98, hspace=0.22, wspace=0.0)
	plt.rcParams["font.family"] = "serif"
	plt.rcParams["mathtext.fontset"] = "dejavuserif"
	plt.savefig(outfileEntropy, dpi=50, format='pdf')
	plt.show()

def GetFigureEntropyRT_vs_gFactor_COMBO(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,parameterName,parameter,numbblocks,numbpass,molecule,ENT_TYPE_LIST,preskip1,postskip1,extra_file_name,src_dir,TypePlot,beadsRef,numbmolecules,purpose):
	'''
	The below 4 parameters are defined arbitrarily 
	as "support.GetFileNamePlot" function needs! 
	'''
	dipolemoment = -1.0
	particleA = int(numbmolecules/2)
	gfact = -1.0
	#-------------------------------------------------#
	# Here some parameters are defined for the Figure!
	#-------------------------------------------------#
	if (purpose == "article"):
		font = 24
	elif (purpose == "ppt"):
		font = 40
		rc('axes', linewidth=2)
	fontlegend = font/2.0
	fig = plt.figure(figsize=(8, 18))
	# plt.grid(True)
	colorList = ['blue', 'red', 'magenta', 'black', 'cyan']
	lsList = ['-', '--', '-.', '-']
	markerList = ['^', 'o', '*', 's', 'D']
	#------------------------------------------------#
	FigLabel = {2: 'a', 4: 'b', 8: 'c'}
	FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,"",preskip1,postskip1,extra_file_name,src_dir,particleA,beadsRef)
	#FilePlotEntropy=FilePlotName.SaveEntropyCOMBO+"-"+FigLabel[numbmolecules]+"-"+purpose+".pdf"
	FilePlotEntropy=FilePlotName.SaveEntropyCOMBO+"-"+purpose+".pdf"
	outfileEntropy=FilePlotEntropy
	print(outfileEntropy)
	call(["rm", FilePlotEntropy])
#
	rotnumb = [2, 4, 8]
	isubplot = 1
	for numbmolecules in rotnumb:
		if (numbmolecules == 2):
			#gFactorList  = [0.5+0.1*i for i in range(76)]
			gFactorList = [0.6+0.2*i for i in range(38)]
			gFactorList += [9.0+1.0*i for i in range(7)]
			gFactorList += [16, 17, 18, 19, 20]
		if (numbmolecules == 4):
			gFactorList = [0.5+0.1*i for i in range(31)]
			#gFactorList += [3.75+0.25*i for i in range(9)]
		if (numbmolecules == 8):
			gFactorList = [0.5+0.1*i for i in range(16)]
		if (numbmolecules == 16):
			gFactorList = [0.5+0.1*i for i in range(11)]
		if (numbmolecules == 32):
			gFactorList = [0.5+0.1*i for i in range(9)]

		entropy1Plot = np.zeros((len(ENT_TYPE_LIST), len(gFactorList)))
		purity1Plot = np.zeros((len(ENT_TYPE_LIST), len(gFactorList)))
		err_entropy1Plot = np.zeros((len(ENT_TYPE_LIST), len(gFactorList)))
		err_purity1Plot = np.zeros((len(ENT_TYPE_LIST), len(gFactorList)))
		entropyFitPlot = np.zeros((len(ENT_TYPE_LIST), len(gFactorList)))
		errorFitPlot = np.zeros((len(ENT_TYPE_LIST), len(gFactorList)))

		FileToBePlotDMRG = src_dir+"rotor_S2"
		iLabel = 0
		for ENT_TYPE in ENT_TYPE_LIST:
			iii = 0
			for gfact in gFactorList:
				gfact = '{:03.2f}'.format(gfact)
				gfact = float(gfact)

				if (ENT_TYPE == "BROKENPATH"):
					preskipDict = {2: 10000, 4: 30000, 8: 0}
					if (gfact > 8.99):
						preskipDict = {2: 90000, 4: 30000, 8: 0}
				else:
					preskipDict = {2: 10000, 4: 10000, 8: 10000}
					if (gfact > 8.99):
						preskipDict = {2: 90000, 4: 30000, 8: 0}

				if ((ENT_TYPE == "BROKENPATH") and (numbmolecules > 2)):
					numbblocks = 50000
				else:
					numbblocks = 20000

				if (gfact > 8.99):
					numbblocks = 100000
				preskip1 = preskipDict[numbmolecules]
				# print(preskip1)
				FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip1,postskip1,extra_file_name,src_dir,particleA,beadsRef)
				FileToBePlotEntropy = FilePlotName.SaveEntropyRT+".txt"

				beads1, var1, purity1, entropy1, err_purity1, err_entropy1 = genfromtxt(FileToBePlotEntropy,unpack=True,usecols=[0, 1, 2, 3, 4, 5],skip_header=0,skip_footer=0)

				entropyFit, errorFit = GetFitPurity(FileToBePlotEntropy,numbmolecules,"g",gfact)
				entropyFitPlot[iLabel, iii] = -math.log(entropyFit)
				errorFitPlot[iLabel, iii] = abs(errorFit/entropyFit)

				if (np.isscalar(entropy1) == True):
					entropy1Plot[iLabel, iii] = entropy1
					err_entropy1Plot[iLabel, iii] = err_entropy1
					purity1Plot[iLabel, iii] = purity1
					err_purity1Plot[iLabel, iii] = err_entropy1

					iii += 1
				else:
					ii = 0
					for i in beads1:
						indexi = int(i+0.5)
						beads = indexi
						if (beads == beadsRef):
							entropy1Plot[iLabel, iii] = entropy1[ii]
							err_purity1Plot[iLabel, iii] = err_purity1[ii]
							purity1Plot[iLabel, iii] = purity1[ii]
							err_entropy1Plot[iLabel, iii] = err_entropy1[ii]

						ii += 1
					iii += 1
			iLabel += 1
	#
		iLabel = 0
		plt.subplot(3, 1, isubplot)
		isubplot = isubplot+1
		gFactorPlot = gFactorList
		if (purpose == "article"):
			textLabel = {'SWAPTOUNSWAP': 'Extended ensemble', 'BROKENPATH': 'Broken path ensemble'}
		elif (purpose == "ppt"):
			textLabel = {'SWAPTOUNSWAP': 'Extended ensemble', 'BROKENPATH': 'Broken path'}
		for ENT_TYPE in ENT_TYPE_LIST:
			labelString = textLabel[ENT_TYPE]
			#plt.errorbar(gFactorPlot, entropy1Plot[iLabel,:], yerr=err_entropy1Plot[iLabel,:], color = colorList[iLabel], ls = 'None', linewidth=1,  marker = markerList[iLabel], markersize = 8, label = labelString)
			plt.errorbar(gFactorPlot,entropyFitPlot[iLabel, :],yerr=errorFitPlot[iLabel, :],color=colorList[iLabel],ls='None',linewidth=1,marker=markerList[iLabel],markersize=8,label=labelString)
			iLabel += 1

		'''
	#Exact results obtained either by diagonalizing the Full Hamiltonian matrix or by DMRG #
		varEDPlot		 = np.zeros(len(gFactorList))
		entropyEDPlot	 = np.zeros(len(gFactorList))

		iii = 0
		for gfact in gFactorList:	
			gfact = '{:03.1f}'.format(gfact)
			gfact = float(gfact)

			FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, "SWAPTOUNSWAP", preskip1, postskip1, extra_file_name, src_dir, particleA, beadsRef)
			FileToBePlotEntropyED	  = FilePlotName.SaveEntropyED+".txt"

			if (numbmolecules <= 4):
				varED, entropyED = genfromtxt(FileToBePlotEntropyED,unpack=True, usecols=[0, 2])
				varEDPlot[iii]		  = 1.0/math.pow(varED,3)
				entropyEDPlot[iii]	  = entropyED
			iii += 1

		labelStringED = 'DMRG'

		#print(varEDPlot)
		#print(entropyEDPlot)
		plt.plot(varEDPlot, entropyEDPlot, color = 'black', ls = '-', linewidth=2,	marker = 'None', markersize = 8, label = labelStringED)
		'''
	# Plot of DMRG data
		FileToBePlotDMRG = src_dir+"rotor_S2"
		labelStringDMRG = 'DMRG'
		iRotors, rFact, EntropyFull = genfromtxt(FileToBePlotDMRG,unpack=True,usecols=[0, 1, 3])

		gFactDMRG = []
		EntropyDMRG = []
		cutDataDMRG = {2: 0, 4: 77, 8: 140}
		#cutDataDMRG = {2:0, 4:0, 8:0}
		for i in range(int(len(iRotors))):
			if (iRotors[i] == numbmolecules):
				gFactDMRG.append(1.0/(rFact[i]*rFact[i]*rFact[i]))
				EntropyDMRG.append(EntropyFull[i])

		plt.plot(gFactDMRG[cutDataDMRG[numbmolecules]:-110],EntropyDMRG[cutDataDMRG[numbmolecules]:-110],color='black',ls='-',linewidth=2,marker='None',markersize=8,label=labelStringDMRG)
	# DMRG section ends here
		xmaxlim = {2: 800.01, 4: 800.01, 8: 200.01}
		xminlim = {2: 0.0, 4: 0.0, 8: 0.0}
		ymaxlim = {2: 9, 4: 8, 8: 8.0}
		yminlim = {2: 0.0, 4: 0.0, 8: 0.0}

		if (purpose == "article"):
			plt.xlabel(r'$g$', fontsize=font, labelpad=0)
			plt.ylabel(r'$S_{2}$', fontsize=font, labelpad=8)
			TextLabel = {2:r'$\mathrm{(a)}$', 4:r'$\mathrm{(b)}$', 8:r'$\mathrm{(c)}$', 16:r'$\mathrm{(d)}$', 32:r'$\mathrm{(e)}$'}
			xboxsize = {2: 0.56, 4: 0.56, 8: 0.56}
			yboxsize = {2: 0.3, 4: 0.3, 8: 0.3}
			stepx = {2: 2.0, 4: 0.5, 8: 0.4}
			stepy = {2: 0.1, 4: 0.1, 8: 0.1}
			plt.xticks(np.arange(xminlim[numbmolecules],xmaxlim[numbmolecules],step=stepx[numbmolecules]),fontsize=font,rotation=0)
			plt.yticks(np.arange(yminlim[numbmolecules],ymaxlim[numbmolecules],step=stepy[numbmolecules]),fontsize=font,rotation=0)
			xmaxlim = {2: 8.05, 4: 3.55, 8: 2.01}
			xminlim = {2: -0.02, 4: -0.02, 8: 0.38}
			ymaxlim = {2: 0.901, 4: 0.901, 8: 0.901}
			yminlim = {2: 0.0, 4: 0.0, 8: 0.0}
		elif (purpose == "ppt"):
			plt.xlabel(r'$g$', fontsize=font, labelpad=-8)
			plt.ylabel(r'$S_{2}$', fontsize=font, labelpad=2)
			TextLabel = {2: '', 4: '', 8: '', 16: '', 32: ''}
			stepx = {2: 1.0, 4: 0.5, 8: 0.5}
			stepy = {2: 0.2, 4: 0.2, 8: 0.2}
			plt.xticks(np.arange(xminlim[numbmolecules],xmaxlim[numbmolecules],step=stepx[numbmolecules]),fontsize=font,rotation=0)
			plt.yticks(np.arange(yminlim[numbmolecules],ymaxlim[numbmolecules],step=stepy[numbmolecules]),fontsize=font,rotation=0)
			xboxsize = {2: 0.5, 4: 0.5, 8: 0.5}
			yboxsize = {2: 0.4, 4: 0.4, 8: 0.4}
			xmaxlim = {2: 8.01, 4: 3.01, 8: 2.01}
			xminlim = {2: 0.0, 4: 0.4, 8: 0.4}
			ymaxlim = {2: 0.84, 4: 0.84, 8: 0.84}
			yminlim = {2: 0.0, 4: 0.0, 8: 0.0}

		plt.xlim(xminlim[numbmolecules], xmaxlim[numbmolecules])
		plt.ylim(yminlim[numbmolecules], ymaxlim[numbmolecules])
		ymin, ymax = plt.ylim()
		xmin, xmax = plt.xlim()

		plt.rcParams["font.family"] = "serif"
		plt.rcParams["mathtext.fontset"] = "dejavuserif"
		plt.text((xmin+(xmax-xmin)*0.01), ymax-(ymax-ymin)*0.08, TextLabel[numbmolecules], fontsize=font)
		plt.text((xmin+(xmax-xmin)*0.4), ymax-(ymax-ymin)*0.1, r'$N = \ $'+str(numbmolecules), fontsize=font)
		plt.legend(bbox_to_anchor=(xboxsize[numbmolecules], yboxsize[numbmolecules]), loc=2, borderaxespad=1., shadow=True, fontsize=0.6*font)

	if (purpose == "article"):
		plt.subplots_adjust(top=0.99, bottom=0.05, left=0.13, right=0.98, hspace=0.20, wspace=0.0)
	elif (purpose == "ppt"):
		plt.subplots_adjust(top=0.97, bottom=0.18, left=0.18, right=0.96, hspace=0.0, wspace=0.0)
		plt.legend(bbox_to_anchor=(xboxsize[numbmolecules], yboxsize[numbmolecules]), loc=2, borderaxespad=1., shadow=True, fontsize=fontlegend)
	plt.savefig(outfileEntropy, dpi=50, format='pdf')
	plt.show()

def GetFigureEntropyRT_vs_gFactor_COMPARE(TypeCal, ENT_TYPE, molecule_rot, TransMove, RotMove, variableName, Rpt, parameterName, parameter, numbblocks, numbpass, molecule, algorithm_LIST, preskip1, postskip1, extra_file_name, src_dir, TypePlot, beadsRef, numbmolecules, purpose):
	'''
	The below 4 parameters are defined arbitrarily 
	as "support.GetFileNamePlot" function needs! 
	'''
	dipolemoment = -1.0
	particleA = int(numbmolecules/2)
	gfact = -1.0
	#-------------------------------------------------#
	# Here some parameters are defined for the Figure!
	#-------------------------------------------------#
	if (purpose == "article"):
		font = 24
	elif (purpose == "ppt"):
		font = 40
		rc('axes', linewidth=2)
	fontlegend = font/2.0
	fig = plt.figure(figsize=(8, 18))
	# plt.grid(True)
	colorList = ['blue', 'red', 'magenta', 'black', 'cyan']
	lsList = ['-', '--', '-.', '-']
	markerList = ['^', 'o', '*', 's', 'D']
	#------------------------------------------------#
	FigLabel = {4: 'a', 8: 'b', 16: 'c'}
	FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules, molecule,"",preskip1,postskip1,extra_file_name,src_dir,particleA,beadsRef)
	#FilePlotEntropy = FilePlotName.SaveEntropyCOMPARE+"-"+FigLabel[numbmolecules]+"-"+purpose+".pdf"
	FilePlotEntropy = FilePlotName.SaveEntropyCOMPARE+"-"+purpose+".pdf"
	outfileEntropy = FilePlotEntropy
	print(outfileEntropy)
	call(["rm", FilePlotEntropy])
#
	rotnumb = [4, 8, 16]
	isubplot = 1
	for numbmolecules in rotnumb:
		if (numbmolecules == 2):
			#gFactorList  = [0.5+0.1*i for i in range(76)]
			gFactorList = [0.6+0.2*i for i in range(38)]
			gFactorList += [9.0+1.0*i for i in range(20)]
		if (numbmolecules == 4):
			#gFactorList  = [0.5+0.1*i for i in range(26)]
			gFactorList = [0.5+0.1*i for i in range(31)]
		if (numbmolecules == 8):
			gFactorList = [0.5+0.1*i for i in range(16)]
		if (numbmolecules == 16):
			gFactorList = [0.5+0.1*i for i in range(10)]
		if (numbmolecules == 32):
			gFactorList = [0.5+0.1*i for i in range(9)]
		print(gFactorList)

		entropy1Plot = np.zeros((len(algorithm_LIST), len(gFactorList)))
		purity1Plot = np.zeros((len(algorithm_LIST), len(gFactorList)))
		err_entropy1Plot = np.zeros((len(algorithm_LIST), len(gFactorList)))
		err_purity1Plot = np.zeros((len(algorithm_LIST), len(gFactorList)))
		entropyFitPlot = np.zeros((len(algorithm_LIST), len(gFactorList)))
		errorFitPlot = np.zeros((len(algorithm_LIST), len(gFactorList)))

		FileToBePlotDMRG = src_dir+"rotor_S2"
		iLabel = 0
		for algorithm in algorithm_LIST:
			print(algorithm)
			iii = 0
			for gfact in gFactorList:
				gfact = '{:03.1f}'.format(gfact)
				gfact = float(gfact)

				if (algorithm == algorithm_LIST[1]):
					preskipDict = {4: 0, 8: 0, 16: 0}
				else:
					preskipDict = {4: 10000, 8: 10000, 16: 10000}

				if ((algorithm == algorithm_LIST[1]) and (numbmolecules == 4)):
					numbblocks = 30000
				elif ((algorithm == algorithm_LIST[1]) and (numbmolecules == 8)):
					numbblocks = 50000
				elif ((algorithm == algorithm_LIST[1]) and (numbmolecules == 16)):
					numbblocks = 80000
				else:
					numbblocks = 20000
				preskip1 = preskipDict[numbmolecules]
				# print(preskip1)
				if (algorithm == algorithm_LIST[1]):
					extra = algorithm_LIST[1]
				else:
					extra = ""
				FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip1,postskip1,extra,src_dir,particleA,beadsRef)
				FileToBePlotEntropy = FilePlotName.SaveEntropyRT+".txt"

				beads1, var1, purity1, entropy1, err_purity1, err_entropy1 = genfromtxt(FileToBePlotEntropy, unpack=True, usecols=[0, 1, 2, 3, 4, 5], skip_header=0, skip_footer=0)

				entropyFit, errorFit = GetFitPurity(FileToBePlotEntropy, numbmolecules, "g", gfact)
				entropyFitPlot[iLabel, iii] = -math.log(entropyFit)
				errorFitPlot[iLabel, iii] = abs(errorFit/entropyFit)
				iii = iii+1

			iLabel += 1
	#
		iLabel = 0
		plt.subplot(3, 1, isubplot)
		isubplot = isubplot+1
		gFactorPlot = gFactorList

		textLabel = {0: '$\mathrm{with \ ratio \ trick}$', 1: '$\mathrm{without \ ratio \ trick}$'}

		for algorithm in algorithm_LIST:
			labelString = textLabel[iLabel]
			#plt.errorbar(gFactorPlot, entropy1Plot[iLabel,:], yerr=err_entropy1Plot[iLabel,:], color = colorList[iLabel], ls = 'None', linewidth=1,  marker = markerList[iLabel], markersize = 8, label = labelString)
			plt.errorbar(gFactorPlot,entropyFitPlot[iLabel, :],yerr=errorFitPlot[iLabel, :],color=colorList[iLabel],ls='None',linewidth=1,marker=markerList[iLabel],markersize=8,label=labelString)
			iLabel += 1

	# Plot of DMRG data
		FileToBePlotDMRG = src_dir+"rotor_S2"
		labelStringDMRG = 'DMRG'
		iRotors, rFact, EntropyFull = genfromtxt(FileToBePlotDMRG, unpack=True, usecols=[0, 1, 3])

		gFactDMRG = []
		EntropyDMRG = []
		cutDataDMRG = {4: 77, 8: 140, 16: 195}
		for i in range(int(len(iRotors))):
			if (iRotors[i] == numbmolecules):
				gFactDMRG.append(1.0/(rFact[i]*rFact[i]*rFact[i]))
				EntropyDMRG.append(EntropyFull[i])

		plt.plot(gFactDMRG[cutDataDMRG[numbmolecules]:-110],EntropyDMRG[cutDataDMRG[numbmolecules]:-110],color='black',ls='-',linewidth=2,marker='None',markersize=8,label=labelStringDMRG)
	# DMRG section ends here
		xmaxlim = {4: 4.01, 8: 2.01, 16: 1.501}
		xminlim = {4: 0.0, 8: 0.0, 16: 0.0}
		ymaxlim = {4: 2.0, 8: 2.0, 16: 2.0}
		yminlim = {4: 0.0, 8: 0.0, 16: 0.0}

		if (purpose == "article"):
			plt.xlabel(r'$g$', fontsize=font, labelpad=0)
			plt.ylabel(r'$S_{2}$', fontsize=font, labelpad=8)
			TextLabel={4:r'$\mathrm{(a)}$', 8:r'$\mathrm{(b)}$', 16:r'$\mathrm{(c)}$'}
			stepx = {4: 0.5, 8: 0.4, 16: 0.25}
			stepy = {4: 0.1, 8: 0.1, 16: .1}
			plt.xticks(np.arange(xminlim[numbmolecules],xmaxlim[numbmolecules],step=stepx[numbmolecules]),fontsize=font,rotation=0)
			plt.yticks(np.arange(yminlim[numbmolecules],ymaxlim[numbmolecules],step=stepy[numbmolecules]),fontsize=font,rotation=0)
			xboxsize = {4: 0.6, 8: 0.6, 16: 0.6}
			yboxsize = {4: 0.5, 8: 0.5, 16: 0.5}
			xmaxlim = {4: 3.55, 8: 2.01, 16: 1.52}
			xminlim = {4: -0.02, 8: 0.38, 16: 0.48}
			ymaxlim = {4: 0.901, 8: 0.901, 16: 0.901}
			yminlim = {4: 0.0, 8: 0.0, 16: 0.0}
		elif (purpose == "ppt"):
			plt.xlabel(r'$g$', fontsize=font, labelpad=-8)
			plt.ylabel(r'$S_{2}$', fontsize=font, labelpad=2)
			TextLabel = {4: '', 8: '', 16: ''}
			stepx = {4: 0.5, 8: 0.5, 16: 0.25}
			stepy = {4: 0.2, 8: 0.2, 16: .2}
			plt.xticks(np.arange(xminlim[numbmolecules],xmaxlim[numbmolecules],step=stepx[numbmolecules]),fontsize=font,rotation=0)
			plt.yticks(np.arange(yminlim[numbmolecules],ymaxlim[numbmolecules],step=stepy[numbmolecules]),fontsize=font,rotation=0)
			xboxsize = {4: 0.45,  8: 0.45,	 16: 0.45}
			yboxsize = {4: 0.35,  8: 0.35,	 16: 0.35}
			xmaxlim = {4: 3.01,  8: 2.0001, 16: 1.53}
			xminlim = {4: 0.4,  8: 0.4,	16: 0.4}
			ymaxlim = {4: 0.82,  8: 0.82,	16: 0.82}
			yminlim = {4: 0.0,	 8: 0.0,	16: 0.0}
		plt.xlim(xminlim[numbmolecules], xmaxlim[numbmolecules])
		plt.ylim(yminlim[numbmolecules], ymaxlim[numbmolecules])
		ymin, ymax = plt.ylim()
		xmin, xmax = plt.xlim()

		plt.text((xmin+(xmax-xmin)*0.01), ymax-(ymax-ymin)*0.08, TextLabel[numbmolecules], fontsize=font)
		plt.text((xmin+(xmax-xmin)*0.4), ymax-(ymax-ymin)*0.1, r'$N = \ $'+str(numbmolecules), fontsize=font)
		plt.legend(bbox_to_anchor=(xboxsize[numbmolecules],yboxsize[numbmolecules]),loc=2,borderaxespad=1.,shadow=True,fontsize=0.6*font)

	if (purpose == "article"):
		#plt.subplots_adjust(top=0.97, bottom=0.12, left=0.13, right=0.97, hspace=0.2, wspace=0.0)
		plt.subplots_adjust(top=0.99, bottom=0.05, left=0.13, right=0.98, hspace=0.20, wspace=0.0)
	if (purpose == "ppt"):
		plt.subplots_adjust(top=0.97, bottom=0.18, left=0.18, right=0.96, hspace=0.0, wspace=0.0)
		plt.legend(bbox_to_anchor=(xboxsize[numbmolecules],yboxsize[numbmolecules]),loc=2,borderaxespad=1.,shadow=True,fontsize=0.6*font)
	plt.rcParams["font.family"] = "serif"
	plt.savefig(outfileEntropy, dpi=50, format='pdf')
	plt.show()

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

def GetFigureEnergyRot_vs_beta(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,preskip,postskip,extra_file_name,src_dir,TypePlot,purpose):
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
		beads_skip_header_final, beads_skip_footer_final = GetKeyIndices(numbmolecules,molecule_rot,Rpt,FileToBePlotEnergy)	

		beads_vec1, valTau1, valTotalEnergy1, errorTotalEnergy1 = np.genfromtxt(FileToBePlotEnergy, unpack=True, usecols=[1, 2, 4, 6], skip_header=0, skip_footer=0)
		trunc=beads_skip_header_final
		trunce=int(len(beads_vec1))-beads_skip_footer_final
		print(beads_vec1[trunc:trunce])
		valTau = valTau1[trunc:trunce]
		valTotalEnergy = valTotalEnergy1[trunc:trunce]/kcalmolinvKinv
		errorTotalEnergy = errorTotalEnergy1[trunc:trunce]/kcalmolinvKinv

		# Fitting of PIGS data
		fitting_term = "quatric" 
		tvar, EnergyFitPlot, ErrorFitPlot = GetFitEnergy(valTau, valTotalEnergy, errorTotalEnergy, variableName, fitting_term)
		energy_fit[i]=EnergyFitPlot[0]
		energy_fit_error[i]=ErrorFitPlot
		print(Rpt,'  ',energy_fit[i],'  ',energy_fit_error[i])
		i+=1

	#Plot begins
	FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,1.0,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip,postskip,extra_file_name,src_dir,particleA,11)
	FilePlotEnergy = FilePlotName.SaveEnergyFitvsR
	print(FilePlotEnergy+"-tot.pdf")
	call(["rm", FilePlotEnergy+"-tot.pdf"])

	# Plotting of fitting data
	plt.errorbar(RList, energy_fit, yerr=energy_fit_error, color='black', ls=None, linewidth=1,  marker='o', markersize=8, capsize=2, capthick=1, ecolor='black', label="PIGS", zorder=0)

	if (numbmolecules == 2):
		colorList = ['r', 'm', 'g', 'r', 'b']
		markerList = ['s', 'p', '>', 's', '<', '8', 'p']
		lsList = ['-', '--', '-.', ':']
		iplot=0
		rmlist = [5.1+0.2*i for i in range(25)]
		for jrot in range(4, 5):
			gtheta = int(2*jrot+1)
			gphi = int(2*(2*jrot+1))

			FileToBePlotEnergy_exact="/Users/tsahoo/academic-project/outputs/exact-computation/ground-state-energy-vs-Rpt-arpack-2-p-H2O-jmax"+str(jrot)+"-grid-"+str(gtheta)+"-"+str(gphi)+"-qTIP4P.txt"
			print(FileToBePlotEnergy_exact)
			exit()
			data_input = genfromtxt(FileToBePlotEnergy_exact, unpack=True, skip_header=0, skip_footer=0)
			for rm_val in rmlist:
				rm_val = "{:3.1f}".format(rm_val)
				rm_index=np.where(data_input[0]==float(rm_val))[0]
				data_input=np.delete(data_input, rm_index, axis=1)
			Rpt_exact = data_input[0]
			TotalEnergy_exact = data_input[1]/kcalmolinvKinv
			plt.plot(Rpt_exact, TotalEnergy_exact, color=colorList[iplot], ls=lsList[0], linewidth=1, marker=markerList[iplot], markersize=6, label="ARPACK Lanczos: "+r'$J=$'+' '+str(jrot))
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

def GetFigureOrderParam_vs_R(TypeCal,molecule_rot,TransMove,RotMove,variableName,RList,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,preskip,postskip,extra_file_name,src_dir,TypePlot,purpose):
	ENT_TYPE=""
	particleA=1

	indices = {'eiejx':3, 'eiejy':4, 'eiejz':4, 'eiej':6, 'eix':7, 'eiy':8, 'eiz':3}
	err_indices = {'eiejx':10, 'eiejy':11, 'eiejz':6, 'eiej':13, 'eix':14, 'eiy':15, 'eiz':5}
	title_plot = {'eiejx':'-eiejx.pdf','eiejy':'-eiejy.pdf','eiejz':'-eiejz.pdf','eiej':'-eiej.pdf','eix':'-eix.pdf','eiy':'-eiy.pdf','eiz':'-eiz.pdf'}
	ylabel_plot = {'eiejx':r'$e_i^x \cdot e_j^x$','eiejy':r'$e_i^y \cdot e_j^y$','eiejz':r'$\langle e_i^z \cdot e_j^z \rangle$','eiej':r'$e_i \cdot e_j$','eix':r'$e_i^x$','eiy':r'$e_i^y$','eiz':r'$\langle e_i^z \rangle$'}
	if (numbmolecules == 11):
		norm = {'eiejx':6, 'eiejy':6, 'eiejz':6, 'eiej':6, 'eix':7, 'eiy':7, 'eiz':7}
	if (numbmolecules == 2):
		norm = {'eiejx':1, 'eiejy':1, 'eiejz':1, 'eiej':1, 'eix':1, 'eiy':1, 'eiz':1}

	colorList = ['r', 'm', 'g', 'r', 'b']
	markerList = ['s', 'p', '>', 's', '<', '8', 'p']
	lsList = ['-', '--', '-.', ':']

	var_plotList = ['eiejz','eiz']
	iplot=0
	for var_plot in var_plotList:
		fx = np.zeros(len(RList))
		fx_err = np.zeros(len(RList))
		i=0
		for Rpt1 in RList:
			Rpt = "{:3.1f}".format(Rpt1)
			parameter1=parameter
			if ((Rpt1 >= 7.0) and (numbmolecules == 11)):
				parameter1 = 0.1
			if ((Rpt1 >= 8.0) and (numbmolecules == 2)):
				parameter1 = 0.1
			FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,float(Rpt),gfact,dipolemoment,parameterName,parameter1,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip,postskip,extra_file_name,src_dir,particleA,11)

			FileToBePlotEnergy = FilePlotName.SaveCorr+".txt"
			beads_skip_header_final, beads_skip_footer_final = GetKeyIndices(numbmolecules,molecule_rot,Rpt,FileToBePlotEnergy)	

			beads_vec1, valTau1, valCorr1, errorCorr1 = np.genfromtxt(FileToBePlotEnergy, unpack=True, usecols=[1, 2, indices[var_plot], err_indices[var_plot]], skip_header=0, skip_footer=0)
			trunc=beads_skip_header_final
			trunce=int(len(beads_vec1))-beads_skip_footer_final
			print(beads_vec1[trunc:trunce])
			valTau = valTau1[trunc:trunce]
			valCorr = valCorr1[trunc:trunce]/norm[var_plot]
			errorCorr = errorCorr1[trunc:trunce]/norm[var_plot]

			# Fitting of PIGS data
			fitting_term = "quatric" 
			tvar, corrFit, errorFit = GetFitEnergy(valTau, valCorr, errorCorr, variableName, fitting_term)
			fx[i]=corrFit[0]
			fx_err[i]=errorFit
			print(Rpt,'  ',fx[i],'  ',fx_err[i])
			i+=1

		#Plot begins
		#FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,1.0,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip,postskip,extra_file_name,src_dir,particleA,11)
		#FilePlotEnergy = FilePlotName.SaveCorr_vs_R
		#print(FilePlotEnergy+title_plot[var_plot])
		#call(["rm", FilePlotEnergy+title_plot[var_plot]])

		# Plotting of fitting data
		plt.errorbar(RList, fx, yerr=fx_err, color=colorList[iplot], ls=lsList[iplot], linewidth=1,  marker=markerList[iplot], markersize=8, capsize=2, capthick=1, ecolor=colorList[iplot], label=ylabel_plot[var_plot])
		iplot = iplot+1

	#plt.ylim(0.0,1.0)
	ymin, ymax = plt.ylim()
	midpointy = 0.5*(ymax-ymin)
	deltay = midpointy*0.15
	xmin, xmax = plt.xlim()
	midpointx = 0.5*(xmax-xmin)
	deltax = midpointx*0.15
	textpositionx = xmin+midpointx-0.25*midpointx
	textpositiony = ymin+midpointy
	plt.text(xmin+(xmax-xmin)*0.40,ymax-(ymax-ymin)*0.52,r'$N = \ $'+str(numbmolecules))
	plt.yticks(rotation=0)
	plt.xticks(rotation=0)
	#plt.xlim(2.0,10.1)

	plt.ylabel("Order parameter",labelpad=5)
	plt.xlabel(r'$r \ (\mathrm{\AA})$', labelpad=5)
	if (numbmolecules == 2):
		plt.subplots_adjust(top=0.99,bottom=0.12,left=0.11,right=0.99,hspace=0.06,wspace=0.0)
	if (numbmolecules == 11):
		plt.subplots_adjust(top=0.99,bottom=0.12,left=0.11,right=0.99,hspace=0.0,wspace=0.0)

	plt.legend(numpoints=1,loc='upper right')

	#For ticks manipulating
	plt.minorticks_on()
	plt.tick_params(axis="both", direction="in", which="minor", right=False, top=False, length=2)
	plt.tick_params(axis="both", direction="in", which="major", right=True, top=True, length=6)

	#Plot begins
	FilePlotName=support.GetFileNamePlot(TypeCal,molecule_rot,TransMove,RotMove,variableName,1.0,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules,molecule,ENT_TYPE,preskip,postskip,extra_file_name,src_dir,particleA,11)
	FilePlotEnergy = FilePlotName.SaveCorrFitvsR
	print(FilePlotEnergy+title_plot[var_plot])
	call(["rm", FilePlotEnergy+title_plot[var_plot]])

	plt.savefig(FilePlotEnergy+title_plot[var_plot], format='pdf')
	plt.show()
