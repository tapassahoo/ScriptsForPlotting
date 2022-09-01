import decimal
import os
import time
from os import system
from subprocess import call

import numpy as np

import FigureGeneratorForPIGS as generator
import mypkg.pkgMoribs.support_without_parallel as support

TransMove = False
RotMove = True

TypeCal = 'PIGS'
variableName = "tau"
#variableName = "beta"
variableName1="distance"

TypePlot = "Energy"
#TypePlot = "ExactEnergy"
#TypePlot = "OrderParam"
#TypePlot = "Energy"
# TypePlot="ChemPot"

molecule = "-p-H2O"
molecule_rot = "H2O"

purpose = "article"
# purpose="ppt"
if (TypeCal == 'PIGS'):
	cal_method = 'pigs'
final_results_path = os.path.expanduser("~") + "/academic-project/outputs/results-of-" + cal_method + "/"

if (((RotMove == True) and (TransMove == False) and (TypePlot == "Energy")) and ((variableName == "tau") or (variableName == "beta"))):
	beta = 0.1
	tau = 0.002
	numbmolecules = 2
	numbblocks = 20000
	numbpass = 200
	preskipList= [0]
	postskip = 0
	extra_file_name = "qTIP4P-"

	# compulsory parameters
	dipolemoment = -1.0
	gfact = -1.0

	rmin=2.5
	rmax=2.7
	dr=0.02
	nr = int(((rmax-rmin)+dr*0.5)/dr)
	nr += 1
	RList = [rmin+dr*i for i in range(nr)]
	RList += [2.75]
	rmin = 2.8
	rmax = 5.0
	dr = 0.1
	nr = int(((rmax-rmin)+dr*0.5)/dr)
	nr += 1
	RList += [rmin+dr*i for i in range(nr)]
	rmin = 5.2
	rmax = 10.0
	dr = 0.2
	nr = int(((rmax-rmin)+dr*0.5)/dr)
	nr += 1
	print(nr)
	RList += [rmin+dr*i for i in range(nr)]
	print(RList)
	
	for preskip in preskipList:
		for Rpt in RList:
			Rpt1 = "{:3.1f}".format(Rpt)
			if (variableName == "beta"):
				parameterName = "tau"
				parameter = tau
				generator.get_plot_energy_vs_beta(TypeCal, molecule_rot, TransMove, RotMove, variableName, float(Rpt1), gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, preskip, postskip, extra_file_name, final_results_path, TypePlot, purpose)
			if ((variableName == "tau") and (variableName1 != "distance")):
				parameterName = "beta"
				parameter = beta
				generator.get_plot_energy_vs_tau(TypeCal, molecule_rot, TransMove, RotMove, variableName, float(Rpt), gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, preskip, postskip, extra_file_name, final_results_path, TypePlot, purpose)
		if ((variableName == "tau") and (variableName1 == "distance")):
			parameterName = "beta"
			parameter = beta
			generator.GetFigureEnergyRot_vs_R(TypeCal, molecule_rot, TransMove, RotMove, variableName, RList, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, preskip, postskip, extra_file_name, final_results_path, TypePlot, purpose)

if (((RotMove == True) and (TransMove == False) and (TypePlot == "OrderParam")) and ((variableName == "tau") or (variableName == "beta"))):
	beta = 0.1
	tau = 0.002
	numbmolecules = 11
	numbblocks = 20000
	numbpass = 200
	preskipList= [0]
	postskip = 0
	#extra_file_name = ""
	extra_file_name = "qTIP4P-"

	# compulsory parameters
	dipolemoment = -1.0
	gfact = -1.0

	rmin = 2.7
	rmax = 5.0
	dr = 0.1
	nr = int(((rmax-rmin)+dr*0.5)/dr)
	nr += 1
	RList = [rmin+dr*i for i in range(nr)]
	rmin = 5.2
	rmax = 8.0
	dr = 0.2
	nr = int(((rmax-rmin)+dr*0.5)/dr)
	nr += 1
	print(nr)
	RList += [rmin+dr*i for i in range(nr)]
	print(RList)
	
	for preskip in preskipList:
		if ((variableName == "tau") and (variableName1 == "distance")):
			parameterName = "beta"
			parameter = beta

			generator.GetFigureOrderParam_vs_R(TypeCal, molecule_rot, TransMove, RotMove, variableName, RList, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, preskip, postskip, extra_file_name, final_results_path, TypePlot, purpose)
