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
#TypePlot = "Entropy"
#TypePlot = "OrderParam"
#TypePlot = "Energy"
# TypePlot="ChemPot"
# TypePlot="CorrFunc"

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
	numbmolecules = 11
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
	#print(RList)
	
	for preskip in preskipList:
		for Rpt in RList:
			Rpt1 = "{:3.1f}".format(Rpt)
			if (variableName == "beta"):
				parameterName = "tau"
				parameter = tau
				generator.GetFigureEnergyRot_vs_beta(TypeCal, molecule_rot, TransMove, RotMove, variableName, float(Rpt1), gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, preskip, postskip, extra_file_name, final_results_path, TypePlot, purpose)
			if ((variableName == "tau") and (variableName1 != "distance")):
				parameterName = "beta"
				parameter = beta
				generator.GetFigureEnergyRot_vs_tau(TypeCal, molecule_rot, TransMove, RotMove, variableName, float(Rpt), gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, preskip, postskip, extra_file_name, final_results_path, TypePlot, purpose)
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


''''
if ((RotMove == True) and (TransMove == False) and (TypePlot == "Energy")):
	parameterName = "beta"
	beta = 0.1
	parameter = beta
	numbmolecules = 2
	numbblocks = 10000
	numbpass = 100
	preskip = 0
	postskip = 0
	#extra_file_name = "TIP4P-2005-"
	extra_file_name = ""

	# compulsory parameters
	dipolemoment = -1.0
	gfact = -1.0

	rmin = 2.4
	rmax =10.0
	dr = 0.1
	nr = int(((rmax-rmin)+dr*0.5)/dr)
	nr = nr+1
	print(nr)
	RList = [rmin+dr*i for i in range(nr)]
	print(RList)

	generator.GetFigureEnergy_vs_R(TypeCal, molecule_rot, TransRotMove, RotMove, variableName, RList, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, preskip, postskip, extra_file_name, final_results_path, TypePlot, purpose)

if ((TransMove == True)  and (RotMove == True) and (TypePlot == "Energy")):
	beta = 0.05
	tau = 0.001
	numbmolecules = 2
	numbblocks = 10000
	numbpass = 500
	preskipList= [0, 5000, 7000, 8000]
	postskip = 0
	#extra_file_name = "COM-Bisection-Norm-moves-"
	#extra_file_name = "COM-and-Bisection-moves-"
	#extra_file_name = ""
	#extra_file_name = "qTIP4PF-"
	extra_file_name = "qspcfw-"

	# compulsory parameters
	dipolemoment = -1.0
	gfact = -1.0

	RList = [6.0]

	for preskip in preskipList:
		if (variableName == "beta"):
			parameterName = "tau"
			parameter = tau
			generator.GetFigureEnergyTransRot_vs_beta(TypeCal, molecule_rot, TransMove, RotMove, variableName, RList, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, preskip, postskip, extra_file_name, final_results_path, TypePlot, purpose)
		if (variableName == "tau"):
			parameterName = "beta"
			parameter = beta
			generator.GetFigureEnergyTransRot_vs_tau(TypeCal, molecule_rot, TransMove, RotMove, variableName, RList, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, preskip, postskip, extra_file_name, final_results_path, TypePlot, purpose)
'''

if (TypePlot == "ExactEnergy"):
	generator.GetFigureEnergy_vs_R_lanczos()
