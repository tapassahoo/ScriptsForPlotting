import decimal
import os
import time
from os import system
from subprocess import call

import numpy as np

import FigureGenerator
import mypkg.pkgMoribs.support_without_parallel as support

TransMove = False
RotMove = True

TypeCal = 'PIGS'
variableName = "tau"

TypePlot = "Energy-vs-R"
#TypePlot = "OrderParam-vs-R"
#TypePlot = "Energy"
# TypePlot="ChemPot"
# TypePlot="CorrFunc"
# TypePlot="S2"
# TypePlot="GFACTOR"


# TypePlot="COMBINE" #Comparision of "Broken Path ensemble" and "Swap+Unswap grand ensemble" methodologies
# TypePlot="COMPARE" # Importance of "ratio trick" algorithm to "Swap+Unswap grand ensemble"

molecule = "-p-H2O"
molecule_rot = "H2O"

Rpt = 8.0
purpose = "article"
# purpose="ppt"
final_results_path = os.path.expanduser("~") + "/ResultsOf" + TypeCal + "/"

if ((TypeCal == "ENT") and (TypePlot == "GFACTOR")):
    beadsRef = 21
    numbblocks = 20000
    numbpass = 100
    preskip = 10000
    postskip = 0
    extra_file_name = ""
    ENT_TYPE = "SWAPTOUNSWAP"
    parameterName = "beta"
    beta = 0.2
    parameter = beta
    plotNumber = 2
    FigureGenerator.GetFigureEntropyRT_vs_gFactor(
        TypeCal,
        molecule_rot,
        TransMove,
        RotMove,
        variableName,
        Rpt,
        parameterName,
        parameter,
        numbblocks,
        numbpass,
        molecule,
        ENT_TYPE,
        preskip,
        postskip,
        extra_file_name,
        final_results_path,
        TypePlot,
        beadsRef,
        plotNumber,
        purpose,
    )

if ((TypeCal == "ENT") and (TypePlot == "S2")):
	numbmolecules = 32
	numbblocks = 20000
	numbpass = 100
	preskip = 10000
	postskip = 0
	extra_file_name = "Ratio-Trick-"
	# extra_file_name	  = ""
	ENT_TYPE = "SWAPTOUNSWAP"

	if (variableName == "beta"):
		parameterName = "tau"
		# tau = 0.005
		tau = 0.02
		parameter = tau
		print("TAPAS")
		FigureGenerator.GetFigureEntropyRT_vs_beta(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,parameterName,parameter,numbblocks,numbpass,molecule,ENT_TYPE,preskip,postskip,extra_file_name,final_results_path,TypePlot,numbmolecules)

	if (variableName == "tau"):
		parameterName = "beta"
		beta = 0.2
		parameter = beta
		gFactor = 0.6

		FigureGenerator.GetFigureEntropyRT_vs_tau(TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,gFactor,parameterName,parameter,numbblocks,numbpass,molecule,ENT_TYPE,preskip,postskip,extra_file_name,final_results_path,TypePlot,numbmolecules)


if TypeCal == "ENT" and TypePlot == "COMBINE":
    beadsRef = 61
    numbblocks = 20000
    numbpass = 100
    preskip = 10000
    postskip = 0

    extra_file_name = ""
    ENT_TYPE = ["BROKENPATH", "SWAPTOUNSWAP"]
    parameterName = "beta"
    beta = 0.2
    parameter = beta
    numbmolecules = 2
    FigureGenerator.GetFigureEntropyRT_vs_gFactor_COMBO(
        TypeCal,
        molecule_rot,
        TransMove,
        RotMove,
        variableName,
        Rpt,
        parameterName,
        parameter,
        numbblocks,
        numbpass,
        molecule,
        ENT_TYPE,
        preskip,
        postskip,
        extra_file_name,
        final_results_path,
        TypePlot,
        beadsRef,
        numbmolecules,
        purpose,
    )
if TypeCal == "ENT" and TypePlot == "COMPARE":

    beadsRef = 61
    numbblocks = 20000
    numbpass = 100
    preskip = 10000
    postskip = 0

    extra_file_name = ""
    algorithm = ["Ratio-Trick-", "NoRatio-Trick-"]
    parameterName = "beta"
    beta = 0.2
    parameter = beta
    numbmolecules = 16
    if numbmolecules < 4:
        print("No. of rotors must be = " + str(numbmolecules))
        print()
        print("Minimum value of No. of rotors must be >= 4 ")
        exit()
    ENT_TYPE = "SWAPTOUNSWAP"
    FigureGenerator.GetFigureEntropyRT_vs_gFactor_COMPARE(TypeCal, ENT_TYPE, molecule_rot, TransMove, RotMove, variableName, Rpt, parameterName, parameter,
                                                          numbblocks, numbpass, molecule, algorithm, preskip, postskip, extra_file_name, final_results_path, TypePlot, beadsRef, numbmolecules, purpose)
"""
if (TypeCal == "ENT" and TypePlot == "GFACTOR" or TypePlot == "S2"):
	beadsRef = 21
	numbblocks			= 10000
	numbmolecules		= 2
	numbpass			= 200
	particleA			= int(numbmolecules/2)
	preskip				= 6000
	postskip			= 0

	FigureGenerator.FigureENT(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, '1.0', parameterName, parameter, numbblocks,
							  numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, particleA,TypePlot, beadsRef)
"""

"""
if (TypeCal == "ENT" and TypePlot == "COMBINE"):
	beadsRef = 101
	numbblocks			= 50000
	numbmolecules		= 2
	numbpass			= 200
	particleA			= int(numbmolecules/2)
	preskip				= 1000
	postskip			= 0
	FigureGenerator.FigureENTCOMBINE(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks,
									 numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, particleA,TypePlot, beadsRef)

if (TypePlot == "CorrFunc"):
	beadsRef = 101
	RefPoint = [3]
	FigureGenerator.FigureCorrelation(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks,
									  numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, particleA, beadsRef, RefPoint)
"""

if (TypePlot == "Energy"):
    parameterName = "beta"
    beta = 0.1
    parameter = beta
    numbmolecules = 2
    numbblocks = 10000
    numbpass = 100
    preskip = 5000
    postskip = 0
    #extra_file_name = "TIP4P-2005-"
    extra_file_name = ""

    # compulsory parameters
    dipolemoment = -1.0
    gfact = -1.0

    RList = [2.4+0.1*i for i in range(57)]
    for Rpt in RList:
        if (numbmolecules == 11):
            if (parameterName == "beta"):
                if (Rpt > 6.0):
                    beta = 0.1
                elif ((Rpt > 2.0) and (Rpt <= 6.0)):
                    beta = 0.05
            parameter = beta
        Rpt = "{:3.1f}".format(Rpt)
        FigureGenerator.FigureEnergyPIGS(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, gfact, dipolemoment, parameterName, parameter,
                                         numbblocks, numbpass, numbmolecules, molecule, "", preskip, postskip, extra_file_name, final_results_path, 1, TypePlot, 11, purpose)

if (TypePlot == "Energy-vs-R"):
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

    RList = [2.4+0.1*i for i in range(57)]

    FigureGenerator.GetFigureEnergy_vs_R(TypeCal, molecule_rot, TransMove, RotMove, variableName, RList, gfact, dipolemoment, parameterName,
                                         parameter, numbblocks, numbpass, numbmolecules, molecule, preskip, postskip, extra_file_name, final_results_path, TypePlot, purpose)

if (TypePlot == "OrderParam-vs-R"):
    parameterName = "beta"
    beta = 0.1
    parameter = beta
    numbmolecules = 11
    numbblocks = 10000
    numbpass = 100
    preskip = 0
    postskip = 0
    #extra_file_name = "TIP4P-2005-"
    extra_file_name = ""

    # compulsory parameters
    dipolemoment = -1.0
    gfact = -1.0

    RList = [2.4+0.1*i for i in range(67)]

    FigureGenerator.GetFigureOrderParam_vs_R(TypeCal, molecule_rot, TransMove, RotMove, variableName, RList, gfact, dipolemoment, parameterName,
                                             parameter, numbblocks, numbpass, numbmolecules, molecule, preskip, postskip, extra_file_name, final_results_path, TypePlot, purpose)

# End plotting ---energy
"""
if (TypePlot == "ChemPot"):
	beadsRef = 61
	FigureGenerator.FigureChemicalPotentialPIGS(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter,
												numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, particleA, beadsRef)
# End plotting ---Chemical Potential
if (TypeCal == "PIGS" and TypePlot == "GFACTOR"):
	beadsRef   = 21
	numbblocks = 10000
	numbmolecules  = 2
	numbpass   = 50
	preskip    = 5000
	postskip   = 0
	particleA			= int(numbmolecules/2)
	FigureGenerator.FigureAngleDistributionGfactor(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter,
												   numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, particleA,TypePlot, beadsRef)
"""

"""
if (TypeCal == "PIGS" and TypePlot == "GFACTOR"):
	beadsRef   = 41
	numbblocks = 10000
	numbpass   = 50
	FigureGenerator.FigureAngleDistribution(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, particleA,TypePlot, beadsRef)
"""
