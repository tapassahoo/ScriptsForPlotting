import decimal
import os
import time
from os import system
from subprocess import call
import numpy as np

#import generator_pigs_plots as generator

import mypkg
import mypkg.moribs_runner as runner
#import mypkg.moribs_runner.support as support

module_path = mypkg.moribs_runner.__file__
module_path = module_path.replace('__init__.py', '')

purpose = "article"
# purpose="ppt"

translational_move = False
rotational_move = True

method = "PIGS"
parameter_name = "tau"
molecular_system = "HF"
rotor = "HF"

second_parameter_name="distance"
plot_type = "Energy"

if (method == "PIGS"):
	file_name_modifier = "pigs"
elif (method == "PIMC"):
	file_name_modifier = "pimc"

final_results_path = os.path.expanduser("~") + "/academic-project/outputs/final-"+file_name_modifier+"-outputs-for-plotting/"

if (((rotational_move == True) and (translational_move == False) and (plot_type == "Energy")) and ((parameter_name == "tau") or (parameter_name == "beta"))):
	beta = 0.1
	numb_molecule = 2
	numb_block = 20000
	numb_pass = 200
	preskip_list= [0]
	postskip = 0
	extra_file_name = "qTIP4P-"

	# compulsory parameters
	dipole_moment = -1.0

	for preskip in preskip_list:
		for Rpt in RList:
			Rpt1 = "{:3.1f}".format(Rpt)
			if (parameter_name == "beta"):
				parameterName = "tau"
				parameter = tau
				generator.get_plot_energy_vs_beta(method, molecule_rot, translational_move, rotational_move, parameter_name, float(Rpt1), gfact, dipole_moment, parameterName, parameter, numb_block, numb_pass, numb_molecule, molecule, preskip, postskip, extra_file_name, final_results_path, plot_type, purpose)
			if ((parameter_name == "tau") and (second_parameter_name != "distance")):
				parameterName = "beta"
				parameter = beta
				generator.get_plot_energy_vs_tau(method, molecule_rot, translational_move, rotational_move, parameter_name, float(Rpt), gfact, dipole_moment, parameterName, parameter, numb_block, numb_pass, numb_molecule, molecule, preskip, postskip, extra_file_name, final_results_path, plot_type, purpose)
		if ((parameter_name == "tau") and (second_parameter_name == "distance")):
			parameterName = "beta"
			parameter = beta
			generator.GetFigureEnergyRot_vs_R(method, molecule_rot, translational_move, rotational_move, parameter_name, RList, gfact, dipole_moment, parameterName, parameter, numb_block, numb_pass, numb_molecule, molecule, preskip, postskip, extra_file_name, final_results_path, plot_type, purpose)
'''
