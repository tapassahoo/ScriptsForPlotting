import decimal
import os
import time
from os import system
from subprocess import call
import numpy as np

import generator_pigs_plot as generator

purpose = "article"
# purpose="ppt"

method = "PIGS"
translational_move = False
rotational_move = True

molecular_system = "HF"
rotor = "HF"
numb_molecule = 2

parameter_name = "beta"
parameter_value = 0.1
dipole_moment = 1.827

numb_block = 20000
numb_pass = 200
preskip_list = [0, 10000]
postskip = 0
extra_file_name = ""

if (parameter_name == "tau"):
	variable_name = "beta"
if (parameter_name == "beta"):
	variable_name = "tau"

distance_flag=False
plot_type = "energy"

if (method == "PIGS"):
	file_name_modifier = "pigs"
elif (method == "PIMC"):
	file_name_modifier = "pimc"

final_result_path = os.path.expanduser(
	"~") + "/academic-project/outputs/final-" + file_name_modifier + "-outputs-for-plotting/"

rlist = np.arange(2.5, 3.01, 0.1, dtype=float)

if (((rotational_move) and (translational_move == False) and (plot_type ==
															  "energy")) and ((parameter_name == "tau") or (parameter_name == "beta"))):
	for preskip in preskip_list:
		for value in rlist:
			rpt_value = "{:3.1f}".format(value)
			if (parameter_name == "tau"):
				generator.get_plot_rotational_energy_vs_beta(
					final_result_path,
					method,
					molecular_system,
					rotor,
					numb_molecule,
					parameter_name,
					parameter_value,
					float(rpt_value),
					dipole_moment,
					numb_block,
					numb_pass,
					preskip,
					postskip,
					extra_file_name,
					plot_type,
					purpose)
			if ((parameter_name == "beta") and (distance_flag == False)):
				generator.get_plot_rotational_energy_vs_tau(
					final_result_path,
					method,
					molecular_system,
					rotor,
					numb_molecule,
					parameter_name,
					parameter_value,
					float(rpt_value),
					dipole_moment,
					numb_block,
					numb_pass,
					preskip,
					postskip,
					extra_file_name,
					plot_type,
					purpose)
		'''
		if ((parameter_name == "tau") and (distance_flag == True)):
			parameterName = "beta"
			parameter = beta
			generator.GetFigureEnergyRot_vs_R(method, molecule_rot, translational_move, rotational_move, parameter_name, RList, gfact, dipole_moment, parameterName, parameter, numb_block, numb_pass, numb_molecule, molecule, preskip, postskip, extra_file_name, final_results_path, plot_type, purpose)
		'''
