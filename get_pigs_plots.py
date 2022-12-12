import os, getpass, subprocess, argparse
from datetime import datetime
from termcolor import colored
import numpy as np
import generator_pigs_plot as generator

purpose = "article"
# Valuable informations about the simulations
method = "PIGS"
translational_move = False
rotational_move = True
#
molecular_system = "HF"
rotor = "HF"
numb_molecule = 2
#
parameter_name = "beta"
parameter_value = 0.1
dipole_moment = 1.827
#
numb_block = 20000
numb_pass = 200
preskip_list = [0, 10000, 15000]
postskip = 0
extra_file_name = ""

if (parameter_name == "tau"):
	variable_name = "beta"
if (parameter_name == "beta"):
	variable_name = "tau"

distance_flag = True
#plot_type = "energy"
plot_type = "order_parameter"

if (method == "PIGS"):
	file_name_modifier = "pigs"
elif (method == "PIMC"):
	file_name_modifier = "pimc"

rlist = np.arange(3.0, 10.01, 0.2, dtype=float)
energy_per_neighbours = True

final_result_path = os.path.join(os.path.expanduser("~"), "academic-project", "output", "final-" + file_name_modifier + "-outputs-for-plotting")

print("*"*80 + "\n")
print(colored("Developer ::".ljust(30),"blue") + colored("Dr. Tapas Sahoo", "yellow") + "\n")
now = datetime.now() # current date and time
date_time = now.strftime("%d/%m/%Y, %H:%M:%S")
print("date and time ::".capitalize().ljust(29), date_time, "\n")

if (((rotational_move == True) and (translational_move == False)) and ((parameter_name == "tau") or (parameter_name == "beta"))):
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
		if ((parameter_name == "beta") and (distance_flag == True)):
			if (plot_type == "energy"):
				generator.get_plot_rotational_energy_vs_rpt(
						final_result_path,
						method,
						molecular_system,
						rotor,
						numb_molecule,
						parameter_name,
						parameter_value,
						rlist,
						dipole_moment,
						numb_block,
						numb_pass,
						preskip,
						postskip,
						extra_file_name,
						plot_type,
						energy_per_neighbours,
						purpose)
			if (plot_type == "order_parameter"):
				generator.get_plot_order_parameter_vs_rpt(
						final_result_path,
						method,
						molecular_system,
						rotor,
						numb_molecule,
						parameter_name,
						parameter_value,
						rlist,
						dipole_moment,
						numb_block,
						numb_pass,
						preskip,
						postskip,
						extra_file_name,
						plot_type,
						energy_per_neighbours,
						purpose)
