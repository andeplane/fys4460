# mdconfig.py
# Module to simplify MD simulations. 
# If logging is enabled, every terminal command is written to run_log.txt
# 
# Anders Hafreager, 2013

import subprocess
import os
import logging
from datetime import datetime

class MD:
	def __init__(self, compiler = "icpc", dt=0.02, logging_enabled=True):
		"""
		Initializes an object with parameters.
		"""
		current_directory = os.getcwd()
		# if os.path.basename(current_directory) == "base_code":
		# 	# We don't want to modify the base_code folder
		# 	print "Cannot run code from base_code folder. Please copy the program into another folder. Aborting!"
		# 	exit()

		self.dt = dt
		self.FCC_b = 5.72
		self.load_state = True
		self.thermostat_enabled = False
		self.timesteps = 5000
		self.create_movie_files = False
		self.movie_every_n_frame = 1
		self.statistics_interval = 100
		self.nodes_x = 1
		self.nodes_y = 1
		self.nodes_z = 1
		self.unit_cells_x = 5
		self.unit_cells_y = 5
		self.unit_cells_z = 5
		self.many_frozen_atoms = False
		self.gravity_force = 0.002
		self.gravity_direction = 0
		self.temperature = 100
		self.thermostat_relaxation_time = 1
		self.r_cut = 2.5
		self.mass = 39.948

		self.compiler = compiler
		self.constants = dict()
		self.constants['boltzmann'] = 1.3806488e-23

		self.test_mode = False
		self.logging_enabled = logging_enabled
		self.total_timesteps = 0

		logging.basicConfig(filename='run_log.txt',level=logging.INFO)
		self.log("System initialized")
		
	def log(self, message):
		if self.logging_enabled:
			logging.info("#LOG "+str(datetime.now())+" : "+message)

	def load_state(self, path):
		if self.test_mode: return

		self.clean()
		self.log("Loading state from "+path)
		self.run_command("cp -r "+path+"/ ./")

	def save_state(self, path):
		if self.test_mode: return
		self.log("Saving state to "+str(path))
		self.run_command("mkdir -p "+path)
		self.run_command("cp -r log state_files statistics run_log.txt md.ini Tocontinue "+path)
		

	def run_command(self, cmd):
		"""
		Runs command <cmd> in terminal.
		
		:param cmd: Command to run.
		:type cmd: str.
		
		"""
		if self.logging_enabled:
			logging.info("#CMD "+str( datetime.now() )+" : "+cmd)
		subprocess.call(cmd, shell=True)

	def clean(self):
		if self.test_mode: return

		self.run_command('rm log')
		self.run_command('rm -rf state_files/*')
		self.run_command('rm MD.ini')
		self.run_command('rm Tocontinue')
		
	def reset(self):
		"""
		Deletes files created by the md program.
		
		"""
		if self.test_mode: return

		self.log('Resetting root folder')
		self.clean()
		self.run_command('mkdir state_files')
		self.run_command('mkdir movie_files')
		self.run_command('mkdir statistics')
		self.run_command('echo 0.0 0 > Tocontinue')

	def compile(self, path = "./program/md", name="md", skip_compile = False):
		if not skip_compile:
			current_directory = os.getcwd()
			os.chdir(path)
			self.run_command('qmake')
			#self.run_command('make clean')
			self.run_command('make')
			os.chdir(current_directory)
			move_command = 'mv '+path+'/md ./%s' % name
			self.run_command(move_command)

		if not os.path.isfile("./"+name):
			print "Executable ./"+name+" is not compiled, aborting!"
			exit()
		return './%s' % name

	def create_config_file(self, config_file='md.ini.original'):
		"""
		Creates config file with current settings.
		
		"""
		
		self.log("Creating config file")
		original_file = open("program/000_config_files/"+config_file,'r');
		output_file   = open('md.ini','w');

		for line in original_file:
			line = line.replace('__FCC_b__',str(self.FCC_b) )
			line = line.replace('__load_state__',str(self.load_state).lower() )
			line = line.replace('__thermostat_enabled__',str(self.thermostat_enabled).lower() )
			line = line.replace('__create_movie__',str(self.create_movie_files).lower() )
			line = line.replace('__timesteps__',str(self.timesteps) )
			line = line.replace('__temperature__',str(self.temperature) )
			line = line.replace('__movie_every_n_frame__',str(self.movie_every_n_frame) )
			line = line.replace('__dt__',str(self.dt) )
			line = line.replace('__statistics_interval__',str(self.statistics_interval) )
			line = line.replace('__nodes_x__',str(self.nodes_x) )
			line = line.replace('__nodes_y__',str(self.nodes_y) )
			line = line.replace('__nodes_z__',str(self.nodes_z) )
			line = line.replace('__unit_cells_x__',str(self.unit_cells_x) )
			line = line.replace('__unit_cells_y__',str(self.unit_cells_y) )
			line = line.replace('__unit_cells_z__',str(self.unit_cells_z) )
			line = line.replace('__many_frozen_atoms__',str(self.many_frozen_atoms).lower() )
			line = line.replace('__gravity_force__',str(self.gravity_force) )
			line = line.replace('__gravity_direction__', str(self.gravity_direction) )
			line = line.replace('__thermostat_relaxation_time__', str(self.thermostat_relaxation_time) )
			line = line.replace('__r_cut__', str(self.r_cut) )
			line = line.replace('__mass__', str(self.mass) )
			
			output_file.writelines(line)

		original_file.close()
		output_file.close()
	
	def run(self, executable="md"):
		"""
		Runs specified executable, puts data into a folder, named <project_name>-<name>.
		
		:param executable: Specifies which executable to run.
		:type executable: str.
		:param name: Specifies which output data to use. Default is the last data output.
		:type name: str.
		
		"""

		self.total_timesteps += self.timesteps
		if self.test_mode: return
		
		if not os.path.isfile(executable):
			print "Executable "+executable+" is not compiled, aborting!"
			exit()
		
		self.log("Running executable "+executable)
		now = datetime.now()
		num_nodes = self.nodes_x*self.nodes_y*self.nodes_z
		self.run_command("mpirun -n %d ./%s | tee log" % (num_nodes, executable))
		t1 = (datetime.now() - now).seconds
		steps_per_second = self.timesteps / max(t1,1)

		self.log("Process used %d seconds (%d timesteps per second)" % ( t1, steps_per_second ))

	def create_cylinder(self, radius, cylinder_dimension, center_x, center_y, center_z):
		num_nodes = self.nodes_x*self.nodes_y*self.nodes_z
		self.run_command("%s -O3 program/create_cylinder/create_cylinder.cpp -o create_cylinder" % self.compiler)
		self.run_command("./create_cylinder %d %f %d %f %f %f" % (num_nodes, radius, cylinder_dimension, center_x, center_y, center_z))

	def create_spheres(self, num_spheres, r_min, r_max, x_max, y_max, z_max):
		num_nodes = self.nodes_x*self.nodes_y*self.nodes_z
		self.run_command("%s -O3 program/create_spheres/create_spheres.cpp -o create_spheres" % self.compiler)
		self.run_command("./create_spheres %d %d %f %f %f %f %f" % (num_nodes, num_spheres, r_min, r_max, x_max, y_max, z_max))

	def create_spheres(self, num_spheres, r_min, r_max, x_max, y_max, z_max):
		num_nodes = self.nodes_x*self.nodes_y*self.nodes_z
		self.run_command("%s -O3 program/create_spheres/create_spheres.cpp -o create_spheres" % self.compiler)
		self.run_command("./create_spheres %d %d %f %f %f %f %f" % (num_nodes, num_spheres, r_min, r_max, x_max, y_max, z_max))
	def create_movie(self, frames):
		num_nodes = self.nodes_x*self.nodes_y*self.nodes_z
		self.run_command("%s -O3 program/create_movie/create_movie.cpp -o create_movie" % self.compiler)
		self.run_command("./create_movie %d %d" % (num_nodes, frames) )

	def prepare_new_system(self):
		self.timesteps = 1
		self.load_state = False
		self.create_config_file()
		self.load_state = True