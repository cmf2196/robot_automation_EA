"""
* Author: Connor Finn
* This script details a physics simulator class.
* The simulator object will be used to evaluate robots learning to move 
* Multiple robots can be simulated
* The simulator has the option of visualizing the moving robot using the vpython graphics library
* the robot can be depicted in a solid form or as a mass / spring representation
"""

import numpy as np
from vpython import *
import matplotlib.pyplot as plt
import time



## constants we will be using

T = 0					# global time variable
g = np.array((0 , -9.8 , 0))	# gravity
dt = 0.0001				# time step
damp = 0.999			# add dampening to make it more realistic
mu_s = .61				# coeficient of static friction: oak on leather
mu_d = .52				# coeficient of dynamic frriction: oak on leather
K_g = -100000			# spring constant of the ground
w = 62.8				# frequency gives 1000 time steps per robot cycle


# make the robot class - going to have all kinds of useful functions
class simulator:
	
	def __init__(self):
		self.T_stop = 1 
		self.springs_shown = [] 		
		self.spheres_shown = []
		self.triangles_shown = []
		self.robot_list = []
		self.mass_list = []
		self.spring_list = []



	def set_robot(self , robot):
		"""
		* this function will allow us to simulate more than one mass at a time
		* it takes in a list of robots, and sets the simulator's mass and spring lists 
		  to be the sum of the masses and springs of the robots

		"""
		mass_list = self.mass_list						# get the simulator's mass list
		spring_list = self.spring_list					# get the simulator's spring list
		self.robot_list += [robot]						# add robot to the list of robots
		for mass in robot.mass_list:
			mass_list += [mass]					# add the masses and springs to the simulator
		for spring in robot.spring_list:
			spring_list += [spring]



	def remove_robot(self , robot):
		"""
		* this takes in a robot, and removes it from the robot list
		* will return an error statement if the robot is not in the list
		"""


		if robot in self.robot_list:																	# check if robot is in the list
			self.mass_list = [item for item in self.mass_list if item not in robot.mass_list]			# remove masses
			self.spring_list = [item for item in self.spring_list if item not in robot.spring_list]		# remove springs
			self.robot_list.remove(robot)																# remove the robot
		else:													
			print( "this robot is not in the simulator!")												# if it isn't return an error message

# _____ Graphics methods: 

	def show_solid(self):

		""" 
	* this function takes in the mass array and will show them using vpython
		"""

		num_robots = len(self.robot_list)  					# get the number of robots in the simulator

		# make the floor: 
			# wood -> use interest_box
			# without internet -> use offline_box
		
		interest_box=box(pos=vector(0,-.005,0),size=vector(10,.005,10), texture = 'https://s3.amazonaws.com/glowscript/textures/wood_texture.jpg')
		#offline_box=box(pos=vector(0,-.005,0),size=vector(10,.005,10), color = color.white)

		# a marker to indicate where the robot started it's path
		rad = np.sqrt(2) * .05
		rod = cylinder(pos=vector(0,-.005,0),  axis=vector(0,.005,0), radius=rad , color = color.red)

		# this fills a list of the sides which are on the exterior of the robot
			# we want to shade these in to give the robot a solid appearance
			# make use of the robot's open side list that we used for developing robots from a 'seed tetrahedron'
			# triangles will be filled with vpython triangle objects 

		triangles = []
		for robot in self.robot_list:
			# start with the bottom of the seed tetrahedron (we did not allow expansion off of here)
			triangles_to_make = [[robot.mass_list[0] , robot.mass_list[1] , robot.mass_list[2]]]	
			# add all exterior sides to the list ( we only neec their list of masses)
			for side in robot.sides:
				triangles_to_make += [side.mass_list]

			# build a triangle object for set of three vertices
			for j in range(len(triangles_to_make)):
				curr_triangle = triangles_to_make[j]										

				# get the masses from our current side
				m1 = curr_triangle[0]
				m2 = curr_triangle[1]				
				m3 = curr_triangle[2]
				# get their positions ( state is a matrix = [pos , vel , accel])
				m1_pos = m1.state[0]
				m2_pos = m2.state[0]
				m3_pos = m3.state[0]
				
				v1 = vector( m1_pos[0] , m1_pos[1] , m1_pos[2])
				v2 = vector( m2_pos[0] , m2_pos[1] , m2_pos[2])
				v3 = vector( m3_pos[0] , m3_pos[1] , m3_pos[2])

			# color:
				# 1) color based on material

			#	a = vertex( pos=v1 , color = m1.color)
			#	b = vertex( pos=v2 , color = m1.color)
			#	c = vertex( pos=v3 , color= m1.color )

				# 2) color for appearance
			
				a = vertex( pos=v1 , color = color.red)
				b = vertex( pos=v2 , color = color.green)
				c = vertex( pos=v3 , color= color.blue )

				# add the triangle objects 

				T = triangle( vs = [a , b , c])
				triangles += [T]

		self.triangles_shown = triangles


	def update_solid_image(self ):
		""" 
	* This function updates the solid solid depiction of a robot in real time.
		"""

		# get the list of vpython triangle objects
		triangles = self.triangles_shown		

		# reevaluate the locations of the sides
		# this needs to be done robot by robot
		for i in range(len(self.robot_list)):
			robot = self.robot_list[i]
			# start with the bottom of the seed tetrahedron (we did not allow expansion off of here)
			triangles_to_make = [[robot.mass_list[0] , robot.mass_list[1] , robot.mass_list[2]]]
			
			for side in robot.sides:
				# add a list of masses associated with each external side to a list
				triangles_to_make += [side.mass_list]
			
			for j in range(len(triangles_to_make)):
				curr_triangle = triangles_to_make[j]										
				
				m1 = curr_triangle[0]					# get each point mass object
				m2 = curr_triangle[1]					# 		 ""    ""				
				m3 = curr_triangle[2]					# 		 ""    ""

				m1_pos = m1.state[0]					# get the positions of each mass
				m2_pos = m2.state[0]					# 		 ""    ""
				m3_pos = m3.state[0]					# 		 ""    ""
				
				# convert the position arrays into vectors
				vec1 = vector( m1_pos[0] , m1_pos[1] , m1_pos[2])
				vec2 = vector( m2_pos[0] , m2_pos[1] , m2_pos[2])
				vec3 = vector( m3_pos[0] , m3_pos[1] , m3_pos[2])

				# update the position for each vertex
				# the 32 is hardcoded and is specific for a robot comprised of 15 tetrahedron
				t = triangles[ 32 * i + j]
				t.v0.pos.x = m1_pos[0]
				t.v0.pos.y = m1_pos[1]
				t.v0.pos.z = m1_pos[2]
				t.v1.pos.x = m2_pos[0]
				t.v1.pos.y = m2_pos[1]
				t.v1.pos.z = m2_pos[2]
				t.v2.pos.x = m3_pos[0]
				t.v2.pos.y = m3_pos[1]
				t.v2.pos.z = m3_pos[2]	

	
	def show_masses(self ):
		""" 
	* This function will visualize a robot as a combination of masses and springs. ( useful for visualizing a robots representation) 
	* Point masses are visualized as vpython spheres
	* springs are visualized as vpython curves
		* different material springs are represented by different colors
		"""

		# get the information for the masses / springs in the simulator object (could be one or more robots)
		num_masses = len(self.mass_list)
		num_springs = len(self.spring_list)
		spring_list = self.spring_list
		mass_list = self.mass_list
		
		# make the floor: 
			# wood -> use interest_box
			# without internet -> use offline_box
		interest_box=box(pos=vector(0,-.005,0),size=vector(10,.005,10), texture = 'https://s3.amazonaws.com/glowscript/textures/wood_texture.jpg')
		#offline_box=box(pos=vector(0,-.005,0),size=vector(10,.005,10), color = color.white)

		# a marker to indicate where the robot started it's path
		rad = np.sqrt(2) * .05
		rod = cylinder(pos=vector(0,-.005,0),  axis=vector(0,.005,0), radius=rad , color = color.red)

		r =.01				# radius of masses
		spheres = []		# this will be filled with vpython sphere objects
		springs = []		# this will be filled with vpython curve objects

		# Point masses
		for i in range(num_masses):
			m = self.mass_list[i]				# point mass object
			curr_pos = m.state[0]				# position of the mass
			
			# build the vpython ball object and add it to our list
			ball = sphere(pos = vector(curr_pos[0] , curr_pos[1] , curr_pos[2]),radius = r,color = color.black)
			spheres += [ball]

		# Springs 
		for j in range(num_springs):
			spr = self.spring_list[j]			# get the spring object
			m1 = spr.mass_1						# get the two masses connected by the spring
			m2 = spr.mass_2						# 				""    "" 

			m1_pos = m1.state[0]				# get the position arrays for each mass
			m2_pos = m2.state[0]				# 				""    ""

			# convert position arrays to vectors
			v1 = vector( m1_pos[0] , m1_pos[1] , m1_pos[2])
			v2 = vector( m2_pos[0] , m2_pos[1] , m2_pos[2])

			# build the curve object
			c = curve(pos=[v1,v2], color=spr.color, radius=0.0025)		
			# add it to our list
			springs += [c]			

		# assign the lists of spheres and curves to our simulator object
		self.spheres_shown = spheres
		self.springs_shown = springs



	def reset(self ):
		"""
		* this function resets to a blank canvas
		* After deleting a robot, this function is called to remove the robot from the vpython visualization.
		* NOTE: if you want to set a robot back to starting location and orientation: The robot needs to be rebuilt.

		"""

		# Solid Representation 
		if self.triangles_shown != []:			
			solid = True
			for i in range(len(self.triangles_shown)):	# this loop will remove all triangles from window
				triangle = self.triangles_shown[i]
				triangle.visible = False
				del triangle	

		# Mass Spring Representation		
		else:
			for i in range(len(self.spheres_shown)):	# this loop will remove all spheres from window
				sphere = self.spheres_shown[i]
				sphere.visible = False
				del sphere						

			for i in range(len(self.springs_shown)):	# this loop will remove all springs from window
				spring = self.springs_shown[i]
				spring.visible = False
				del spring	

		# reset the simulator object's lists of spheres / triangles / curves  to empty lists
		self.triangles_shown = []	
		self.triangles_to_make = []
		self.spheres_shown = []
		self.springs_shown = []



	def update_image(self):
		""" 
	* This function updates the mass / spring representation of the robot in vpython canvas
		"""
		r =.01# radius of masses

		spheres = self.spheres_shown					# get the sphere objects
		springs = self.springs_shown					# get the curve objects

		num_masses = len(self.mass_list)				# how many masses in simulator
		num_springs = len(self.spring_list)				# how many springs in simulator

		# update the sphere objects
		for i in range(num_masses):					
			m = self.mass_list[i]											# get point mass
			curr_pos = m.state[0]											# get mass position
			ball = spheres[i]												# get a sphere object
			ball.pos = vector(curr_pos[0] , curr_pos[1] , curr_pos[2])		# update the object's position
			
		# update the curve objects
		for j in range(num_springs):
			spr = self.spring_list[j]					# get the spring object
			m1 = spr.mass_1								# get the masses connected by the spring
			m2 = spr.mass_2								# 			""   ""

			m1_pos = m1.state[0]						# get the position arrays for each mass
			m2_pos = m2.state[0]						# 			""   ""

			v1 = vector( m1_pos[0] , m1_pos[1] , m1_pos[2])			# convert pos arrays to pos vectors
			v2 = vector( m2_pos[0] , m2_pos[1] , m2_pos[2])			# 			""   ""
			
			c = springs[j]								# get a curve object	
			c.clear()									# clear its pos attributes
			c.append([v1,v2])							# update its location with the new pos vectors
	


# _____ The Simulation: 


	def simulate(self ,  T_stop , show):
		# simulation
		# show is an int (0 = no visualization  , 1 = mass spring visualization , 2 = solid visualization)

		
		global T

		# visualize the robots if indicated
		if show == 1:
			self.show_masses()
					
		if show == 2:				
			self.show_solid()

		# putting this here allows for it to be input dependant
		self.T_stop = T_stop
		self.kinetic_energy =np.zeros(int(T_stop / dt) + 2)
		self.potential_energy =np.zeros(int(T_stop / dt) + 2) 
		self.spring_energy=np.zeros(int(T_stop / dt) + 2) 
		self.total_energy =np.zeros(int(T_stop / dt) + 2) 
		self.time =np.zeros(int(T_stop / dt) + 2) 


		num_springs= len(self.spring_list)
		num_masses = len(self.mass_list)


		# If you want to record the robotL: Use this delay. 
		# IMPORTANT: Comment out the delay if you are running a learning algorithm

		#time.sleep(2)
		
		# reset time to 0
		T =0 

		force_list = [0] * num_masses												# need a list of forces for each mass
		ground_force_list = [0] * num_masses										# place holder
		ground_energy_list = [0] * num_masses										# for plotting energy of system

		tally = 0																	# stop watch of sorts
		
		start = time.time()															
		while T < T_stop:	
			self.time[tally] = time.time() - start 									# keep track of elapsed time
			temp_force_list = [np.zeros(3)] * num_masses 							# temperary place holder

#_____________ Spring Forces ____________
			for k in range(num_springs): 											# evaluate spring forces
				spring_temp = self.spring_list[k]									# get the current spring	
					
				m1 = spring_temp.mass_1												# get associated masses
				m2 = spring_temp.mass_2												# 		"  "
																					
				curr_length = m1.state[0] - m2.state[0]								# get the length of the spring 
				curr_length2 = m2.state[0] - m1.state[0] 
				dist = np.linalg.norm(curr_length)
				unit_vec1 = curr_length / dist 										# get the unit vector
				unit_vec2 = -1 * unit_vec1 											# need unit vector in both directions

				rest = spring_temp.rest_length 										# get the rest length
				Lo = rest[0] 														# 		"  "
				rest_dist = np.linalg.norm(Lo) + rest[1] * np.sin(rest[2] + T * w) 	# 		"  "

				force =  -1 * spring_temp.spring_constant * (dist - rest_dist)		# get the force scalar
 
				F_spring1 = force * unit_vec1   		# force * unit vector is force vector on each mass 
				F_spring2 = force * unit_vec2 			# 			"   "

				# add spring force to the foce on each masss
				index_1 = self.mass_list.index(m1)									# get indices for the masses
				index_2 = self.mass_list.index(m2)
				
				# add it in!
				temp_force_list[index_1] = temp_force_list[index_1] +  F_spring1 	# update the force list at the proper locations
				temp_force_list[index_2] = temp_force_list[index_2] +  F_spring2

				#add to the potential energy 
				self.spring_energy[tally] = self.spring_energy[tally] + spring_temp.spring_constant * 0.5 * (dist - rest_dist)**2				

# __________ Gravity, Friction and Ground Force __________
			static = [0] * num_masses
			for i in range(num_masses):
			# ______ Gravity ______
				mass_temp = self.mass_list[i]		         					     # get the point mass object
				temp_force_list[i] = temp_force_list[i] + g * mass_temp.mass 		 # add gravity
				
				# If it is below ground ...
				if mass_temp.state[0][1] <= 0: 

					# ______ Friction ______
					normal = abs(temp_force_list[i][1])							# get the normal force
					horizontal = np.delete(temp_force_list[i] , 1)				# get the horizontal forces
					horizontal_magnitude = np.linalg.norm(horizontal)

					# _____ static ____
					if horizontal_magnitude < normal * mu_s:	
						"""
						" No horizontal motion and no horziontal forces "
						"""
						static[i] = 1 							# need to know which masses don't move on ground
						temp_force_list[i][0] = 0
						temp_force_list[i][2] = 0

					# _____ dynamic ______ 
					else:
						horizontal_unit = horizontal / horizontal_magnitude
						friction_force = normal * mu_d										# dynamic friction force
						friction = horizontal_unit * friction_force * -1 					# in the oposite direction of the horziontal force
						friction_array = np.array((friction[0] , 0 , friction[1]))			# need to put y force of zero back in
						temp_force_list[i] = temp_force_list[i] - friction_array 			# add to force for the mass


				# ____ Ground force _____
					# model the ground as a spring with a large K
					bounce_force = mass_temp.state[0][1] * K_g			 					# calculate the spring foce
					adjust_force = np.array((0 , bounce_force , 0))	 						# only in the upward direction
					temp_force_list[i] = temp_force_list[i] + adjust_force 					# add it to our force counters
					ground_force_list[i] = bounce_force
					self.spring_energy[tally] = self.spring_energy[tally] + -.5 * K_g * mass_temp.state[0][1] ** 2  
	 				
				
# ___________ Update Kinematics __________
			for p in range(num_masses):
				force_list[p] = temp_force_list[p] 			# update the official force list

			for j in range( num_masses):
				curr_force = force_list[j]					# get the force
				curr_mass = self.mass_list[j] 				# get the mass
				ground_energy = ground_energy_list[j]		# get the ground force

				# Update energy for plots
				self.kinetic_energy[tally] = self.kinetic_energy[tally] + .5 * curr_mass.mass * np.linalg.norm(curr_mass.state[1])**2
				self.potential_energy[tally] = self.potential_energy[tally] + curr_mass.mass * -1 *  g[1] * curr_mass.state[0][1]  
				
				# update the state				
				# methodology also brought in from intro to robotics course work

				mass = curr_mass.mass 					# get mass in kg
				arr1 = np.array(((1 , dt , 0) , ( 0 , 1 , 0) , ( 0 , 0 , 0)) )
				arr2 = np.array((dt * dt / (mass) , dt / mass  , 1/ mass ) ).reshape(3 , 1)

				curr_mass.state = np.dot(arr1 , curr_mass.state ) + np.dot(arr2 , curr_force.reshape(1 , 3)) 	# apply manipulation due to forces
				curr_mass.state[1] = curr_mass.state[1] * damp 													# take into account damping ration

				if static[j] == 1:					 		# in the event we can't overcome static friction				
					curr_mass.state[1][0] = 0				# x velocity = 0
					curr_mass.state[1][2] = 0				# z velocity = 0

			# update energy for plotting
			self.total_energy[tally] = self.kinetic_energy[tally] + self.potential_energy[tally] + self.spring_energy[tally]
			tally += 1	
			T += dt
			
			# visualize if we indicated it
			if show == 1:					# this is for mass spring representation
				if T*100000 %20:
					self.update_image()
					
			if show == 2:					# this is for solid representation
				if T*100000 %20:				
					self.update_solid_image()


# _____ Plotting methods: 

	def plot_energy(self):
		"""
		* plot the energy in the system
		* if material is chosen such that no muscles are selected, and damping is 1 the energy should be constant
		* ( no breathing )
		"""

		x = np.arange(int(self.T_stop / dt) + 2)
		plt.plot( x[:20000] , self.kinetic_energy[:20000], color = 'blue', linewidth = 2 ,   label="kinetic_energy")
		plt.plot( x[:20000] , self.potential_energy[:20000],  color='red', linewidth=2 , label="potential_energy")
		plt.plot( x[:20000], self.total_energy[:20000],   color='black', linewidth=2, label="total_energy")
		plt.xlabel('Simulation Time')
		plt.ylabel('Energy (Newtons)')
		plt.title('Energy')
		plt.legend()
		plt.show()

	def plot_time(self):
		"""
		* this plots computational efficiency
		* varies per computer; not that useful ; was done for class purposes

		"""
		x = np.arange(int(self.T_stop / dt) + 2) * 28
		plt.plot( x[:-2] , self.time[:-2], color = 'blue', linewidth = 2 ,   label="Time Elapsed (seconds) ")
		plt.xlabel('Spring Evaluations')
		plt.ylabel('Elapsed Time (seconds)')
		plt.title('Computational Efficiency')
		plt.show()
