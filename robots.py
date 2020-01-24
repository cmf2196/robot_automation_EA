"""
This script bounces a cube off of a floor using a custom built physics simulator and v python


"""
import numpy as np
from vpython import *


## constants we will be using

K = 10000				# spring constant N/M 
musc_1_const = 11000	# muscle type 1
musc_2_const = 9000		# muscle type 2
tissue_const = 500		# soft support
bone_const = 1000000	# hard support
tall = 0.08164965809277261




			
class point_mass:
    
    def __init__(self, mass , state ):
        self.mass = mass 					# flaot in kg
        self.state = state					# state inclues position, velocity, and accelerate ( 3 by 3 np array)
        self.color = vector(0, 0, 0) 

class spring:
	def __init__(self,  m1 , m2):
		self.spring_constant = K 						#
		self.rest_length = [m1.state[0] - m2.state[0] , 0 , 0]				# this is an array [Lo, b , c] feed into eq Lo + b * sin(wt + c)
		self.mass_1 = m1 							# mass object
		self.mass_2 = m2  							# mass object
		self.color = vector(0, 0, 0)				# use http://www.rgbtool.com/ to select color for this

"""
 * the following will all be subclasses of the spring class
 * they represent the types of springs that a robot could have
"""

class muscle_1(spring):
	
	def __init__(self, m1 , m2):
		super().__init__(m1 , m2)								# get typical spring initiation
		self.color = vector(0 , 252 , 30)						# change the color
		self.rest_length[1] = 0.025								# change the b value
		self.rest_length[2] = 0									# change the c value
		self.spring_constant = musc_1_const
		


class muscle_2(spring):
	
	def __init__(self, m1 , m2):
		super().__init__(m1 , m2)							# get typical spring initiation
		self.color = vector(252, 0 , 30) 					
		self.rest_length[1] = 0.01							# change the b value
		self.rest_length[2] = 1.25 * np.pi					# change the c value
		self.spring_constant = musc_2_const


class tissue(spring):

	def __init__(self, m1 , m2):
		super().__init__(m1 , m2)							# get typical spring initiation
		self.color = vector(252, 252 , 0)					# Yellow
		self.rest_length[1] = 0								# change the b value
		self.rest_length[2] = 0								# change the c value
		self.spring_constant = tissue_const

class bone(spring):
	
	def __init__(self, m1 , m2):
		super().__init__(m1 , m2)							# get typical spring initiation
		self.color = vector(252, 252 , 252)
		self.rest_length[1] = 0  							# change the b value
		self.rest_length[2] = 0								# change the c value
		self.spring_constant = bone_const


class robot:
	"""
	* this is our Parent Class

	"""

	def __init__(self):
		self.mass_list = [0] * self.masses_num				# masses_num , spring_num set in subclass
		self.spring_list = [0] * self.spring_num
		self.center_of_mass = 0
		self.breathing_list = np.zeros(self.spring_num)		# list of rest lengths


	
	def get_com(self, mass_list):
		""" 
		* this function takes in a mass list
		* it returns the com

		"""
		num_points = len(mass_list)		# get the number of masses
		pos_sum = np.zeros(3)			# place holder
		
		for mass in mass_list:			# add to the place holder
			pos = mass.state[0] 
			pos_sum +=pos 
		com = pos_sum / num_points		# divide by number of masses

		return com


	def center_object(self, com):
		"""
		* this function takes in the com of the object, and moves every mass so that the com is hovering about (0 , val , 0)

		"""
		movement = -1* com 		# want to get com to (0 , val , 0)
		movement[1] = 0			# dont change height

		masses = self.mass_list
		for mass in masses:					# move all masses the necessary shift
			mass.state[0] += movement

	def set_breathing_aray(self):
		"""
		* this function sets the spring values to those in self.breathing_list
		* self.breathing_list is an array that contains elements ( b , c) in the formula rest_length = Lo + b * sin(w*T + c)

		"""
		num_springs = self.spring_num										

		for i in range(num_springs):												# set the values of breathing array equal to the inputs
			self.spring_list[i].rest_length[1:] = self.breathing_list[i]			# the first value of rest lenght does not change

	def make_breathing_array(self):
		"""
		* this function takes in a robot and creates a new list of (b , c) values (in arrays)
		* The array will be used for the breath function
		* by definition of the period for a sin function c is in (0 , 2pi)
		* b can be as large or small as we want. Let's just make it a random percentage of the rest length Lo so b in (0 , 1)
		"""
		num_springs = self.spring_num 									# get spring number
		breath_list = [0] * num_springs			# make an empty array of that size
		for i in range(num_springs):									# fill empty array with random values from apropriate domains
			b = np.random.uniform(0 , .03) 								# 					"  "
			c = np.random.uniform(0 , 2 * np.pi)						# 					"  "
			breath_list[i] = np.array((b,c))							# 					"  "

		self.breathing_list = breath_list								# Set the generated array


	def mutate(self):
		""" 
		* this function makes a small mutation to our breathing cube
		* the mutation is a simple adjustment to one of the springs.
		"""

		num_springs = self.spring_num
		spring_index = np.random.randint(num_springs)			# select a random spring index

		spring_temp = self.spring_list[spring_index]			# get the spring
		mutate_items = spring_temp.rest_length[1:]				# get the things we will mutate

		length = len(mutate_items)								# this is a list

		pick_one = np.random.randint(length)					# pick either b , c ect..
		val = mutate_items[pick_one]

		adjust = np.random.uniform(0.99 , 1.01) 				# get adjustment 
		spring_temp.rest_length[pick_one + 1] = val * adjust 	# update the spring

		return spring_index										# this will be useful for updating self.breathing_array


	def get_diameter(self , mass_list):
		""" 
		* this function gets the size of the smalles possible diameter that encases the robot
		* useful for evaluating the robots motion
		"""
		
		di = 0
		for mass in mass_list:						# iterate over all the masses
			for i in range(len(mass_list)):			# get all combinations
				alt = mass_list[i]					
				pos_sum = np.zeros(3)				# place holder
				
				pos1 = mass.state[0]				# get distance between two masses
				pos = alt.state[0]					# 		"  "
				pos_sum = pos - pos1				# 		"  "
				dist = np.linalg.norm(pos_sum)		# 		"  "
				
				if dist > di:						# Get the maximum distance 
					di = dist

		return di



class cube(robot):

	"""
	* this is a subclass of robot
	* the cube was used a unit test for the physics simulator

	"""
	
	def __init__(self):
		self.spring_num = 28		
		self.masses_num = 8
		super().__init__()
		self.type = 0					# this is a pointer for use in a switch statement 



	def set_cube(self , location , theta , spin ):
		"""  
	 * This method will make the mass and spring objects we need for a cube
	 * location is a 1 by 3 numpy array, it gives the adjustmetn from 0 , 0
	 * theta is the angle rotated in radians
	 * this cube will not have breathing components.

		"""

		# build the transformation matrix. This will set the cube to the euler angles provided.
		theta_x = theta[0]
		theta_y = theta[1]
		theta_z = theta[2]

		# this are definitions from intro to robotics course

		trans_x = np.array(((1 , 0 , 0 , location[0] ), (0 , np.cos(theta_x) , -1 * np.sin(theta_x) , 0) , (0 , np.sin(theta_x) , np.cos(theta_x) , 0) ,(0 , 0 , 0 , 1)) )
		trans_y = np.array(((np.cos(theta_y) , 0 , np.sin(theta_y) , 0 ), (0 , 1 , 0 , location[1]) , (-1 * np.sin(theta_y) , 0 , np.cos(theta_y) , 0) ,(0 , 0 , 0 , 1)) )
		trans_z = np.array(((np.cos(theta_z) , -1 * np.sin(theta_z) , 0 , 0 ), (np.sin(theta_z) ,  np.cos(theta_z), 0 , 0) , (0 , 0 , 1 , location[2]) ,(0 , 0 , 0 , 1)) )

		# the tranformation matrix of interest is the dot product of the x , y , z transformation matices
		transformation_matrix = np.dot(trans_x , np.dot(trans_y , trans_z))

# _______ Build Masses ________
		index = 0															  		 # used to keep track of loc
		
		for p in range(2):																				# build the robot's position
			for l in range(2):
				for n in range(2):
										
					pos = np.array((p / 10 -.05 , l/ 10 -.05  , n / 10 -.05 ))							# center the cube at origin
					pos_m = np.hstack((pos , 1)) 														# build the position vector (in R^4)
					pos_m_alt = np.dot( transformation_matrix , pos_m)									# dot product with transformation matrix
					pos_alt = pos_m_alt[0:3]															# This is the desired location
					
					if spin:																			# we can initiate the cube with a rotation
						if p == 0:																		# add a motiton to the two opposite faces
							kinematics = np.vstack((pos_alt , np.array((0 , 1 , 0)) , np.zeros(3)))		# 				"  "
							m = point_mass(.4 , kinematics)         									# 	build point mass
						if p ==1:																		# 				"  "
							kinematics = np.vstack((pos_alt , np.array((0 , 0 , 0)) , np.zeros(3)))		# 				"  "
							m = point_mass(.4 , kinematics)         									# 	build point mass
																										# note (8 masses: all acocunted for here)
					else:
						kinematics = np.vstack((pos_alt , np.zeros(3) , np.zeros(3)))					# if no spin, set initial velocity and accell to 0
						m = point_mass(.4 , kinematics)         										# build point  mass
					self.mass_list[index] = m 															# add to robot's mass list
					index += 1																			# increase the counter

# _______ Build Springs ________
		num = self.masses_num
		count = 0

		for j in range(num):								
			for i in range(j + 1 , num):					# Build a spring between every combo of masses. 
				m1 = self.mass_list[i] 						# not concerned about interfearance
				m2 = self.mass_list[j]
															
				if self.type == 0:							# indicate which type of material will be used
					s = spring( m1 , m2)					

				elif self.type == 1:
					s =muscle_1(m1 ,m2)

				elif self.type == 2:
					s = muscle_2(m1 ,m2)

				elif self.type == 3:
					s = bone(m1 , m2)

				elif self.type ==4:
					s = tissue(m1 , m2)

				self.spring_list[count] = s 				# add spring to spring_list
				count += 1




class tetrahedron(robot):

	"""
	* This subclass starts with a seed being a tetrahedron
	* more complex robots are build off of this seed

	"""

	
	def __init__(self):
		self.spring_num = 6
		self.masses_num = 4
		super().__init__()
		self.type = 0					# this is a pointer for use in a switch statement 
		self.sides = [] 				# this list tells which sides are availible to build off 
		self.genome = [] 				# this is used for the adding tetrahedron function



	def set_genome(self , num_tetr_added):
		"""
		* this function builds a random genome
		* the genome is a 2 by n ( i used 14) array 
		* It indicates the side a new tetrahedron will be added to, and what material will be used
		"""

		genome = np.zeros((num_tetr_added , 2))					# first column is side to add tetr to, 
																# second column is what type of spring
		for i in range(num_tetr_added):
			side_select = np.random.randint(2 * i + 3)			# start with 3 sides, add two more each time ( add three subtract 1)
			type_select = np.random.randint(4)					# bone, tissue, muscle 1 , muscle 2
			genome[i] = (side_select , type_select)

		self.genome = genome



	def mutate_genome(self):
		"""
		* this function makes a small mutation to the genome
		* it will return the indices of the adjustment - that way we can keep or revert the change after evaluating the fitness
		"""

		genome = self.genome 								# get the genome and its dimensions
		num_rows = np.size(genome , 1)						# 			"  "
		num_cols = np.size(genome ,0)						# 			"  "

		column = np.random.randint(num_cols)				# get a random position in the genome
		row = np.random.randint(num_rows)					# 			"  "
		old_val = genome[column , row]						# 			"  "

		if row ==0:											# reselecting the side to build off
			val = np.random.randint(column * 2 + 3)			# number of available sides

		else:												# reselecting the material for this item
			val = np.random.randint(4)						# only have four materials 
	
		genome[column , row] = val 							# update the genome

		return [column , row] , old_val						# return the indices of the genome that we mutated 
															# this way we can revert back if mutation was bad


	def remove_old(self):
		"""
		* this function removes all the old masses springs and sides on a robot
		* used after a genome mutation to rebuild the robot 
		"""

		self.mass_list = [0] * self.masses_num			
		self.spring_list = [0] * self.spring_num		
		self.sides = []			


	def build_robot(self):	
		""" 
		* this function is ran after set_genome
		* it will add the tetrahedron to the seed robot in accordance to the genome
		* the function then ensures no masses start below zero
		* finally it sets the center of mass to have a 0 x and z value
		"""
		genome = self.genome 							# collect the genome
		
		for i in range(np.size(genome , 0)):			# go through the genome
			v1 = int(genome[i , 0 ])					# get the side to build on
			v2 = int(genome[i , 1])						# get the material to be used
			self.add_tetr(v1 , v2)						# build the tetrahedron on top of the current robot

		min_y = 0										# need to adjust so no masses start below the floor
		for mass in self.mass_list:						# check all masses
			if mass.state[0][1] < min_y:
				min_y = mass.state[0][1] 				# calculate the lowest point on the robot
		
		if min_y < 0:									# if the robot's lowest point is below zero, move everything up
			for mass in self.mass_list:
				mass.state[0][1] += -1* min_y

	def set_tetr(self , location , theta ):
		"""  
	 * This method will make the mass and spring objects we need for a tetrahedron
	 * location is a 1 by 3 numpy array, it gives the adjustmetn from 0 , 0
	 * theta is the angle rotated in radians
		"""

# ____ Build Transformation Matrix _____
"""
* This is done in accordance to definition from intro to robotics course
"""

		# get the rotational components
		theta_x = theta[0]
		theta_y = theta[1]
		theta_z = theta[2]
		# build for transfomraiton about x , y , z axis
		trans_x = np.array(((1 , 0 , 0 , location[0] ), (0 , np.cos(theta_x) , -1 * np.sin(theta_x) , 0) , (0 , np.sin(theta_x) , np.cos(theta_x) , 0) ,(0 , 0 , 0 , 1)) )
		trans_y = np.array(((np.cos(theta_y) , 0 , np.sin(theta_y) , 0 ), (0 , 1 , 0 , location[1]) , (-1 * np.sin(theta_y) , 0 , np.cos(theta_y) , 0) ,(0 , 0 , 0 , 1)) )
		trans_z = np.array(((np.cos(theta_z) , -1 * np.sin(theta_z) , 0 , 0 ), (np.sin(theta_z) ,  np.cos(theta_z), 0 , 0) , (0 , 0 , 1 , location[2]) ,(0 , 0 , 0 , 1)) )
		# matrix of interest is the dot product about the principle axis
		transformation_matrix = np.dot(trans_x , np.dot(trans_y , trans_z))
		
#_______ Calculations for tetr centered at origin _____
		side_length1 = 0.1
		depth = np.sqrt(0.1**2 - 0.05**2)
		mid = np.sqrt((depth/3)**2 + 0.05**2)
		height = np.sqrt(0.1**2 - mid**2)
		mass_positions = [ np.array((- 0.5 * side_length1 , 0.0 , -depth / 3)) , np.array((0.5 * side_length1 , 0.0 , - depth / 3) ), np.array((  0.0 , 0.0 , 2*depth / 3)) , np.array((0.0, height , 0.0))]
		
# __________ Build the masses ________
		index = 0															   # used to keep track of loc
		# these nested for loops are used to fill a position vertor 
		for p in range(self.masses_num):										# go through the number of masses
			pos = mass_positions[p]
			pos_m = np.hstack((pos , 1)) 										# build position vector for each mass in R^4
			pos_m_alt = np.dot( transformation_matrix , pos_m)					# dot product with the transformaiton matrix
			pos_alt = pos_m_alt[0:3]	

			kinematics = np.vstack((pos_alt , np.zeros(3) , np.zeros(3)))		# set accell and velocity to zero
			m = point_mass(.4 , kinematics)                						# build the point mass
			self.mass_list[index] = m 											# add it to the mass list
			index += 1

# ________ Build the springs ________
		num = self.masses_num 
		count = 0

		for j in range(self.masses_num):					# make a spring connecting every combination of masses
			for i in range(j + 1 , num):
				m1 = self.mass_list[i]
				m2 = self.mass_list[j]

				s = muscle_1( m1 , m2)						# build the spring object
				self.spring_list[count] = s 				# add to the spring list
				count += 1

#_______ build side objects ______

# only concerned about the one we build on

		mass_list = self.mass_list
		# we would like to provide the center of mass as well / will help with calculating normal
		com = self.get_com(mass_list)
		s1 = side([mass_list[0] , mass_list[1] , mass_list[3]] , com)				# only three open sides. Hard code this
		s2 = side([mass_list[0] , mass_list[2] , mass_list[3]] , com)				# bottom is not opne
		s3 = side([mass_list[1] , mass_list[2] , mass_list[3]] , com)

		self.sides = [s1 , s2 , s3]																




	def get_triangle(self, points , com ):
		"""
		* this function gets the coordinates for the centroid and normal of a triangle (side)
		* it takes in a list of 3 point masses
		* (it could have just taken a side object)

		"""
		# get the x,y,z coordinates for each of the point masses
		p1 = points[0].state[0]
		p2 = points[1].state[0]
		p3 = points[2].state[0]

		# get the centroid coordinates
		Cx = (p1[0] + p2[0] + p3[0] ) / 3 
		Cy = (p1[1] + p2[1] + p3[1] ) / 3 
		Cz = (p1[2] + p2[2] + p3[2] ) / 3 

		# get two sides of the triangle
		s1 = p2 - p1
		s2 = p3 - p1

		# take the cross product
		n = np.cross(s1 , s2)
		size = np.linalg.norm(n)
		n = n / size

# build the normal vector
		# need to make sure direction is ok.  
		if Cx > com[0]:						# want to go away from the tetrahedron's com
			n[0] = abs(n[0])
		else:
			n[0] = -1 * abs(n[0])

		if Cy > com[1]:						# want to go away from the tetrahedron's com
			n[1] = abs(n[1])
		else:
			n[1] = -1 * abs(n[1])

		if Cz > com[2]:						# want to go away from the tetrahedron's com
			n[2] = abs(n[2])
		else:
			n[2] = -1 * abs(n[2])

		return Cx , Cy , Cz , n 					# return centroid and normal vecotr



	def add_tetr(self , side_num , type_spring):
		"""
		* This function adds a tetrahedron to our current robot
		* side_num: index for the side to build onto
		* type_spring: pointer indicating which material type will be used
		"""

		side1 = self.sides[side_num]				# get the side
		masses = side1.mass_list					# get the masses
		com = side1.com 								
		Cx , Cy , Cz , n = self.get_triangle(masses , com)		# get centroid and normal vector

# _______ Add new mass _____
		height = tall * n 			### IMPORTANT ###  This height is empirically determined: NEED to be more exact if we want to code collision prevention later
		height += np.array((Cx , Cy , Cz))		# start at centroid
		
		kinematics = np.vstack((height , np.zeros(3) , np.zeros(3))) 	# get the mass' state 
		m = point_mass(.4 , kinematics)									# build point mass
		self.mass_list += [m]											# add mass to our robot

		for i in range(len(masses)):									# build springs from the three point masses on the side to the new mass
			m1 = masses[i]
			
			if type_spring == 0:										# chose the indicated material type
				s = tissue(m1 , m)	

			elif type_spring == 1:										# 			"  "
				s =muscle_1(m1 ,m)

			elif type_spring == 2:										# 			"  "
				s = muscle_2(m1 ,m)

			elif type_spring == 3:										# 			"  "
				s = bone(m1 , m)
			
			self.spring_list += [s]										# add the spring to the robot

# _______ Update availible sides _____
		
		self.sides.pop(side_num)										# remove the currently used side							

		masses +=[m]													# get the com for new tetrahedron												
		new_com = self.get_com(masses)

		for mass in masses:												# Asthetics
			mass.color = s.color

		s1 = side([masses[0] , masses[1] , m ] , new_com) 				# build the three new sides (again, can simply hard code)
		s2 = side([masses[0] , masses[2] , m ], new_com)
		s3 = side([masses[1] , masses[2] , m ] , new_com)

		self.sides += [s1 , s2 , s3]									# add to robnot


class side:									
	"""	
	* new class for adding tetrahedron
	* Typcally a side is defined by three points and a normal vecor
	* here it is three points and the COM of the base tetrahedron
	* rectify this by using the get_triangle function to calculate the normal vector
	"""

	def __init__(self , masses , com):		# takes in the list of masses, as well as the center of mass
		self.mass_list = masses
		self.com = com 						# this will help calculate the normal ( com of the tetrahedron attatched)
		self.spring_list = []

