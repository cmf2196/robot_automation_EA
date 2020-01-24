
import numpy as np
from physics_simulator import simulator
from robots import *


# ___________________ Operations library

def crossover(r1 , r2):
	"""
	* this function takes in two robots and performs a crossover
	* the cross over will occur between the lists of (b , c) values
	* a new robot object will be generated.	
	"""	
	child_robot = tetrahedron()			# make a cube object
	child_robot.set_tetr( np.array((0 , 0 , 0)) , np.array((0 , 0 , 0)))	

	genome1 = r1.genome 	# get the two lists of (b , c) values
	genome2 = r2.genome

	length = len(genome2)				# get the length of the genome

	cut1 = np.random.randint(length)	# make two unique cuts
	cut2 = np.random.randint(length)

	while cut1 == cut2:					# make sure the cuts are different
		cut2 = np.random.randint(length)

	first = min(cut1 , cut2)			# get the smaller and larget cuts
	second = max(cut1 , cut2)

	new_list = np.vstack((genome1[:first] , genome2[first : second] , genome1[second :]))	# make new list with cuts
	child_robot.genome = new_list			# set as the new list

	child_robot.build_robot()				# actuate our list to the springs
	# center the robot about the origin
	mass_list = child_robot.mass_list
	com = child_robot.get_com(mass_list)
	child_robot.center_object(com)

	return child_robot					# return the robot with the crossed over list


def make_population(num_robots , num_tetr_added):
	"""
	* this function makes a population of robots 
	* it takes the population size as an input
	* it returns a list of these robots
	"""	
	num_robots = int(num_robots)					# these need to be integers
	num_tetr_added = int(num_tetr_added)

	robot_list = [0] * num_robots					 # empty list
	for i in range(num_robots):
		robot = tetrahedron()
		robot.set_tetr( np.array((0 , 0 , 0)) , np.array((0 , 0 , 0)))			# make the seed
		robot.set_genome(num_tetr_added)							# build the genome
		robot.build_robot()			# builds the robot and sets it about (0, val , 0)
		# center the robot about the origin
		mass_list = robot.mass_list
		com = robot.get_com(mass_list)
		robot.center_object(com)

		robot_list[i] = robot 		# fill list with breathing robots
	
	return robot_list


def calc_fitness( r , sim):
	"""	
	* this function calculates the fitness of the robot
	* the fitness is defined by how far it travels in 3 robot cycles
	* this function takes in a robot object (robots.py) and a simulator object (physics_simulator.py)

	"""

	sim.set_robot(r)						# set the robot in the simulator
	com = r.get_com(r.mass_list)			# center the object
	r.center_object(com)					# 		.... 
	com1 = r.get_com(r.mass_list)
	
	sim.simulate(.3 , 0)					# run simulator for three robot cycles | simulate(time, visualization_key)
	com2 = r.get_com(r.mass_list)			# see how far it traveled 
	distance = com2 - com1					# 		....
	horizontal = np.delete(distance , 1)	# 		....	
	fit = np.linalg.norm(horizontal)		# 		....
	
	sim.remove_robot(r)						# remove the robot from the simulator object


				# remove the robot from the simulator

	return fit

# this bubble sort function was modeled off of " https://www.geeksforgeeks.org/bubble-sort/"
def bubbleSort(population , fitness_list): 
    n = len(population) 
  
    # Traverse through all array elements 
    for i in range(n): 
  
        # Last i elements are already in place 
        for j in range(0, n-i-1): 
  
            # traverse the array from 0 to n-i-1 
            # Swap if the element found is greater 
            # than the next element 
            if fitness_list[j] < fitness_list[j+1] : 
                fitness_list[j], fitness_list[j+1] = fitness_list[j+1], fitness_list[j] 
                population[j] , population[j+1] = population[j+1] , population[j]





def diversity_val(population):
	"""
	* this function takes in a population and returns an overall diversity score by evaluating every combination of robots
	* the diversity value between two robots is decided by the function div_distance
	* the diversity score has little stand alone meaning, but is good for comparison between populations
	"""
	divScore = 0										# init the diversity score
	for i in range(len(population)):					# iterate over every possible pairing
		for j in range(len(population)):				# 			....
			if i ==j:									# compare the genomic elements
				pass 									# if the same, do nothing
			else:										# if they are different, add the distance betweeen the numbers
				divScore = divScore + div_distance(population[i] ,population[j])	
	return divScore

def div_distance(r1 , r2):
	"""
	* this function generates a score representing the genomic diversity of two robots.
	* b is a value in  (0 , .03) , c is a value in (0 , 2 * np.pi).
	* because the variation will be larger for c we cannot use (c2-c1) + (b2-b1) as a determinant for diversity bc diff in c will outweigh diff in b
	* instead use  a normalized version (b2 - b1) / .03 + (c2 - c1) / 2*pi

	"""

	diversity_score = 0
	array1 = r1.genome				# get the genome of robot 1
	array2 = r2.genome				# get the genome of robot 2
	length = np.size(array1 , 0)	# get the number of columns ( genome is an array of values (b,c))

	for i in range(length):
		b1 = array1[i][0]					# get the shape ( b) and type (c ) vals for each robot
		c1 = array1[i][1]					#               ....
		b2 = array2[i][0]					#               ....
		c2 = array2[i][1]					#               ....

		val = abs(b2 -  b1) / 0.03 + abs(c2 - c1) / (2 *np.pi)		# get value of interest for each spring
		diversity_score += val 									# update the diversity score


	return diversity_score



# _______________ Optimization algorithms 


def hill_climber():
	"""
	* this function performs a parametric optimization hill climber algorithm on a single robot.
	* it will save the top performing robot, as well as learning curve 

	"""
	sim = simulator()											# create the simulator for this algorithm
	num_climb = 9999											# the totol number of robot evaluations is num_climb + 1
	num_tetr_added = 14											# this is the size of our robots
	hc_fit = np.zeros(num_climb+1)								# create an array to store the learning curve data
	# keep track of progress
	best_fit = 0												# init the best fit as zero
	best_array = np.empty(8)									# init a random array as the best robot
	# make a random robot	
	robot = tetrahedron()													# make the seed 
	robot.set_tetr( np.array((0 , 0 , 0)) , np.array((0 , 0 , 0)))			# 		....
	robot.set_genome(num_tetr_added)										# generate the genome
	robot.build_robot()														# build the robot 
	
	mass_list = robot.mass_list					# get the robot masses
	com = robot.get_com(mass_list)				# center the robot
	robot.center_object(com)					# 		....

	best_fit = calc_fitness( robot , sim)		# get the fitness of the first robot
	hc_fit[0] = best_fit						# place it as the first value in the learning curve array

	for i in range(1 ,num_climb + 1):			# perform mutations equal to num_Climb
		print('mut' , i)						
		mut_loc , old_val = robot.mutate_genome()		# Mutation:  Keep track of mut location and previous vals
												# rebuild the robot with the new (mutated) genome
		robot.remove_old()						# remove old springs and masses
		robot.set_tetr( np.array((0 , 0 , 0)) , np.array((0 , 0 , 0)))		# set the robot to pos(0 , 0 , 0)
		robot.build_robot()						# builds the robot
												
		mass_list = robot.mass_list				# center the robot
		com = robot.get_com(mass_list)			# 		....
		robot.center_object(com)				# 		....

		fit_new = calc_fitness( robot , sim)	# calculate new fitness
		if fit_new > best_fit:					# if better update the top fitness, and top performing genome
			best_fit = fit_new					
			best_array = np.copy(robot.genome)

		else:									# if worse, revert the genome back to its previous state. 
			robot.genome[mut_loc[0] , mut_loc[1]] = old_val				# only need to adjust genome
		hc_fit[i] = best_fit											# robot is rebuilt prior to fitness calculation each time
		
	np.savetxt("hc2_robot.csv", best_array, delimiter=",")
	np.savetxt("hc2_learning.csv", hc_fit, delimiter=",")



def random_search():
	""" 
	* this function will perform a random search for the best robot
	* this is done by generating a new random robot, comparing it to the current best, and keeping it only if it achieves higher fitness
	"""		
	sim = simulator()							# simulator object, to be used during this algorithm
	num_search = 20								# number of robots evaluated
	num_tetr_added = 14							# this is the size of our robots ( 15 tetrahedron total )

	rs_fit = np.zeros(num_search)				# learning curve array
	best_fit = 0								# best fitness place holder
	best_array = np.empty(8)					# top performing genome place holder
	
	for i in range(num_search):					# do the random search	
		print('counter: ', i)

		robot = tetrahedron()													# make a random robot
		robot.set_tetr( np.array((0 , 0 , 0)) , np.array((0 , 0 , 0)))			# make the seed
		robot.set_genome(num_tetr_added)										# generate a random genome
		robot.build_robot()														# builds the robot 
		
		mass_list = robot.mass_list				# center the robot
		com = robot.get_com(mass_list)			# 		....
		robot.center_object(com)				# 		....

		fitness = calc_fitness( robot , sim)	# evaluate the robot's fitness

		if fitness > best_fit:					# if the fitness is better than current best
			best_fit = fitness 					# update placeholders
			best_array = np.copy(robot.genome)	# 		....

		rs_fit[i] = best_fit					# add the best fit to the learning curve

	np.savetxt("rs_robot.csv", best_array, delimiter=",")
	np.savetxt("rs_learning.csv", rs_fit, delimiter=",")

def simple_evolution():
	"""
	* this function takes in an initial population of robots in the form of a list
	* it then performs evolutionary algorithm consisting of simple mutations and crossovers
	"""
	sim = simulator()							# create the simulator object to be used with this algorithm
	num_tetr_added = 14
	population_size = 4							# how many robots 
	num_mutations = 4							# mutations between generation
	num_generations = 4							# how many generations

	evals = population_size * num_generations * num_mutations		# this is the number of robot evaluations
	diversity_vals = np.zeros(num_generations)						# diversity score array

	overall_fitness = np.zeros(evals + 1)		# learning curve array
	ticker = 0									# keeps track of which robot eval we are on
	parent_size = int(population_size / 2)		# number of robots kept for crossover ( using truncation selection so top 50% are kept )

	parents = make_population(parent_size , num_tetr_added )		# make the population of parents
	population_fitness = [0] * population_size 						# array used to keep track of the fitness of each robot in the population
	
	for m in range(num_generations):			# perform the EA
		print('generation: ' , m)

# _______  crossover ________

		population = parents.copy() 							# get an array for the entire population ( needs to have parent robots in it )
		for k in range(parent_size):							# select 2 random parents
			p1 = np.random.randint(population_size/2)			#		....
			p2 = np.random.randint(population_size/2)			#		....
											
			while p2 == p1:										# make sure they are not the same robot					
				p2 = np.random.randint(population_size/2)		# note: if pop size = 2; stuck here forever
			
			pr1 = parents[p1]									# perform a 2 cut crossover
			pr2 = parents[p2]									# 		....
			c = crossover(pr1 , pr2)							# 		....
			population = population + [c]						# add the child robot to the population

		for z in range(len(population_fitness)):				# update the population fitness array 
			if population_fitness[z] ==0:						# only do this if there are zeros 
				robot = population[z]														# select the robot
				robot.remove_old()															# build according to genome
				robot.set_tetr( np.array((0 , 0 , 0)) , np.array((0 , 0 , 0)))				# 		....
				robot.build_robot()															# 		....
				
				mass_list = robot.mass_list						# center the robot about the origin
				com = robot.get_com(mass_list)					# 			....
				robot.center_object(com)						# 			....

				population_fitness[z] = calc_fitness(robot , sim)		# calc fitness and add to array 
	
# ______ Mutation _________

		diversity_vals[m] = diversity_val(population)			# get population diversity before mutation | add to diversity array
		
		for j in range(num_mutations):
			print('mut' , j)
			for p in range(population_size):					# mutate each robot in the population
				robot = population[p]							#			....

				if ticker == 0:									# if this is the first evaluation, add fitness as first element to learning curve array
					overall_fitness[ticker] = population_fitness[p]			#			....
					ticker += 1												#Keep track of evaluation number
				
				mut_loc , old_val = robot.mutate_genome()		# Mutation:  Keep track of mut location and previous vals
																# rebuild the robot with the mutated genome
				robot.remove_old()								# remove old springs and masses
				robot.set_tetr( np.array((0 , 0 , 0)) , np.array((0 , 0 , 0)))			# set coordinates and initial velocity
				robot.build_robot()								# builds the robot
				
				mass_list = robot.mass_list						# center the robot about the origin
				com = robot.get_com(mass_list)					# 			....
				robot.center_object(com)						# 			....

				fit_new = calc_fitness( robot , sim )			# get the robot's fintess
					
				if fit_new > population_fitness[p]:				# if higher fitness then this robot's previous best : update
					population_fitness[p]= fit_new				# 			....
					
				else:
					robot.genome[mut_loc[0] , mut_loc[1]] = old_val			# if lower, revert the genome to it's previous state
				
				if overall_fitness[ticker - 1] < population_fitness[p]:		# if the best overall robot thus far
					overall_fitness[ticker] = population_fitness[p]			# add this as next value for the learning curve
					best_array = np.copy(robot.genome)						# update the best robot's genome
				else:														# if lower, push the previous best score forward in learning curve
					overall_fitness[ticker] = overall_fitness[ticker - 1]
				ticker +=1													# end of mutation, add one to robot evaluation ticker

# __________ Selection ____________

		bubbleSort(population , population_fitness)			# sort the population according to fitness rank ( sorts both population and population_fintess in tandem)

		parents = population[: parent_size] 				# perform truncation selection (choose top half)
		second_half = [0] * parent_size						# 	array of zeros
		population_fitness = population_fitness[: parent_size] + second_half		# 	 keep first half (second half are zeros)
		

	best_robot = population[0]			# extract the best performing robot
	a = best_robot.genome 				# get it's genome

	np.savetxt("ea2_robot.csv", a, delimiter=",")
	np.savetxt("ea2_learning.csv", overall_fitness, delimiter=",")
	np.savetxt("ea2_diversity.csv", diversity_vals, delimiter=",")





def diverse_evolution():
	"""
	* this function takes in an initial population of robots in the form of a list
	* it then performs evolutionary algorithm consisting of simple mutations and crossovers
	* this evolutionary algorithm is enhanced by adding diversity maintenance ( deterministic Crowding )
	"""
	sim = simulator()							# create physics simulator object
	num_tetr_added = 14							# robots are constructed with 15 tetrahedron
	population_size = 50						# how many robots 
	num_mutations = 25							# mutations between generation
	num_generations = 40						# how many generations

	evals = population_size * num_generations * num_mutations		# number of robot evaluations
	overall_fitness = np.zeros(evals + 1)							# learning curve array
	diversity_vals = np.zeros(num_generations)						# diversity score array
	
	ticker = 0													# keeps track of which robot eval we are on
	parent_size = int(population_size / 2)						# half the population will be kept for crossover
	parents = make_population(parent_size , num_tetr_added)		# make the population
	
	population_fitness = [0] * population_size  				# init fitness score evals array ( kept for each robot in the populaiton )
	
	for m in range(num_generations):							# start the EA
		print('generation: ' , m)

# __________ crossover ____________: 

		population = parents.copy() 
		for k in range(parent_size):
			p1 = k												# want to ensure, each robot is the primary parent for a child robot
			p2 = np.random.randint(population_size/2)			# select a random second robot
			
			while p2 == p1:										# make sure robot's are different
				p2 = np.random.randint(population_size/2)		# note: if pop size = 2; stuck here forever
			
			pr1 = parents[p1]					# make a two cut crossover 
			pr2 = parents[p2]					# 			....
			c = crossover(pr1 , pr2)			# 			....
			population = population + [c] 		# add the child robot to population

		# update the population fitness array ( only need to if there are zeros)
		#print(" before crossover" , population_fitness)
		for z in range(len(population_fitness)):
			if population_fitness[z] ==0:
				robot = population[z]					# select the robot where the fitness is zero
				
				# clear and rebuild this robot
				robot.remove_old()
				robot.set_tetr( np.array((0 , 0 , 0)) , np.array((0 , 0 , 0)))			# build seed
				robot.build_robot()														# build the robot according to its current genome
																						
				mass_list = robot.mass_list												# center the robot 				
				com = robot.get_com(mass_list)											# 		"  "
				robot.center_object(com)												# 		"  "
				population_fitness[z] = calc_fitness(robot , sim)						# find new fitness

# __________ Mutation ____________
		diversity_vals[m] = diversity_val(population)									# get population diversity before mutation
		for j in range(num_mutations):
			print('mut' , j)							
			for p in range(population_size):									# mutate each robot in the population
				robot = population[p]											#				"  "
		
				if ticker == 0:													# when this is the first robot
					overall_fitness[ticker] = population_fitness[p]				# set the best fitness to this score
					ticker += 1
				
				mut_loc , old_val = robot.mutate_genome()								# mutate the genome
				
				robot.remove_old()														# clear robor's springs and point masses
				robot.set_tetr( np.array((0 , 0 , 0)) , np.array((0 , 0 , 0)))			# build the robot according to new genome
				robot.build_robot()														#				"  "
			
				mass_list = robot.mass_list										# center the robot 
				com = robot.get_com(mass_list)									#		"  "
				robot.center_object(com)										#		"  "
				fit_new = calc_fitness(robot , sim)								# find new fitness
	
				if fit_new > population_fitness[p]:								# better fitness than this robot previously 
					population_fitness[p]= fit_new								# update

				else:
					robot.genome[mut_loc[0] , mut_loc[1]] = old_val				# revert genome back to previous state
															
				if overall_fitness[ticker - 1] < population_fitness[p]:			# if best overall fitness is found
					overall_fitness[ticker] = population_fitness[p]				# update learning curve vals
					best_array = np.copy(robot.genome)							# save the best genome
				
				else:															# if not, the best fitness at this # of evals was the previous one
					overall_fitness[ticker] = overall_fitness[ticker - 1]
				ticker +=1														# add to the number of evaluations counter

# __________ Selection ____________
		# use a selection process with deterministic crowding

		for i in range(parent_size):
			if population_fitness[parent_size + i] > population_fitness[i]:			# compare child to primary parent
				parents[i] = population[parent_size +i]								# if child is better, replace parent 
				population_fitness[i] = population_fitness[parent_size +i]			#				"  "
				population[i] =  population[parent_size +i]							#				"  "
			
			population_fitness[parent_size +i] = 0									# set fitness of second half to be zero (we create these with crossover)

	bubbleSort(population , population_fitness)								# this just to keep track of progress	
	best_robot = population[0]
	a = best_robot.genome 													# return the best robot genome


	# store locally
	np.savetxt("div_ea_3k_best_robot2.csv", a, delimiter=",")
	np.savetxt("div_ea_3k_learning2.csv", overall_fitness, delimiter=",")
	np.savetxt("div_ea_3k_diversity2.csv", diversity_vals, delimiter=",")



def main():


# RUN PROGRAMS

	#simple_evolution()
	#diverse_evolution()
	#hill_climber()
	#random_search()
	
# RUN SIMULATION OF 2 BUILT ROBOTS, STARTING AT THE ORIGIN

	r1 = tetrahedron()
	r1.set_tetr(np.array((0 , 0 , 0)) , np.array((0 , 1.4 * np.pi, 0)))
	print(len(r1.spring_list))
	# r1.genome = np.genfromtxt(PATH TO GENOME, delimiter=',')
	# r1.build_robot()

	# r2 = tetrahedron()
	# r2.set_tetr(np.array((0 , 0 , 0)) , np.array((0 , 1.4 * np.pi, 0)))
	# r2.genome = np.genfromtxt(PATH TO GENOME, delimiter=',')
	# r2.build_robot()


	# sim = simulator()
	# sim.set_robot(r1)
	# sim.simulate(.5 , 2)
	# sim.remove_robot(r1)
	# sim.reset()
	# sim.set_robot(r2)
	# sim.simulate(.5 , 2)	


main()