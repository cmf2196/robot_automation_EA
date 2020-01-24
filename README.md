Authors: Connor Finn, Abhinit Kothari
Date: Dec 11th, 2019

Project Overview: 
  This project involved using evoluationary algorithms to create robots which learn to move. A physics simulator was developed from scratch for evaluation of the robots. All robots were created as a combination of point masses and springs, and the robots were build by stacking tetrahedron on top of each other in a way detailed by the robot's genome.  The simulator itself modeled standard kinematics, considering gravity, friction, damping, a ground force, and spring forces between each mass. Parametric optimization techniques, such as a random search and hill climber were used as a baseline for comparison with the evolutionary algorithms. The two EA's considered for this project were a simple algorithm using crossover, mutation, and truncation selection, and a more sophisticated algorithm which used deterministtic crowding to maintain diversity in the population of robots. Both EA's outperformed the baseline algorithms, and the learning curves and diversity plots are shown in 
MECS4510_Final_report.pdf.

Files Included: 
                
                EA_library.py	
                physics_simulator.py	
                robots.py
                MECS4510_Final_report.pdf
               
                  
Model Generation Instructions:

                First run an optimization algorithm of your choosing in the main of EA_library.py
                Once complete, the robot can be visualized using the method commented out in EA_library.py              
                
Credits:  

This is a semester project for Columbia University Class MECS 4510, Evolutionary Computation and Automotive Design.
Thank you to Dr. Hod LIpson for guidance with this project, and teaching on the topic of evolutiony algorithms.

Additional Comments:

The simulator does not include collision prevention. This would be an excellent next step for improving the class.
