#Copyright © 2019 ALEXANDER GROH
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#This program produces a phase transition plot, using a function Simulate()
    # the parameters of the plot (region of sigma and steps in sigma) can be changed
    # in section #parameters
# before using this program for the first time, we have to set a visual inspectionn variable,
# which can be obtained from visual_inspection.py
# initially visual_inspection = 20

import numpy as np
from numpy import linalg as la
import numpy.random as rd
import random
import scipy
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as ani




#parameters
N = 40          # number of particles/bird
L_x = 3.1       # boundary in x
L_y = 3.1       # boundary in y
R = 1           # radius of birds' control
v = 0.1         # speed of each particle/bird
sigma_min = 0.1   # starting standard deviation (sigma)
sigma_max = 4     # ending sigma
sigmaStep = 0.2   #step in sigma
steps = 50  # steps to complete
visual_inspection = 20  #after how many steps does order parameter converge
                            #can be changed if the plot gives different result

#generate_positions, given size x,y it generates uniform distribution of n elements on two dimensional grid
#returns an array of position vectors for each particle
def generate_positions(n, x, y):
    positions = np.zeros((2,n)) #positions is an array of 2D vectors
                                #(each vector corresponds to position of one particle)
    for i in range(len(positions[0,:])):
        positions[0,i] = random.uniform(0,x)    #using a uniform random generation
        positions[1,i] = random.uniform(0,y)
    return positions

#generates a random unit vector
#random unit vector (in all directions):
def random_unit_vector():
    x = random.uniform(-1.0,1.0)
    y = np.sqrt(1-x**2)*random.choice([-1,1])
    unit_vector = [x,y]
    return unit_vector

#returns distance of two vectors
def distance(a,b):
    dist = np.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2)
    return dist

#finds the distance between two 2D vectors in a periodic grid
#returns the smallest distance between equivalent point of the grid
def distance_periodic(a,b):
    x = a[0]
    y = a[1]
    distance_array = np.zeros((1,9))
    equivalent_a = np.array([[x,x,x+L_x,x+L_x,x+L_x, x, x-L_x, x-L_x, x-L_x ],[y, y + L_y,y+L_y,y,y-L_y, y-L_y, y-L_y, y, y+L_y]])

    for i in range(len(distance_array[0,:])):
        distance_array[0,i] =  distance(equivalent_a[:,i], b)

    dist = min(distance_array[0,:])
    return dist

#boolean ftion, arguments are 2 vectors a,b and radius R
#returns 1 if the two vectors are within distance R (smaller than R)
#works in periodical grid
def inRadius(a,b,R):
    dist = distance_periodic(a,b)
    if dist < R:
        return True
    else:
        return False

#checks if the position of the vector is within boundaries, if not, this ftion applies periodic boundaries
#returns vector within periodic bounaries
def periodic_check(vector, boundary_x, boundary_y):
    x = vector[0]
    y = vector[1]

    if x < 0:
        x = x+boundary_x
    elif x>boundary_x:
        x = x-boundary_x

    if y < 0:
        y = y+boundary_y
    elif y>boundary_x:
        y = y-boundary_y
    vector = np.array([x,y])
    return vector

#sums all the unit vectors in radius R in the position vector and return the unit vector in that direction
def direction_sum(direction_vec, position_vec, position ):
    sum = [0,0]
    for i in range(len(direction_vec[0,:])):
        if inRadius(position_vec[:,position], position_vec[:,i], R):
            sum = sum + direction_vec[:,i]

    norm = la.norm(sum)

    if norm == 0:       #if the sum of the vectors gives 0, return the currecnt direction (don't change anything)
        #print("old_dir")
        return np.array(direction_vec[:,position])
    #print (sum/norm)
    return sum/norm     #otherwise returns unit vector in the mean direction

#returns random angle in radians from a normal distribution with variance = sigma
def random_angle(sigma):
    angle = rd.normal(0,sigma)
    return angle

# computes next position and direction of a particle i
# returns two vectors - one with new position and 2nd with new direction
def step(i):

    #update position vector
    position = positions[:,i]
    direction = directions[:,i]
    new_position =  position + v*direction
    new_position = periodic_check(new_position,L_x, L_y)

    #update direction vector
    new_direction = direction_sum(directions, positions, i)         #gives me unit vector in the mean direction
    random = random_angle(sigma)                                    #random angle generation
    angle = np.arctan(new_direction[1]/new_direction[0])            #get the angle of the sums
    new_direction = [np.cos(angle + random),np.sin(angle + random)] #get the components of the new direction

    return new_position, new_direction

#returns the measure of order in the directions - for perfectly aligned system returns 1
def order_parameter(direction_vec):
    sum = [0,0]
    for i in range(N):
        sum = sum + direction_vec[:,i]
        result = np.sqrt(np.dot(sum,sum))
    return result/N

#arguments are the array you want to find mean for and visual_inspection - how many initial values you chop set_offsets
#returns mean of the sliced array
def mean_order_parameter(array, visual_inspection) :
    mean = -1 #initialisation
    cut_array = array[visual_inspection:] #slicing the array, excluding first n values, n=visual_inspection
    mean = np.mean(cut_array)

    return mean

#this functions completes the simulation and returns resulting mean order parameters
#argument is st. deviation (sigma)
def Simulate(s):
    global sigma
    sigma = s

    new_positions = np.zeros((2,N))
    new_directions = np.zeros((2,N))
    order_parameters = np.zeros((1,steps))

    #setting up the grid
    global positions,directions
    positions = generate_positions(N,L_x,L_y)   #generates positions matrix #                                                                        #
    directions = np.zeros((2, N))               #generates directions matrix#
    for i in range(N):
        directions[:,i] = random_unit_vector()

    for i in range(steps):

        print(i) ### testing of the position - I know which step the computer is at
        for k in range(N):
            new_pos, new_dir = step (k)
            new_positions[:,k] = new_pos[:]
            new_directions[:,k] = new_dir[:]
        order_parameters[0,i] = order_parameter(directions)
        #print (order_parameters) ### test

        positions = new_positions.copy()        # copy the new vector over the old one
        directions = new_directions.copy()      # (this copying is so that I can do the whole process
                                                # of updating angles in one go, without changing the
                                                # position and direction vector while going in the loop)
    step_vector = np.linspace(0,steps,steps)
    #print(mean_order_parameter(order_parameters[0,:], visual_inspection))
    mean = mean_order_parameter(order_parameters[0,:], visual_inspection)

    return mean


######################################################################
######################################################################


sigma_vector = np.arange(sigma_min, sigma_max,sigmaStep)
mean_vector = []
for i in range(len(sigma_vector)):
    mean_vector.append(Simulate(sigma_vector[i]))
print(mean_vector)

fit = np.polyfit(sigma_vector, mean_vector, 3) #fitting the data, polynomial of the order 3
polynomial = np.poly1d(fit)
polynomial_y = polynomial(sigma_vector)




plt.plot(sigma_vector ,mean_vector, marker = ".", label = "simulation data")

#uncomment bellow if you want a fit to be displayed:
#plt.plot(sigma_vector, polynomial_y, linestyle = "-", label = "polynomial fit")

plt.legend()
plt.title("Phase transition plot (variance vs. resulting mean order parameter)") #+ title)
plt.xlabel("Variance")
plt.ylabel("Order Parameter")

plt.show()
