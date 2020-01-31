#Copyright © 2019 ALEXANDER GROH
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

import numpy as np
from numpy import linalg as la
import numpy.random as rd
import random
import scipy
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as ani


#parameters
N = 100     # number of particles/bird
L_x = 10      # boundary in x
L_y = 10       # boundary in y
R = 1           # radius of birds' control
v = 0.5       # speed of each particle/bird
sigma_min = 0.1   # starting variance
sigma_max = 3     # ending variance
sigmaStep = 0.6   #step in variance
steps = 100   # steps to complete
#visual_inspection = 10  #after how many steps does order parameter converge
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
    if dist < R and dist != 0:
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
        return np.array([0,0])
    return sum/norm     #otherwise returns unit vector in the mean direction

#returns random angle in radians from a normal distribution with variance = sigma
def random_angle(sigma):
    angle = rd.normal(0,sigma)
    return angle

# computes next position and direction of a particle i
# returns two vectors - one with new position and 2nd with new direction
def step(i):#positions_vector, directions_vector, sigma):
    direction = directions[:,i]

    #update position
    new_position =  positions[:,i] + v*direction
    new_position = periodic_check(new_position,L_x, L_y)

    #update direction
    new_direction = direction_sum(directions, positions, i)      #gives me unit vector in the mean direction
    random = random_angle(sigma)
    old_angle = np.arctan(direction[1]/direction[0])
    new_direction = new_direction  + [np.cos(old_angle+ random),np.sin(old_angle+ random)]           #normalize the new vector
    norm = la.norm(new_direction)
                                        #normalize the added noise vector
    new_direction = new_direction/norm

    return new_position, new_direction

#returns the measure of order in the directions - for perfectly aligned system returns 1
def order_parameter(direction_vec):
    sum = [0,0]
    for i in range(N):
        sum = sum + direction_vec[:,i]
        result = np.sqrt(np.dot(sum,sum))
    return result/N


"""
#########################################################################
#########################################################################
#first SETUP                                                            #
positions = generate_positions(N,L_x,L_y)   #generates positions matrix #
                                                                        #
directions = np.zeros((2, N))               #generates directions matrix#
for i in range(N):                                                      #
    directions[:,i] = random_unit_vector()                              #
#########################################################################
"""
"""
#one full step
new_positions = np.zeros((2,N))
new_directions = np.zeros((2,N))
for i in range(N):
    new_pos, new_dir = step (i)
    new_positions[:,i] = new_pos
    new_directions[:,i] = new_dir
"""
#print(positions)
#print(new_positions)
#print(directions)
#print(new_directions)



#no animation
"""
#animation:
animation_positionX = np.zeros((steps, N))
animation_positionY = np.zeros((steps, N))
animation_directionX = np.zeros((steps, N))
animation_directionY = np.zeros((steps, N))
"""
sigma_vector = np.arange(sigma_min, sigma_max, sigmaStep)
for i in range(len(sigma_vector)):
    sigma = sigma_min + i *sigmaStep
    print(sigmaStep, sigma)

    #I have to zero the positions and set up the grid
        #n steps:
    new_positions = np.zeros((2,N))
    new_directions = np.zeros((2,N))
    order_parameters = np.zeros((1,steps))

    #setting up the grid
    positions = generate_positions(N,L_x,L_y)   #generates positions matrix #
                                                                            #
    directions = np.zeros((2, N))               #generates directions matrix#
    for i in range(N):                                                      #
        directions[:,i] = random_unit_vector()

    for i in range(steps):
        print(i) ### testing of the position - I know which step the computer is at
        for k in range(N):
            new_pos, new_dir = step (k)
            new_positions[:,k] = new_pos
            new_directions[:,k] = new_dir
        order_parameters[:,i] = order_parameter(directions)
        print (order_parameters) ### test

        positions = new_positions.copy()        # copy the new vector over the old one
        directions = new_directions.copy()      # (this copying is so that I can do the whole process
                                                # of updating angles in one go, without changing the
                                                # position and direction vector while going in the loop)
    step_vector = np.linspace(0,steps,steps)
    # graphs
    #title = "my title"
    #plt.plot(positions[0,:], positions[1,:], '.')
    plt.plot(step_vector ,order_parameters[0,:], label = ("sigma=" + str(sigma)),marker = ".")
    plt.legend()
    plt.title("Changes in Vicsek Order Parameter each step") #+ title)
    plt.xlabel("step")
    plt.ylabel("Vicsek Order Parameter")


plt.show()
