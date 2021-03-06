#Copyright © 2019 ALEXANDER GROH
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#this is little bit more fancy animation with current order parameter displayed
#for more detailed comments of code, refer to basic_animation.py
import numpy as np
from numpy import linalg as la
import numpy.random as rd
import random
import scipy
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as ani


#parameters
N = 30       # number of birds
L_x = 5      # boundary in x
L_y = 5      # boundary in y
v = 0.3      # speed of each particle/bird
sigma = 0.1  # st. deviation of random angle
steps = 100   # steps to complete
R = 0.4      # radius of birds' control

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
    if dist < R :
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
        return np.array(direction_vec[:,position])
    else:
        return sum/norm     #otherwise returns unit vector in the mean direction

#returns random angle in radians from a normal distribution with variance = sigma
def random_angle(sigma):
    angle = rd.normal(0,sigma)
    return angle

# computes next position and direction of a particle i
# returns two vectors - one with new position and 2nd with new direction
def step(i):#positions_vector, directions_vector, sigma):
    #update position vector
    position = positions[:,i]
    direction = directions[:,i]

    #update position
    new_position =  position + v*direction
    new_position = periodic_check(new_position,L_x, L_y)

    #update direction
    new_direction = direction_sum(directions, positions, i)      #gives me unit vector in the mean direction
    random = random_angle(sigma)
    angle = np.arctan(new_direction[1]/new_direction[0])

    new_direction = [np.cos(angle + random),np.sin(angle + random)]
    # no need to normalize the new vector because (cos(x)^2 + sin(x)^2 = 1)

    return new_position, new_direction

#returns the measure of order in the directions - for perfectly aligned system returns 1
def order_parameter(direction_vec):
    sum = [0,0]
    for i in range(N):
        sum = sum + direction_vec[:,i]
        result = np.sqrt(np.dot(sum,sum))
    return result/N



#########################################################################
#########################################################################
#initial SETUP                                                            #
positions = generate_positions(N,L_x,L_y)   #generates positions matrix #
                                                                        #
directions = np.zeros((2, N))               #generates directions matrix#
for i in range(N):                                                      #
    directions[:,i] = random_unit_vector()                              #
#########################################################################

#n steps:
new_positions = np.zeros((2,N))
new_directions = np.zeros((2,N))
order_parameters = np.zeros((1,steps))

#animation:
animation_positionX = np.zeros((steps, N))
animation_positionY = np.zeros((steps, N))
animation_directionX = np.zeros((steps, N))
animation_directionY = np.zeros((steps, N))

for i in range(steps):
    print(i) ### testing of the position - I know which step the computer is at
    for k in range(N):
        new_pos, new_dir = step (k)
        new_positions[:,k] = new_pos[:]
        new_directions[:,k] = new_dir[:]
        #print(new_directions[:,k])

    order_parameters[:,i] = order_parameter(directions)
    positions = new_positions.copy()        # copy the new vector over the old one
    directions = new_directions.copy()      # (this copying is so that I can do the whole process
                                            # of updating angles in one go, without changing the
                                            # position and direction vector while going in the loop)

    #matrices for the animation
    animation_positionX[i] = positions[0]
    animation_positionY[i] = positions[1]
    animation_directionX[i] = directions[0]
    animation_directionY[i] = directions[1]

#print(directions)
#print(order_parameters)

#animation code
fig, ax = plt.subplots()
plt.title("Sigma = " + str(sigma) + ", R = " + str(R) )
time_text = ax.text(0.45,0.95,"",transform = ax.transAxes, ha="right") #adding order parameter
t = 0

x = []
y = []

X = animation_positionX[0]
Y = animation_positionY[0]
U = 0
V = 0
Q = ax.quiver(X, Y, U, V, pivot='mid',units='inches', color='black', headwidth=20, headlength=10, linewidth=0,headaxislength=10, minlength = 0.1 , scale = 9)

mat, = ax.plot(x, y, '.', color = 'black')


def animate(i):
    x = animation_positionX [i,:]
    y = animation_positionY[i,:]

    U = animation_directionX[i, :]
    V = animation_directionY[i, :]

    Q.set_offsets(np.transpose(np.array([x,y])))
    Q.set_UVC(U,V)
    time_text.set_text("Vicsek order parameter = " + str(round(order_parameters[0,i+1],2)))

    mat.set_data(x, y)
    return mat, time_text

ax.axis([0,L_x,0,L_y])
ani = ani.FuncAnimation(fig, animate, interval=80, frames=np.arange(0, steps-1, 1))
ani.save("final_output.gif" , writer="imagemagick")
plt.show()
