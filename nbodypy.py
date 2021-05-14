#Christopher Butler
#This program calculates trajectories of n bodies through 3d cartesian space.

import numpy
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import math

#DECLARATION OF CONSTANTS 
G = 1
#G = 6.674e-11 #for SI units
#G = 0.0002959122083 # square of the gaussian gravitational constant, for AU units

#open input data file
inputFile = open(r"C:\\Users\Christopher\Desktop\Computational Physicc\nbodyinput.txt", "r")

line1 = inputFile.readline()

dt = float(line1.split()[1])

line2 = inputFile.readline()

timesteps = int(line2.split()[1])

line3 =inputFile.readline() #skip this header

n = len(open(r"C:\\Users\Christopher\Desktop\Computational Physicc\nbodyinput.txt", "r").readlines()) - 3

#Initialize matrices to store position data:

#X-position for each body at each time, initialized to zero for now:
PosArrayX = [[0.0 for i in range(n)] for j in range(timesteps+1)]
#likewise for y and z:
PosArrayY = [[0.0 for i in range(n)] for j in range(timesteps+1)]
PosArrayZ = [[0.0 for i in range(n)] for j in range(timesteps+1)]

#initialize arrays for x, y and z velocities of each body at each time:
Velx = [[0.0 for i in range(n)] for j in range(timesteps+1)]
Vely = [[0.0 for i in range(n)] for j in range(timesteps+1)]
Velz = [[0.0 for i in range(n)] for j in range(timesteps+1)]

#masses of the bodies:
Mass = [0.0]*n


i = 0 #body iterator
for line in inputFile:
    currentLineValues = line.split()
    PosArrayX[0][i] = float(currentLineValues[0])
    PosArrayY[0][i] = float(currentLineValues[1])
    PosArrayZ[0][i] = float(currentLineValues[2])
    Velx[0][i] = float(currentLineValues[3])
    Vely[0][i] = float(currentLineValues[4])
    Velz[0][i] = float(currentLineValues[5])
    Mass[i] = float(currentLineValues[6])
    
    i+=1

#print(PosArrayX)
#print(PosArrayY)
#print(PosArrayZ)

#function to find the net force on body j from all other n bodies, at one given point in time;
def netForce(PosX, PosY, PosZ, Mass, numberOfBodies, jthBody, time):
    netForces = [0.0]*3 #index 0 is X-comp of force, index 1 is Y-comp and index 2 is Z-comp
    
    netForceX = 0.0
    netForceY = 0.0
    netForceZ = 0.0

    for i in range (numberOfBodies): #loop through each other body
        if (i != jthBody): #don't calculate attraction from itself!

            rjiX = PosX[time][i] - PosX[time][jthBody] # X distance between body j and i
            rjiY = PosY[time][i] - PosY[time][jthBody] # Y distance between body j and i
            rjiZ = PosZ[time][i] - PosZ[time][jthBody] # Z distance between body j and i
            
            #print('rjiX = ',rjiX)
            #print('rjiY = ',rjiY)
            #print('rjiZ = ',rjiZ)

            rji = math.sqrt(rjiX*rjiX + rjiY*rjiY + rjiZ*rjiZ + 1e-8) #distance between body j and i. The 1e-8 is added as a "fudge factor" to avoid "collision" issues.

            #print('rji = ',rji)

            #Now for the application of Newton's inverse square law. Done on a per-coordinate basis.
            netForceX = netForceX + G*Mass[jthBody]*Mass[i]*(rjiX / (rji*rji*rji))
            netForceY = netForceY + G*Mass[jthBody]*Mass[i]*(rjiY / (rji*rji*rji))
            netForceZ = netForceZ + G*Mass[jthBody]*Mass[i]*(rjiZ / (rji*rji*rji))

            #print('netForceX = ', netForceX)
            #print('netForceY = ', netForceY)
            #print('netForceZ = ', netForceZ)
    
    netForces[0] = netForceX
    netForces[1] = netForceY
    netForces[2] = netForceZ

    #print('netForces[0] = ',netForces[0])
    #print('netForces[1] = ',netForces[1])
    #print('netForces[2] = ',netForces[2])

    return netForces


for i in range(timesteps): #This loop traverses the time dimension, doing simulation calculations at each time step.

    for j in range(n): #loops through each body
        
        # **calculate the net force on the j-th body**

        netForceOnj = netForce(PosArrayX, PosArrayY, PosArrayZ, Mass, n, j, i)
        #print('netForceOnj = ', netForceOnj)
        #**calculate net acceleration of j-th body using the the net force of body j and mass of body j**
		#Using F = m*a, we know that acceleration of body j is the net force upon it, divided by its mass;
        netAccelOnjX = (netForceOnj[0]) / Mass[j]
        netAccelOnjY = (netForceOnj[1]) / Mass[j]
        netAccelOnjZ = (netForceOnj[2]) / Mass[j]

        #print('netAccelonjX = ',netAccelOnjX)
        #print('netAccelonjY = ',netAccelOnjY)
        #print('netAccelonjZ = ',netAccelOnjZ)

        # **calculate velocity of body j by integrating acceleration of body j with the timestep**

        Velx[i+1][j] = netAccelOnjX*dt + Velx[i][j]
        Vely[i+1][j] = netAccelOnjY*dt + Vely[i][j]
        Velz[i+1][j] = netAccelOnjZ*dt + Velz[i][j]

        #print('Velx[i+1][j] = ',Velx[i+1][j])
        #print('Vely[i+1][j] = ',Vely[i+1][j])
        #print('Velz[i+1][j] = ',Velz[i+1][j])

        #calculate position of body j by integrating velocity
        PosArrayX[i+1][j] = Velx[i+1][j]*dt + PosArrayX[i][j]
        PosArrayY[i+1][j] = Vely[i+1][j]*dt + PosArrayY[i][j]
        PosArrayZ[i+1][j] = Velz[i+1][j]*dt + PosArrayZ[i][j]


xMax = float('-inf')
yMax = float('-inf')
zMax = float('-inf')
xMin = float('inf')
yMin = float('inf')
zMin = float('inf')

for j in range(n):
    for i in range(timesteps):
        if(PosArrayX[i][j] > xMax):
            xMax = PosArrayX[i][j]
        if(PosArrayY[i][j] > yMax):
            yMax = PosArrayY[i][j]
        if(PosArrayZ[i][j] > zMax):
            zMax = PosArrayZ[i][j]
        if(PosArrayX[i][j] < xMin):
            xMin = PosArrayX[i][j]
        if(PosArrayY[i][j] < yMin):
            yMin = PosArrayY[i][j]
        if(PosArrayZ[i][j] < zMin):
            zMin = PosArrayZ[i][j]


fig = plt.figure()
ax = fig.add_subplot( projection='3d')

#for i in range(timesteps):
#    for j in range(n):
#        ax.scatter(PosArrayX[i][j], PosArrayY[i][j], PosArrayZ[i][j], marker='.')


def animate(i):
    ax.cla()
    ax.set_xlim(xMin-0.1, xMax+0.1)
    ax.set_ylim(yMin-0.1, yMax+0.1)
    ax.set_zlim(zMin-0.1, zMax+0.1)
    for j in range(n):
        ax.scatter(PosArrayX[i][j], PosArrayY[i][j], PosArrayZ[i][j], marker='o')

anim = animation.FuncAnimation(fig, animate, frames = timesteps, interval = 1, blit = False)

#ax.scatter(PosArrayX[0:timesteps][0], PosArrayY[0:timesteps][0], PosArrayZ[0:timesteps][0])
#ax.scatter(PosArrayX[0:timesteps][1], PosArrayY[0:timesteps][1], PosArrayZ[0:timesteps][1])
#anim.save('fig8.gif', fps = 60, dpi = 200, writer='imagemagick')
plt.show()
#print(PosArrayX)
#print(PosArrayY)
#print(PosArrayZ)


#figure 8 input file, assuming SI units
"""
timestep:	0.01
timesteps:	2000
PosX	PosY	PosZ	VelX	VelY	VelZ	Mass
0	0	0	0.93240737	0.86473146	0	1
-0.97000436	0.24308753	0	-0.466203685 -0.43236573	0	1
0.97000436   -0.24308753   0   -0.466203685    -0.43236573   0   1


"""