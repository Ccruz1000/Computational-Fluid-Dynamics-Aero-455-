import numpy as np
import matplotlib.pyplot as plt

dt = 0.1  # insert appropriate time step
Nsteps = 10  # insert appropriate total number of timesteps

## Setting up initial conditions (vortex centres and circulation)
# Vortex rings
# y_v = np.array([-2,2,-2,2])#insert the y-positions of the 4 vortices])
# x_v = np.array([1,1,3,3])#insert the x-positions of the 4 vortices])
# k_v = np.array([-1,1,-1,1])#insert the line vortex constant k of the 4 vortices])
y_v = [-1, 1, -1, 1]
x_v = [1, 1, 3, 3]
k_v = [-1, 1, -1, 1]
# Setting up the plot
plt.ion()
fig, ax = plt.subplots(1, 1)
# mark the initial positions of vortices
p, = ax.plot(x_v, y_v, '*k', markersize=10)
# play around with the marker size and type as you see fit

# draw the initial velocity streamline
ngrid = 2  # insert the dimension of your simulation grid
# Y, X = np.mgrid[-ngrid:ngrid:360j, -ngrid:ngrid:360j]
Y, X = np.mgrid[-2:2:360j, 0:44:360j]
# 360j sets the resolution of the cartesian grid; play around with it as you see fit
vel_x = np.zeros(np.shape(X))  # this holds x−velocity
vel_y = np.zeros(np.shape(Y))  # this holds y−velocity

# masking radius for better visualization of the vortex centres
r_mask = 1  # insert the radius of the mask around the vortex centres
# within this mask you will not plot any streamline
# so that you can see more clearly the movement of the vortex centres

for i in range(len(x_v)):  # looping over each vortex
    # insert lines for computing the total velocity field
    for j in range(len(vel_x[0])):  # looping through the field columns
        for k in range(len(vel_x)):  # looping through the field rows
            diff_y = y_v[i] - Y[k][j]
            diff_x = x_v[i] - X[k][j]
            # and adding the advection velocity for each point on the field
            # if Y[]
            if (y_v[i] - Y[k][j]) ** 2 + (x_v[i] - X[k][j]) ** 2 != 0:
                vel_x[k][j] += k_v[i] * diff_y / (diff_x ** 2 + diff_y ** 2)
                vel_y[k][j] += -k_v[i] * diff_x / (diff_x ** 2 + diff_y ** 2)

            # set up the boundaries of the simulation box
ax.set_xlim([0, 44])
ax.set_ylim([-ngrid, ngrid])

# initial plot of the streamlines
ax.streamplot(X, Y, vel_x, vel_y, density=[1, 1])
# play around with density as you see fit;
# see the API doc for more detail

fig.canvas.draw()

# Evolution
count = 0
vel_x_vortex = [0, 0, 0, 0]
vel_y_vortex = [0, 0, 0, 0]

while count < Nsteps:
    ##Compute and update advection velocity
    # insert lines to re-initialize the total velocity field
    for i in range(len(x_v)):
        # insert lines to compute the total advection velocity on each vortex
        # vel_x = np.zeros(np.shape(X))
        # vel_y = np.zeros(np.shape(Y))
        for j in range(len(x_v)):
            diff_y = -y_v[i] + y_v[j]
            diff_x = -x_v[i] + x_v[j]
            if diff_y ** 2 + diff_x ** 2 != 0:
                vel_x_vortex[i] += k_v[j] * diff_y / (diff_x ** 2 + diff_y ** 2)
                vel_y_vortex[i] += -k_v[j] * diff_x / (diff_x ** 2 + diff_y ** 2)
                # insert lines to update the positions of vortices
    for i in range(len(x_v)):
        x_v[i] = x_v[i] + vel_x_vortex[i] * dt
        y_v[i] = y_v[i] + vel_y_vortex[i] * dt

    print(x_v)
    print(y_v)
    # insert lines to re-initialize the total velocity field
    for i in range(len(vel_x[0])):
        for j in range(len(vel_x)):
            vel_x[j][i] = 0
            vel_y[j][i] = 0
    # vel_x = np.zeros(np.shape(X))
    # vel_y = np.zeros(np.shape(Y))
    for i in range(len(x_v)):
        # insert lines to update the streamlines
        for j in range(len(vel_x[0])):  # looping through the field columns
            for k in range(len(vel_x)):  # looping through the field rows
                vortex_x = x_v[i]
                vortex_y = y_v[i]
                diff_y = y_v[i] - Y[k][j]
                diff_x = x_v[i] - X[k][j]
                # and adding the advection velocity for each point on the field
                if (y_v[i] - Y[k][j]) ** 2 + (x_v[i] - X[k][j]) ** 2 != 0:
                    vel_x[k][j] += k_v[i] * diff_y / (diff_x ** 2 + diff_y ** 2)
                    vel_y[k][j] += -k_v[i] * diff_x / (diff_x ** 2 + diff_y ** 2)

    ##update plot
    # the following two lines clear out the previous streamlines
    ax.collections = np.zeros(len())
    ax.patches = []

    p.set_xdata(x_v)
    p.set_ydata(y_v)

    ax.streamplot(X, Y, vel_x, vel_y, density=[1, 1])

    fig.canvas.draw()
    plt.pause(0.001)  # play around with the delay time for better visualization

    count += 1