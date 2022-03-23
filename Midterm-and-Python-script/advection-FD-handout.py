#!/usr/bin/python
"""A Finite Difference code to solve the 1D transient advection equation,
            dq/dt + a dq/dx = 0
"""

import numpy as np
import matplotlib.pyplot as plt


def step(t, t_end, step_value):
    """Returns a step function with maximum equal to step_value
    up to time t_end.
    """
    if t < t_end:
        value = step_value
    else:
        value = 0

    return value


def d1_o1_b2(dx, q, N):
    """Function to calculate the first derivative of a field q in numpy
    array form. Uses 1st order backward 2-point stencil.
    """
    N = q.shape[0]
    dqdx = np.zeros_like(q)

    # use 1st order forward differencing for inlet grid point
    dqdx[0] = (q[1] - q[0])/dx

    # loop over all other points with 1st order backward differencing
    for i in np.arange(1, N-1):
        dqdx[i] = (q[i] - q[i-1])/dx
    return dqdx


def d1_o2_c2(dx, q, N):
    """Function to calculate the first derivative of a field q in numpy
    array form. Uses 2nd order central 2-point stencil.
    """
    return None


def d1_o2_b3(dx, q, N):
    """Function to calculate the first derivative of a field q in numpy
    array form. Uses 2nd order backward 3-point stencil.
    """

    return None


# physical parameters
L = 1.0  # length of line domain
a = 1.0  # advection velocity a
# beta = 1e-4  # diffusion coefficient beta

# time discretization variables
dt = 1e-4  # time step
t_final = 0.75  # final time
t = 0.  # time variable
# t = np.zeros(int(t_final/dt))

# space discretization variables
N = 200  # number of 1D cells
dx = L/N  # grid spacing

# blending factor
# p = 0.25

# print("Number time steps: ", t.shape[0])
print("Grid spacing: ", dx)
print("Time step: ", dt)
# print("Blending factor: ", p)

# create array of x-coordinates x[N+1] and solution variables q[N+1]
x = np.zeros(N + 1)

# create solution arrays initialized to 0
q1 = np.zeros(N + 1)
q1_old = np.zeros(N + 1)
dq1dx = np.zeros(N + 1)

q2 = np.zeros(N + 1)
q2_old = np.zeros(N + 1)
dq2dx = np.zeros(N + 1)  # this is the central scheme solution

q3 = np.zeros(N + 1)
q3_old = np.zeros(N + 1)
dq3dx = np.zeros(N + 1)  # this is the 2nd order backward scheme solution

# initialize x[]
for i in np.arange(N):
    x[i + 1] += x[i] + L/N

# time loop
while t < t_final:
    # increment time
    t += dt

    q1[0] = step(t, 0.2, 1)

    # calculate first derivatives of previous time
    dq1dx = d1_o1_b2(dx, q1_old, N)

    # update solution at current time
    q1[1:-1] = q1_old[1:-1] - a*dt*dq1dx[1:-1]  # stops at -2
    # q1[-1] = q1[-2]  # zero gradient outlet BC
    q1_old = q1

plt.plot(x, q1, label='1st order backward')
# plt.plot(x, q2, label='2nd order central')
# plt.plot(x, q3, label='2nd order backward')
plt.legend()
plt.xlim(-.1, 1.1)
plt.ylim(-.5, 1.5)
plt.title('dt = ' + str(dt) + ', N = ' + str(N))
plt.grid()
plt.show()





