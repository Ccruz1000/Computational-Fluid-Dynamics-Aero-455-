#!/usr/bin/python
"""A Finite Difference code to solve the 1D transient advection equation,
            dq/dt + a dq/dx = 0
"""

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.pyplot import figure
figure(figsize=(8, 6))
num = [100, 150, 200, 250, 400]
courant = 0.05
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
    N = q.shape[0]
    dqdx = np.zeros_like(q)

    # Use 2nd order forward differencing for inlet grid point
    dqdx[0] = (q[1] - q[0]) / (2 * dx)

    # Use 2nd order central differencing for all other points
    for i in np.arange(1, N-1):
        dqdx[i] = (q[i+1] - q[i-1]) / (2 * dx)
    # Use 2nd order backward differencing for outlet grid point
    dqdx[-1] = (q[-2] - q[-1]) / (2 * dx)
    return dqdx


def d1_o2_b3(dx, q, N):
    """Function to calculate the first derivative of a field q in numpy
    array form. Uses 2nd order backward 3-point stencil.
    """
    N = q.shape[0]
    dqdx = np.zeros_like(q)

    # Use 2nd order forward differencing with 3 point stencil for first 2 inlet grid points
    dqdx[0] = (-3*q[0] + 4*q[1] - q[2]) / (2 * dx)
    dqdx[1] = (-3*q[1] + 4*q[2] - q[3]) / (2 * dx)
    # Use 2nd order backwards scheme with 3 point stencil for all other points
    for i in np.arange(2, N-1):
        dqdx[i] = (3 * q[i] - 4 * q[i-1] + q[i-2]) / (2 * dx)
    return dqdx

def linear_advection (N):
    # physical parameters
    L = 1.0  # length of line domain
    a = 1.0  # advection velocity a
    # beta = 1e-4  # diffusion coefficient beta

    # time discretization variables
    dt = (courant * L/N) / a  # time step
    # dt = 1e-4
    t_final = 0.75  # final time
    t = 0.  # time variable
    # t = np.zeros(int(t_final/dt))

    # space discretization variables
    # N = 200  # number of 1D cells
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
    dq1dx = np.zeros(N + 1) # This is the 1st order backward scheme solution

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
        q2[0] = step(t, 0.2, 1)
        q3[0] = step(t, 0.2, 1)

        # calculate first derivatives of previous time
        dq1dx = d1_o1_b2(dx, q1_old, N)
        dq2dx = d1_o2_c2(dx, q2_old, N)
        dq3dx = d1_o2_b3(dx, q3_old, N)
        # update solution at current time
        q1[1:-1] = q1_old[1:-1] - a*dt*dq1dx[1:-1]  # stops at -2
        q2[1:-1] = q2_old[1:-1] - a*dt*dq2dx[1:-1]
        q3[1:-1] = q3_old[1:-1] - a*dt*dq3dx[1:-1]
        # q1[-1] = q1[-2]  # zero gradient outlet BC
        q1_old = q1
        q2_old = q2
        q3_old = q3

    plt.plot(x, q1, label='1st order backward')
    plt.plot(x, q2, label='2nd order central')
    plt.plot(x, q3, label='2nd order backward')
    plt.legend()
    plt.title('dt = ' + str(dt) + ', N = ' + str(N) + ', dx = ' + str(dx))
    plt.grid()
    plt.show()
    return None


def first_order_backward (N):
    # physical parameters
    L = 1.0  # length of line domain
    a = 1.0  # advection velocity a
    # beta = 1e-4  # diffusion coefficient beta

    # time discretization variables
    dt = (courant * L/N) / a  # time step
    # dt = 1e-4
    t_final = 0.75  # final time
    t = 0.  # time variable
    # t = np.zeros(int(t_final/dt))

    # space discretization variables
    # N = 200  # number of 1D cells
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
    dq1dx = np.zeros(N + 1) # This is the 1st order backward scheme solution


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

    # plt.plot(x, q1, label='1st order backward')
    # plt.legend()
    # plt.title('dt = ' + str(dt) + ', N = ' + str(N) + ', dx = ' + str(dx))
    # plt.grid()
    # plt.show()
    return q1, x, dt, dx


def central_second_order (N):
    # physical parameters
    L = 1.0  # length of line domain
    a = 1.0  # advection velocity a
    # beta = 1e-4  # diffusion coefficient beta

    # time discretization variables
    dt = (courant * L/N) / a  # time step
    # dt = 1e-4
    t_final = 0.75  # final time
    t = 0.  # time variable
    # t = np.zeros(int(t_final/dt))

    # space discretization variables
    # N = 200  # number of 1D cells
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
    dq1dx = np.zeros(N + 1) # This is the 1st order backward scheme solution


    # initialize x[]
    for i in np.arange(N):
        x[i + 1] += x[i] + L/N

    # time loop
    while t < t_final:
        # increment time
        t += dt

        q1[0] = step(t, 0.2, 1)

        # calculate first derivatives of previous time
        dq1dx = d1_o2_c2(dx, q1_old, N)
        # update solution at current time
        q1[1:-1] = q1_old[1:-1] - a*dt*dq1dx[1:-1]  # stops at -2
        # q1[-1] = q1[-2]  # zero gradient outlet BC
        q1_old = q1

    # plt.plot(x, q1, label='1st order backward')
    # plt.legend()
    # plt.title('dt = ' + str(dt) + ', N = ' + str(N) + ', dx = ' + str(dx))
    # plt.grid()
    # plt.show()
    return q1, x, dt, dx


def backward_second_order (N):
    # physical parameters
    L = 1.0  # length of line domain
    a = 1.0  # advection velocity a
    # beta = 1e-4  # diffusion coefficient beta

    # time discretization variables
    dt = (courant * L/N) / a  # time step
    # dt = 1e-4
    t_final = 0.75  # final time
    t = 0.  # time variable
    # t = np.zeros(int(t_final/dt))

    # space discretization variables
    # N = 200  # number of 1D cells
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
    dq1dx = np.zeros(N + 1) # This is the 1st order backward scheme solution


    # initialize x[]
    for i in np.arange(N):
        x[i + 1] += x[i] + L/N

    # time loop
    while t < t_final:
        # increment time
        t += dt

        q1[0] = step(t, 0.2, 1)

        # calculate first derivatives of previous time
        dq1dx = d1_o2_b3(dx, q1_old, N)
        # update solution at current time
        q1[1:-1] = q1_old[1:-1] - a*dt*dq1dx[1:-1]  # stops at -2
        # q1[-1] = q1[-2]  # zero gradient outlet BC
        q1_old = q1

    # plt.plot(x, q1, label='1st order backward')
    # plt.legend()
    # plt.title('dt = ' + str(dt) + ', N = ' + str(N) + ', dx = ' + str(dx))
    # plt.grid()
    # plt.show()
    return q1, x, dt, dx
dt = []
q = []
dx = []
x = []
for i in range(len(num)):
    q_sol, x_sol, dt_sol, dx_sol = backward_second_order(num[i])
    q.append(q_sol)
    x.append(x_sol)
    dt.append(dt_sol)
    dx.append(dx_sol)

plt.figure(1)
for i in range(len(num)):
    plt.plot(x[i], q[i], label='2nd Order Backward, ' + str(num[i]) + ' Divisions')
plt.title('2nd Order Central Backward')
plt.grid()
plt.legend(loc='best')
plt.show()


