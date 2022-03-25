#!/usr/bin/python
"""A Finite Difference code to solve the 1D transient advection equation,
            dq/dt + a dq/dx = 0
"""
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.pyplot import figure

figure(figsize=(10, 8))


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


def bonus_d1_o1_b2(dx, q, N):
    """Function to calculate the first derivative of a field q in numpy
    array form. Uses 1st order backward 2-point stencil.
    """
    N = q.shape[0]
    dqdx = np.zeros_like(q)
    dq2dx = np.zeros_like(q)

    # use 1st order forward differencing for inlet grid point
    dqdx[0] = (q[1] - q[0])/dx

    # loop over all other points with 1st order backward differencing
    for i in np.arange(1, N-1):
        dqdx[i] = (q[i] - q[i-1])/dx

    # Use forward differencing to find inlet grid point for second derivative
    dq2dx[0] = (q[2] + q[0] - 2 * q[1]) / (2 * dx ** 2)

    # Use backward differencing to find inlet grid point for second derivative
    dq2dx[-1] = (q[-1] - 2 * q[-2] + q[-3]) / (dx ** 2)

    # Use 2nd order central difference with 3 point stencil for all other points in second derivative
    for i in np.arange(1, N-1):
        dq2dx[i] = (q[i-1] - 2 * q[i] + q[i+1]) / (2 * dx ** 2)
    return dqdx, dq2dx


def bonus_d1_o2_c2(dx, q, N):
    """Function to calculate the first derivative of a field q in numpy
    array form. Uses 2nd order central 2-point stencil.
    """
    N = q.shape[0]
    dqdx = np.zeros_like(q)
    dq2dx = np.zeros_like(q)

    # Use 2nd order forward differencing for inlet grid point for first derivative
    dqdx[0] = (q[1] - q[0]) / (2 * dx)

    # Use 2nd order central differencing for all other points for first derivative
    for i in np.arange(1, N-1):
        dqdx[i] = (q[i+1] - q[i-1]) / (2 * dx)

    # Use 2nd order backward differencing for outlet grid point for first derivative
    dqdx[-1] = (q[-2] - q[-1]) / (2 * dx)

    # Use forward differencing to find inlet grid point for second derivative
    dq2dx[0] = (q[2] + q[0] - 2 * q[1]) / (dx ** 2)

    # Use 2nd order central difference with 3 point stencil for all other points in second derivative
    for i in np.arange(1, N - 1):
        dq2dx[i] = (q[i - 1] - 2 * q[i] + q[i + 1]) / (dx ** 2)

    # Use backward differencing to find outlet grid point for second derivative
    dq2dx[-1] = (q[-1] - 2 * q[-2] + q[-3]) / (dx ** 2)
    return dqdx, dq2dx


def bonus_d1_o2_b3(dx, q, N):
    """Function to calculate the first derivative of a field q in numpy
    array form. Uses 2nd order backward 3-point stencil.
    """
    N = q.shape[0]
    dqdx = np.zeros_like(q)
    dq2dx = np.zeros_like(q)

    # Use 2nd order forward differencing with 3 point stencil for first 2 inlet grid points
    dqdx[0] = (-3*q[0] + 4*q[1] - q[2]) / (2 * dx)
    dqdx[1] = (-3*q[1] + 4*q[2] - q[3]) / (2 * dx)

    # Use 2nd order backwards scheme with 3 point stencil for all other points
    for i in np.arange(2, N-1):
        dqdx[i] = (3 * q[i] - 4 * q[i-1] + q[i-2]) / (2 * dx)

        # Use forward differencing to find inlet grid point for second derivative
    dq2dx[0] = (q[2] + q[0] - 2 * q[1]) / (2 * dx ** 2)

    # Use backward differencing to find inlet grid point for second derivative
    dq2dx[-1] = (q[-1] - 2 * q[-2] + q[-3]) / (dx ** 2)

    # Use 2nd order central difference with 3 point stencil for all other points in second derivative
    for i in np.arange(1, N - 1):
        dq2dx[i] = (q[i - 1] - 2 * q[i] + q[i + 1]) / (2 * dx ** 2)
    return dqdx, dq2dx


def linear_advection_diffusion(N, L = 1, a = 1, beta = 1e-3, dt=1e-4):
    # physical parameters
    L = 1.0  # length of line domain
    a = 1.0  # advection velocity a

    # time discretization variables
    # dt = (courant * L/N) / a  # time step
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
    dq1dx = np.zeros(N + 1)  # This is the 1st order backward scheme solution
    dq1dx2 = np.zeros(N + 1)  # This is the 2nd order central scheme for second derivative

    q2 = np.zeros(N + 1)
    q2_old = np.zeros(N + 1)
    dq2dx = np.zeros(N + 1)  # this is the central scheme solution
    dq2dx2 = np.zeros(N + 1)  # This is the 2nd order central scheme for second derivative

    q3 = np.zeros(N + 1)
    q3_old = np.zeros(N + 1)
    dq3dx = np.zeros(N + 1)  # this is the 2nd order backward scheme solution
    dq3dx2 = np.zeros(N + 1)  # This is the 2nd order central scheme for second derivative

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
        dq1dx, dq1dx2 = bonus_d1_o1_b2(dx, q1_old, N)
        dq2dx, dq2dx2 = bonus_d1_o2_c2(dx, q2_old, N)
        dq3dx, dq3dx2 = bonus_d1_o2_b3(dx, q3_old, N)
        # update solution at current time
        q1[1:-1] = q1_old[1:-1] - a * dt * dq1dx[1:-1] + beta * dt * dq1dx2[1:-1]  # stops at -2
        q2[1:-1] = q2_old[1:-1] - a * dt * dq2dx[1:-1] + beta * dt * dq2dx2[1:-1]
        q3[1:-1] = q3_old[1:-1] - a * dt * dq3dx[1:-1] + beta * dt * dq3dx2[1:-1]
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


def d_backward_first_order(N, L = 1, a = 1, beta = 1e-4, dt=1e-4):
    # time discretization variables
    # dt = (courant * L/N) / a  # time step
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
    print("Beta: ", beta)
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
        dq1dx, dq2dx = bonus_d1_o1_b2(dx, q1_old, N)
        # update solution at current time
        q1[1:-1] = q1_old[1:-1] - a*dt*dq1dx[1:-1] + beta * dt * dq2dx[1:-1]  # stops at -2
        q1[-1] = q1[-2]  # zero gradient outlet BC
        q1_old = q1

    # plt.plot(x, q1, label='1st order backward')
    # plt.legend()
    # plt.title('dt = ' + str(dt) + ', N = ' + str(N) + ', dx = ' + str(dx))
    # plt.grid()
    # plt.show()
    return q1, x, dt, dx


def d_central_second_order(N, L = 1, a = 1, beta = 1e-4, dt = 1e-4):
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
        dq1dx, dq2dx = bonus_d1_o2_c2(dx, q1_old, N)
        # update solution at current time
        q1[1:-1] = q1_old[1:-1] - a * dt * dq1dx[1:-1] + beta * dt * dq2dx[1:-1]  # stops at -2
        q1[-1] = q1[-2]  # zero gradient outlet BC
        q1_old = q1

    # plt.plot(x, q1, label='1st order backward')
    # plt.legend()
    # plt.title('dt = ' + str(dt) + ', N = ' + str(N) + ', dx = ' + str(dx))
    # plt.grid()
    # plt.show()
    return q1, x, dt, dx


def d_backward_second_order(N, L = 1, a = 1, beta = 1e-3, dt = 1e-4):
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
    dq1dx = np.zeros(N + 1)  # This is the first derivative
    dq2dx = np.zeros(N + 1)  # This is the second derivative

    # initialize x[]
    for i in np.arange(N):
        x[i + 1] += x[i] + L/N

    # time loop
    while t < t_final:
        # increment time
        t += dt

        q1[0] = step(t, 0.2, 1)

        # calculate first derivatives of previous time
        dq1dx, dq2dx = bonus_d1_o2_b3(dx, q1_old, N)
        # update solution at current time
        q1[1:-1] = q1_old[1:-1] - a * dt * dq1dx[1:-1] + beta * dt * dq2dx[1:-1]  # stops at -2
        q1[-1] = q1[-2]  # zero gradient outlet BC
        q1_old = q1

    # plt.plot(x, q1, label='1st order backward')
    # plt.legend()
    # plt.title('dt = ' + str(dt) + ', N = ' + str(N) + ', dx = ' + str(dx))
    # plt.grid()
    # plt.show()
    return q1, x, dt, dx
