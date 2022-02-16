# Import packages
import numpy as np
import matplotlib.pyplot as plt
import math as m
import os
# Import User Defined Functions
from Geometric_Integral import geometric_integral

# Initiate panel solver
geometry1 = 'fx76mp140_selig'
geometry = 'Circle'
num_panels = np.arange(4, 20, 2, dtype=int)
alpha_0 = 0
V_inf = 1
number = 2


def source_panel(num, geom, alpha):
    if geom == 'Circle':
        x_top_lim = 1.5
        x_bot_lim = -1.5
        y_top_lim = 1.5
        y_bot_lim = -1.5
        delta_theta = 360 / (2 * num)  # Angle offset to calculate boundary points

        # Calculate angles for boundary points
        theta = np.linspace(0, 360, num + 1)  # Angle for boundary points [deg]
        theta = theta + delta_theta  # Angle offset for boundary points
        theta = theta * (np.pi / 180)  # Convert angle to rad

        # Calculate boundary points
        data = np.zeros((num + 1, 2))
        for i in range(num + 1):
            data[i][0] = (np.cos(theta[i]))
            data[i][1] = (np.sin(theta[i]))

    else:
        x_top_lim = 1.1
        x_bot_lim = -0.1
        y_top_lim = 0.5
        y_bot_lim = -0.5
        data = np.loadtxt(geom, skiprows=1, dtype='float')


    panel_num = len(data) - 1
    # Check direction of panels
    edge = np.zeros(panel_num)  # Initialize edge array
    for i in range(panel_num):
        edge[i] = (data[i + 1][0] - data[i][0]) * (data[i + 1][1] - data[i][1])  # Check value for each panel
    edge_sum = np.sum(edge)

    # Flip if panels are the wrong direction
    if edge_sum < 0:
        data = np.flipud(data)

    # Calculate control points
    control_data = np.zeros((panel_num, 2))
    panel_length = np.zeros(panel_num)
    phi = np.zeros(panel_num)  # Initialize panel orientation angle phi

    # Calculate geometric quantities for panels
    for i in range(panel_num):
        control_data[i][0] = 0.5 * (data[i][0] + data[i + 1][0])  # Calculate x control points
        control_data[i][1] = 0.5 * (data[i][1] + data[i + 1][1])  # Calculate y control points
        delta_x = data[i + 1][0] - data[i][0]  # Panel x length
        delta_y = data[i + 1][1] - data[i][1]  # Panel y length
        panel_length[i] = (delta_y ** 2 + delta_x ** 2) ** 0.5  # Calculate panel length
        phi[i] = m.atan2(delta_y, delta_x)  # Calculate panel orientation angle in rad
        if phi[i] < 0:
            phi[i] = phi[i] + 2 * np.pi  # Add 2pi rad to convert negative angle to positive
    # Compute panel angle w.r.t horizontal and w.r.t alpha (rad)
    delta = phi + (np.pi / 2)   # Panel normal angle (rad)
    beta = delta - (alpha * (np.pi/180))    # Angle between alpha and panel normal (rad)
    beta[beta > 2*np.pi] = beta[beta > 2 * np.pi] - 2*np.pi # Convert angles over 2pi to between 0 and 2pi

    # Calculate I and J integral
    I, J = geometric_integral(control_data, data, phi, panel_length)

    # Create A matrix
    A = np.zeros([panel_num, panel_num])
    for i in range(panel_num):
        for j in range(panel_num):
            if i == j:
                A[i, j] = np.pi
            else:
                A[i, j] = I[i, j]
    # Create b array in Ax + b = 0
    b = np.zeros(panel_num)
    for i in range(panel_num):
        b[i] = - V_inf * 2 * np.pi * np.cos(beta[i])

    # Compute source panel strength
    lam = np.linalg.solve(A, b)
    # Check source panel strength (Should be 0 assuming closed shape)
    print("Sum of L: ", sum(lam*panel_length))

    # Calculate velocities and pressure coefficient
    vt = np.zeros(panel_num)
    cp = np.zeros(panel_num)

    for i in range(panel_num):
        cntr = 0
        for j in range(panel_num):
            cntr = cntr + (lam[j] / (2 * np.pi)) * J[i, j]

        vt[i] = V_inf * np.sin(beta[i]) + cntr
        cp[i] = 1 - (vt[i] / V_inf)**2



    # Plot panels
    x = data[:, 0]
    y = data[:, 1]
    plt.figure(0)
    plt.fill(x, y, color='black', label='Panels')  # Plot geometry
    plt.ylabel('Y Axis')
    plt.ylim((y_bot_lim, y_top_lim))
    plt.xlabel('X Axis')
    plt.xlim((x_bot_lim, x_top_lim))
    if geom == 'Circle':
        # Add circle specific plot data
        plt.figure(0)
        circle = plt.Circle((0, 0), 1, fill=False, color='black', linestyle='--', label=' Unit Circle')  # Create circle
        plt.gca().add_patch(circle)  # Plot circle
        plt.title('Panel Geometry ' + str(num) + ' Panels')
    else:
        plt.title(geometry1)
    # Plot Boundary Points
    plt.scatter(x, y, label='Boundary Points', color='r')
    # Plot Control Points
    plt.scatter(control_data[:, 0], control_data[:, 1], label='Control Points', color='b')
    plt.legend(loc='upper center', ncol=2)
    plt.tight_layout()
    # Save files in folder
    current_path = os.getcwd()
    plot_folder = current_path + '/' + 'Panels'
    if not os.path.exists(plot_folder):
        os.makedirs(plot_folder, exist_ok=True)
    plt.savefig(plot_folder + '/' + geometry + str(num) + '.png', bbox_extra_artists='legend_outside')
    # plt.show()
    plt.close()
    plt.figure(1)
    plt.scatter(beta * (180/np.pi), cp)
    plt.show()


# for number in num_panels:
#     source_panel(number, geometry, alpha_0)
source_panel(400, geometry1 + '.txt', alpha_0)
