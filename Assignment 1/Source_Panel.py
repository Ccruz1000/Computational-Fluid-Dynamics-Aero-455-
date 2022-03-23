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
# num_panels = np.arange(4, 20, 2, dtype=int)
# num_panels = [6, 12, 24, 48, 96, 192, 384, 768]
num_panels = [4,6, 8, 10]
# for i in range(2, 14):
#     num_panels.append(2**i)

alpha_0 = 0
V_inf = 1
number = 2


def source_panel(num, geom, alpha):
    print(num)
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
        edge[i] = (data[i + 1][0] - data[i][0]) * (data[i + 1][1] + data[i][1])  # Check value for each panel
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
    print(str(num) + ' Panels')
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

    # Calculate analytical value for pressure coefficient and velocity
    angle = np.linspace(0, 360, 360)
    cp_analytical = 2*np.cos(2*angle*np.pi/180)-1
    u_theta = V_inf * 2 * np.sin(angle * np.pi / 180)

    residual = sum(lam*panel_length)

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
    plt.scatter(x, y, label='Boundary Points', color='r', s=8)
    # Plot Control Points
    plt.scatter(control_data[:, 0], control_data[:, 1], label='Control Points', color='b', s=8)
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
    # Plot Pressure Coefficient
    plt.figure(1)
    # Plot calculated Cp
    plt.scatter(beta * (180/np.pi), cp, s=12, label='Source Panel Method', color='b')
    # Plot analytical Cp
    plt.plot(angle, cp_analytical, label='Analytical Solution', color='r')
    plt.legend(loc='best')
    plt.xlabel('Theta (Deg)')
    plt.ylabel('Cp')
    plt.title('Pressure Coefficient Comparison ' + str(num) + ' Panels')
    # Save files in folder
    plot_folder = current_path + '/' + 'Pressure Coefficient'
    if not os.path.exists(plot_folder):
        os.makedirs(plot_folder, exist_ok=True)
    plt.savefig(plot_folder + '/' + geometry + str(num) + '.png', bbox_extra_artists='legend_outside')
    # plt.show()
    plt.close()
    # Plot velocity
    plt.figure(2)
    # Plot calculated velocity
    plt.scatter(beta * (180/np.pi), vt, s=12, label='Source Panel Method', color='b')
    # Plot analytical velocity
    plt.plot(angle, u_theta, label='Analytical Solution', color='r')
    plt.legend(loc='best')
    plt.xlabel('Theta (Deg)')
    plt.ylabel('Velocity (m/s)')
    plt.title('Velocity Comparison ' + str(num) + ' Panels')
    plot_folder = current_path + '/' + 'Velocity'
    if not os.path.exists(plot_folder):
        os.makedirs(plot_folder, exist_ok=True)
    plt.savefig(plot_folder + '/' + geometry + str(num) + '.png', bbox_extra_artists='legend_outside')
    # plt.show()
    plt.close()
    return residual


residual1 = []
for number in num_panels:
    check = source_panel(number, geometry, alpha_0)
    residual1.append(check)
# source_panel(76, geometry1 + '.txt', 0)

plt.figure(3)
plt.plot(num_panels, residual1, color='r', label='Residual')
plt.title('Source Panel Method Residual ')
plt.xlabel('Number of Panels')
plt.ylabel('Residual')
plt.legend(loc='best')
current_path = os.getcwd()
plot_folder = current_path + '/' + 'Residual'
if not os.path.exists(plot_folder):
    os.makedirs(plot_folder, exist_ok=True)
plt.savefig(plot_folder + '/Residual Question 1.png', bbox_extra_artists='legend_outside')
plt.show()
plt.close()
source_panel(16138, geometry, alpha_0)
