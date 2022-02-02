# Import packages
import numpy as np
import matplotlib.pyplot as plt
import math as m
import os
# Import User Defined Functions

# Initiate panel solver
geometry = 'Circle'
# geometry = 'fx76mp140_selig.txt'
num_panels = np.arange(4, 40, 2, dtype=int)
alpha_0 = 0


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
        edge[i] = (data[i + 1][0] - data[i][0] * (data[i + 1][1] - data[i][1]))  # Check value for each panel
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
    delta = phi + (np.pi / 2)
    beta = delta - (alpha * (np.pi/180))

    # Plot panels
    x = data[:, 0]
    y = data[:, 1]
    plt.figure(0)
    plt.fill(x, y, color='black', label='Panels')  # Plot geometry
    plt.title('Panel Geometry ' + str(num) + ' Panels')
    plt.ylabel('Y Axis')
    plt.ylim((y_bot_lim, y_top_lim))
    plt.xlabel('X Axis')
    plt.xlim((x_bot_lim, x_top_lim))
    if geom == 'Circle':
        # Add circle specific plot data
        plt.figure(0)
        circle = plt.Circle((0, 0), 1, fill=False, color='black', linestyle='--', label=' Unit Circle')  # Create circle
        plt.gca().add_patch(circle)  # Plot circle
    # Plot Boundary Points
    plt.scatter(x, y, label='Boundary Points', color='lime')
    # Plot Control Points
    plt.scatter(control_data[:, 0], control_data[:, 1], label='Control Points', color='cyan')
    plt.legend(loc='upper center', ncol=2)
    plt.tight_layout()
    # Save files in folder
    current_path = os.getcwd()
    plot_folder = current_path + '/' + 'Panels'
    if not os.path.exists(plot_folder):
        os.makedirs(plot_folder, exist_ok=True)
    plt.savefig(plot_folder + '/Circle' + str(num) + '.png', bbox_extra_artists='legend_outside,')
    # plt.show()
    plt.close()


for number in num_panels:
    source_panel(number, geometry, alpha_0)
