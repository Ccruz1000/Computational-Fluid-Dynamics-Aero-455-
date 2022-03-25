from linear_advection import *
from linear__advection_diffusion import *
import matplotlib.pyplot as plt
# Define values
L = 1.  # Length
a = 1.  # Advection Value
beta = 1e-4  # Diffusion Coefficient
num = [50, 100, 200, 300, 400, 500]  # Number of divisions
courant = 0.25 # Courant Number

equation = 'Linear Advection'  # Equation
scheme = '2nd Order Backward, '  # Scheme for titles

# Initialize empty variables for plotting
dt = []
q = []
dx = []
x = []

'''
Input which scheme and type you'd like to use here. Place a_ in front for linear advection, or d_ in front for linear 
advection diffusion. Can also call linear_advection(), or linear_advection_diffusion to see all 3 schemes. The 
available schemes are shown. 
1st order backwards - first_order_backward
2nd order backwards - second_order_backward
2nd order central - second_order_central
'''
for i in range(len(num)):
    q_sol, x_sol, dt_sol, dx_sol = a_backward_second_order(num[i], L, a, beta)
    q.append(q_sol)
    x.append(x_sol)
    dt.append(dt_sol)
    dx.append(dx_sol)
plt.figure(1)
for i in range(len(num)):
    plt.plot(x[i], q[i], label=scheme + str(num[i]) + ' Divisions, ' + 'dt=' + str(np.format_float_scientific(dt[i], precision=3)))
plt.title(scheme + equation)
plt.grid()
plt.xlabel('X')
plt.ylabel('q')
plt.legend(loc='best')
plt.show()
# linear_advection(200, L, a, beta)

