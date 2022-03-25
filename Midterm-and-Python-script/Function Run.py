from linear_advection import *
from linear__advection_diffusion import *
import matplotlib.pyplot as plt
# Define values
L = 1.  # Length
a = 1.  # Advection Value
beta = [1e-5, 1e-4, 1e-3, 1e-2]  # Diffusion Coefficient
dt_input = 1e-4
num = 200  # Number of divisions
courant = 0.25  # Courant Number

equation = 'Linear Advection-Diffusion Comparison, '  # Equation
scheme = 'Second Order Central, '  # Scheme for titles

# Initialize empty variables for plotting
dt = []
q = []
dx = []
x = []
# linear_advection_diffusion(200, L, a, beta, dt_input)
'''
Input which scheme and type you'd like to use here. Place a_ in front for linear advection, or d_ in front for linear 
advection diffusion. Can also call linear_advection(), or linear_advection_diffusion to see all 3 schemes. The 
available schemes are shown. 
1st order backwards - first_order_backward
2nd order backwards - second_order_backward
2nd order central - second_order_central
'''
for i in range(len(beta)):
    q_sol, x_sol, dt_sol, dx_sol = d_central_second_order(num, L, a, beta[i], dt_input)
    q.append(q_sol)
    x.append(x_sol)
    dt.append(dt_sol)
    dx.append(dx_sol)

plt.figure(1)
for i in range(len(x)):
    plt.plot(x[i], q[i], label=scheme + ' beta=' + str(np.format_float_scientific(beta[i], precision=3)))
plt.title(equation + str(num) + ' Divisions, dt=' + str(np.format_float_scientific(dt_input, precision=3)))
plt.grid()
plt.xlabel('X')
plt.ylabel('q')
plt.legend(loc='best')
plt.show()

# dt1 = 1e-4
# dt2 = 1e-6
# dt3 = 2.5e-4
# q1, x1, dt1, dx1 = a_backward_first_order(200, 1, 1, 1e-4, dt1)
# q2, x2, dt2, dx2 = a_central_second_order(200, 1, 1, 1e-4, dt2)
# q3, x3, dt3, dx3 = a_backward_second_order(200, 1, 1, 1e-4, dt3)
# dt.append(dt1)
# dt.append(dt2)
# dt.append(dt3)
# q.append(q1)
# q.append(q2)
# q.append(q3)
# x.append(x1)
# x.append(x2)
# x.append(x3)
# dx.append(dx1)
# dx.append(dx2)
# dx.append(dx3)
# linear_advection(200, L, a, beta)

