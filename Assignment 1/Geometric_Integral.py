# Imported packages
import numpy as np
import math as m
# User defined functions

np.seterr('raise')  # Raise floating point error when encountered


# Define function to calculate geometric integral I and J
def geometric_integral(control, boundary, phi, length):

    num_panel = len(control)    # Determine number of panels

    # Initialize array for integrals
    i_int = np.zeros([num_panel, num_panel])
    j_int = np.zeros([num_panel, num_panel])

    # Loop over panels to compute integrals
    for i in range(num_panel):  # Loop over i panels
        for j in range(num_panel):  # Loop over j panels
            if j != i:    # Calculate when i is not j
                A = -(control[i][0] - boundary[j][0]) * np.cos(phi[j]) - (control[i][1] - boundary[j][1]) * np.sin(phi[j])  # A Term
                B = (control[i][0] - boundary[j][0]) ** 2 + (control[i][1] - boundary[j][1]) ** 2   # B Term
                Cn = np.sin(phi[i] - phi[j])    # Normal C Term
                Dn = -(control[i][0] - boundary[j][0]) * np.sin(phi[i]) + (control[i][1] - boundary[j][1]) * np.cos(phi[i]) # Normal D Term
                Ct = -np.cos(phi[i] - phi[j])   # Tangential C Term
                Dt = (control[i][0] - boundary[j][0]) * np.cos(phi[i]) + (control[i][1] - boundary[j][1]) * np.sin(phi[i])  # Tangential D Term
                E = np.sqrt(B - A ** 2) # E Term
                # Set integral to 0 if E is 0 or complex or NAN or INF
                if E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E):
                    i_int[i, j] = 0
                    j_int[i, j] = 0
                else:
                    # Calculate I integral
                    i_term1 = 0.5 * Cn * np.log((length[j] + 2 * A * length[j] + B) / B)
                    i_term2 = ((Dn - A * Cn)/E) * (m.atan2((length[j] + A), E) - m.atan2(A, E))
                    i_int[i, j] = i_term1 + i_term2
                    # Calculate J integral
                    j_term1 = 0.5 * Ct * np.log((length[j] ** 2 + 2 * A * length[j] + B) / B)
                    j_term2 = ((Dn - A * Cn) / E) * (m.atan2((length[j] + A), E) - m.atan2(A, E))
                    j_int[i, j] = j_term1 + j_term2

            # Remove values that cause error
            if np.iscomplex(i_int[i, j]) or np.isnan(i_int[i, j]) or np.isinf(i_int[i, j]):
                i_int[i, j] = 0
            if np.iscomplex(j_int[i, j]) or np.isnan(j_int[i, j]) or np.isinf(j_int[i, j]):
                j_int[i, j] = 0

    return i_int, j_int

