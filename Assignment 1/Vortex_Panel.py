# Import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
# Import User Defined Functions

data_file = 'coordinates.txt'  # Import coordinates file
data = np.loadtxt(data_file, skiprows=1, dtype=float)

# Plot panels
x = data[:, 0]
y = data[:, 1]
plt.figure(0)
plt.plot(x, y)  # Plot geometry
circle = plt.Circle((0, 0), 1, fill=False, color='r')  # Create circle
plt.gca().add_patch(circle)  # Plot circle
plt.title('Panel Geometry')
plt.ylabel('y')
plt.xlabel('y')
# Add axis
plt.axvline(0, color='black', linewidth=1)
plt.axhline(0, color='black', linewidth=1)
plt.grid()
plt.show()


