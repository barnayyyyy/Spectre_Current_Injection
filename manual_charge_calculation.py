# By: John Barney, Lucas Nichols

import numpy as np
from scipy.integrate import simps

# Define the constants
CURRENT_STOP = 0.000425
DAMPING_FACTOR_RISE = 1e-9  
DAMPING_FACTOR_FALL_STOP = 1.5e-9
TRAN_TIME = 4e-9  

# Define the double exponential current function
def I(t):
    return CURRENT_STOP * (np.exp(t/DAMPING_FACTOR_RISE) - np.exp(-t/DAMPING_FACTOR_FALL_STOP))

# Create an array of time values
t_values = np.linspace(0, TRAN_TIME, 1000)

# Calculate the current values at these time points
I_values = I(t_values)

# Perform the integration
Q = simps(I_values, t_values)

print(f'total charge transferred is {Q} Coulombs')