#/usr/bin/env python3
"""
FRAGMENT-MNP: Micro and NanoPlastic FRAGMentation in the EnvironmeNT
"""
import numpy as np
from scipy.integrate import solve_ivp
import plotly.express as px
import plotly.graph_objects as go

# How many size classes and time steps
K = 7
T = 100

# Set up a particle size distribution - not actually needed in this example
d = np.logspace(-3, -9, K)

# Set an arbitrary number concentration of plastic in each size class
n_0 = np.full(K, 42.0)

# Fragmentation rate - set the 1e-3 for each size class. Note that this
# means fragmentation still happens from smallest size class (i.e. the mass
# doesn't balance)
k_frag = np.full(K, 0.01)

# Fragment size distribution matrix - assume fragmentation event results in even
# split between size classes of daughter particles
fsd = np.zeros((K,K))
for k in np.arange(K):
  fsd[k,:] = 1 / (K - k - 1) if (K - k) != 1 else 0
# Get the upper triangle of this matrix, which effectively sets f to zero for
# size classes larger (or equal to) than the current one
fsd = np.triu(fsd, k=1)

# Now we can set up the differential equation and solve using SciPy

# Define the function that satisfies n'(t) = f(t, n)
# i.e. the RHS of our differential eq
def f(t, n):
  # Get number of size classes and create empty result to be filled
  N = n.shape[0]
  dndt = np.empty(N)
  # Loop over the size classes and perform the calculation
  for k in np.arange(N):
    dndt[k] = - k_frag[k] * n[k] + np.sum(fsd[:,k] * k_frag * n)
  # Return the solution for all of the size classes
  return dndt

# Numerically solve this given the initial values for n, over T time steps
soln = solve_ivp(fun=f,
                 t_span=(0, T),
                 y0=n_0,
                 t_eval=np.arange(0, T))

# If k_frag != 0 for the smallest size class, then there will be a loss to the
# system, so keep track of that here
n_loss = np.sum(n_0) - np.sum(soln.y, axis=0)

# Finally we can plot the results

# Create the plot with the first size class and add the size classes
fig = go.Figure()
for i in range(0,K):
  fig.add_trace(go.Scatter(x=soln.t, y=soln.y[i], name=f'{d[i]} m'))
# Plot the loss with a different style
fig.add_trace(go.Scatter(x=soln.t, y=n_loss, name='Loss', line={'width': 3, 'dash': 'dash'}))
fig.update_layout(xaxis_title='Time', yaxis_title='Particle number concentration')
fig.show()