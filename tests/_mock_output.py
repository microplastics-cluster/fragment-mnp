"""
Mock output for testing purposes
"""
import numpy as np

from fragmentmnp.output import FMNPOutput

# Generate some arbitrary data
t = np.arange(0, 100)
psd = np.array([1e-6, 1e-3])
c = np.array([
    np.linspace(0, 42, 100),
    np.linspace(0, 10, 100)
])
n = c / (4.0 / 3.0) * np.pi * psd[:, np.newaxis] ** 3
c_diss = np.zeros((100,))
c_min = np.zeros((100,))
# Create output based on this. `soln` shouldn't be none,
# but we don't want to create a full ODE solution just for
# the test
mock_output = FMNPOutput(t, c, n, c_diss, c_min,
                         None, psd)
