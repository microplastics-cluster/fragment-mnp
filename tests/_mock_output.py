"""
Mock output for testing purposes
"""
import numpy as np
from fragmentmnp.output import FMNPOutput

t = np.arange(0.5,100.5,1.0)
c = np.vstack([np.full(100,2.0), np.full(100,3.0)])
n = np.vstack([np.full(100,1.0), np.full(100,1.5)])
c_diss = np.zeros(100)
c_min = np.zeros(100)
mock_output=FMNPOutput(t=t, c=c, n=n, c_diss=c_diss, c_min=c_min,
                       soln=None, psd=np.array([1e-6,2e-6]))