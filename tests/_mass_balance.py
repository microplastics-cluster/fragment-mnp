"""
Test the model mass balance
"""
import numpy as np


def check_mass_balance(fmnp, output):
    """
    Check the mass balance of the model output is
    equal for each timestep and against initial conc data
    """
    # Calculate the total mass conc on each timestep
    mass = np.sum(output.c, axis=0) + output.c_diss + output.c_min
    # Prepend with the initial concs from data
    initial = np.sum(fmnp.initial_concs) + fmnp.initial_concs_diss
    mass = np.insert(mass, 0, initial)
    # Check all of the timsteps are equal within 5%
    return np.allclose(mass, mass[0], rtol=0.05)
